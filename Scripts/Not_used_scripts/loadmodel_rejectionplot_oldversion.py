#!/usr/bin/env python
########################################################################
#                                                                      #
#                   Pronton Charge Confusion with VGG16                #
#                                 03.05.2019                           #
#                                                                      # 
########################################################################


########################## packages ####################################
from __future__ import division
import numpy as np
import os
import math
import json
import collections
import matplotlib.pyplot as plt 
import argparse
from mpl_toolkits.mplot3d import Axes3D
from keras.layers import Activation, Dropout, Flatten, Dense  
from keras.layers import Conv1D, GlobalAveragePooling1D, MaxPooling1D,ZeroPadding1D,Convolution1D  
from keras.models import Sequential, Model, model_from_json   
from keras.optimizers import SGD, RMSprop, Adam
from keras.preprocessing.image import ImageDataGenerator   
from keras.utils import np_utils
from keras import initializers
from keras.callbacks import LearningRateScheduler
plt.switch_backend('agg')
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist

parser = argparse.ArgumentParser()
parser.add_argument('--inputlists', choices=["mva_variables", "mva_variables+rigiditywiouthit"], help='inputlists')
parser.add_argument('--testsample', choices=["proton", "antiproton", "electron", "negative_events","validation_sample","B1119_ele","B1200proton"], help='validation_data_choice')
parser.add_argument('--energybin',  choices=["38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93.0","93.0_108.0","108.0_125.0","125.0_147.0","147_175", "175_211", "211_259", "259_330", "330_525",], help='energybin')
parser.add_argument('--energyrange',  choices=["147-1000GeV", "147-1000GeV_v2", "38.9-147GeV_v2"], help='energy_range.')
parser.add_argument('--binnumber',  type=int, help='binnumber for training. Default:80')
parser.add_argument('--ecalBDTCut', dest='ecalBDTCut', action='store_true',help='ecalBDTCut. Default:False')
parser.add_argument('--no-ecalBDTCut', dest='ecalBDTCut', action='store_false',help='ecalBDTCut. Default:False')
parser.set_defaults(ecalBDTCut=False)
parser.add_argument('--trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_true',help='trdProtonHeliumCut. Default:False')
parser.add_argument('--no-trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_false',help='trdProtonHeliumCut. Default:False')
parser.set_defaults(trdProtonHeliumCut=False)
arguments = parser.parse_args()

########################## Free Parameters #############################
if (arguments.inputlists):
    if arguments.inputlists == "mva_variables":
        sign_independence = np.arange(4,7)
        inputlist = np.arange(0,4)
        inputlist = np.append(inputlist,np.arange(7,19))
        features = 16
    if arguments.inputlists == "mva_variables+rigiditywiouthit":
        s1 = np.arange(4,7)
        s2 = np.arange(16,25)
        sign_independence = np.append(s1,s2)
        inputlist = np.arange(0,4)
        inputlist = np.append(inputlist,np.arange(7,28))
        features = 25
else:
    print("You need to choose a inputlist! Potential choises are: mva_variables, mva_variables+rigiditywiouthit")
    os._exit(0)

if (arguments.testsample):
    Testsample = arguments.testsample
else:
    print("You need to choose a testdataset! Potential choises are: proton, antiproton, electron")
    os._exit(0)

if (arguments.energyrange):
    Energyrange = arguments.energyrange
else:
    print("You need to choose a energyrange! Potential choises are: 147-1000GeV, 147-1000GeV_v2, 38.9-147GeV_v2")
    os._exit(0)

if (arguments.energybin):
    Energybin = arguments.energybin
else:
    print("You need to choose a energybin! Potential choises are: 147_175, 175_211, 211_250, 250_330，330_525")
    os._exit(0)

if (arguments.binnumber):
    binnumber = arguments.binnumber
else:
    binnumber = 80

EcalBDTCut = arguments.ecalBDTCut

TrdProtonHeliumCut = arguments.trdProtonHeliumCut


########################################################################
EcalBDTCutvalue_select_proton = 0.0  ## -1 denote proton, 1 denote electron
EcalBDTCutvalue_select_electron = 0.0  ## -1 denote proton, 1 denote electron
TrdProtonHeliumCutValue = 0.3

########################## Properties Label ############################
ElectronCCMVABDT = 36
MvAresult = 37
TrdLogLikelihoodRatioProtonHeliumTracker = 38
EcalBDT_EnergyD = 39

######################### Training Sample #############################
positiveprotonB1042 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
CCprotonB1042 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
CCprotonB1042 = CCprotonB1042[0:int(CCprotonB1042.shape[0]*0.7),:]
positiveprotonB1042 = positiveprotonB1042[0:CCprotonB1042.shape[0],:]
if Energybin == "211_259" or Energybin == "259_330" or Energybin == "330_525":
    electronpositive = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy')
    electronnegative = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy')
    electronnegative = electronnegative[0:electronpositive.shape[0],:]
elif Energybin == "175_211":
    electronpositive1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy')
    electronnegative1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy')
    electronpositive1 = electronpositive1[np.where(electronpositive1[:,-1]>200)[0],:]
    electronnegative1 = electronnegative1[np.where(electronnegative1[:,-1]<-200)[0],:]
    electronnegative1 = electronnegative1[0:electronpositive1.shape[0],:]
    electronpositive2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy')
    electronnegative2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy')
    electronpositive2 = electronpositive2[np.where(electronpositive2[:,-1]<200)[0],:]
    electronnegative2 = electronnegative2[np.where(electronnegative2[:,-1]>-200)[0],:]
    electronnegative2 = electronnegative2[0:electronpositive2.shape[0],:]
    electronpositive = np.row_stack((electronpositive1,electronpositive2))
    electronnegative = np.row_stack((electronnegative1,electronnegative2))
else:
    electronpositive = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy')
    electronnegative = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy')
    electronnegative = electronnegative[0:electronpositive.shape[0],:]
if (EcalBDTCut):
    positiveprotonB1042=positiveprotonB1042[np.where((positiveprotonB1042[:,EcalBDT_EnergyD]>-2) & (positiveprotonB1042[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
    CCprotonB1042=CCprotonB1042[np.where((CCprotonB1042[:,EcalBDT_EnergyD]>-2) & (CCprotonB1042[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
    electronpositive = electronpositive[np.where((electronpositive[:,EcalBDT_EnergyD]>-2) & (electronpositive[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
    electronnegative = electronnegative[np.where((electronnegative[:,EcalBDT_EnergyD]>-2) & (electronnegative[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
if (TrdProtonHeliumCut):
    positiveprotonB1042 = positiveprotonB1042[np.where(positiveprotonB1042[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
    CCprotonB1042 = CCprotonB1042[np.where(CCprotonB1042[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
    electronpositive = electronpositive[np.where(electronpositive[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
    electronnegative = electronnegative[np.where(electronnegative[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
correct = np.row_stack((positiveprotonB1042,electronnegative))
confused = np.row_stack((CCprotonB1042,electronpositive))

######################### Test Sample #################################
if Testsample == "proton":
    testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2 = np.row_stack((testset2,np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")))
    testset1 = testset1[-testset2.shape[0]:,:]
    if (EcalBDTCut):
        testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
        testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
    if (TrdProtonHeliumCut):
        testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

elif Testsample == "B1200proton":
    testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1200_pr.pl1.l1.054000.4_00_7.7_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1200_pr.pl1.l1.054000.4_00_7.7_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset1 = testset1[-testset2.shape[0]:,:]
    if (EcalBDTCut):
        testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
        testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
    if (TrdProtonHeliumCut):
        testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

elif Testsample == "antiproton":
    testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset1 = testset1[-testset2.shape[0]:,:]
    if (EcalBDTCut):
        testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
        testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
    if (TrdProtonHeliumCut):
        testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

elif Testsample == "electron":
    if Energybin == "211_259" or Energybin == "259_330" or Energybin == "330_525":
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    elif Energybin == "175_211":
        testset1_below = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2_below = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset1_above = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2_above = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2_above = testset2_above[np.where(testset2_above[:,-1]>200)[0],:]
        testset1_above = testset1_above[np.where(testset1_above[:,-1]<-200)[0],:]
        testset2_below = testset2_below[np.where(testset2_below[:,-1]<200)[0],:]
        testset1_below = testset1_below[np.where(testset1_below[:,-1]>-200)[0],:]        
        testset1 = np.row_stack((testset1_above,testset1_below))
        testset2 = np.row_stack((testset2_above,testset2_below))
    else:
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset1 = testset1[-testset2.shape[0]:,:]
    if (EcalBDTCut):
        testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
        testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
    if (TrdProtonHeliumCut):
        testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
    print('testset1'+str(testset1.shape))
    print('testset2'+str(testset2.shape))
    testset1 = testset1[np.where(testset1[:,ElectronCCMVABDT] > -2.0)[0],:]
    testset2 = testset2[np.where(testset2[:,ElectronCCMVABDT] > -2.0)[0],:]
    print('testset1'+str(testset1.shape))
    print('testset2'+str(testset2.shape))
    
elif Testsample == "B1119_ele":
    if Energybin == "211_259" or Energybin == "259_330" or Energybin == "330_525":
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.2004000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.2004000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    elif Energybin == "175_211":
        testset1_below = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2_below = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset1_above = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.2004000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2_above = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.2004000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2_above = testset2_above[np.where(testset2_above[:,-1]>200)[0],:]
        testset1_above = testset1_above[np.where(testset1_above[:,-1]<-200)[0],:]
        testset2_below = testset2_below[np.where(testset2_below[:,-1]<200)[0],:]
        testset1_below = testset1_below[np.where(testset1_below[:,-1]>-200)[0],:]
        testset1 = np.row_stack((testset1_above,testset1_below))
        testset2 = np.row_stack((testset2_above,testset2_below))
    else:
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1119_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset1 = testset1[-testset2.shape[0]:,:]
    if (EcalBDTCut):
        testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
        testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
    if (TrdProtonHeliumCut):
        testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
    print('testset1'+str(testset1.shape))
    print('testset2'+str(testset2.shape))
    testset1 = testset1[np.where(testset1[:,ElectronCCMVABDT] > -2.0)[0],:]
    testset2 = testset2[np.where(testset2[:,ElectronCCMVABDT] > -2.0)[0],:]
    print('testset1'+str(testset1.shape))
    print('testset2'+str(testset2.shape))

elif Testsample == "negative_events":
    testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset1 = testset1[-testset2.shape[0]:,:]
    if (EcalBDTCut):
        testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
        testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] > EcalBDTCutvalue_select_electron))[0],:]
    if (TrdProtonHeliumCut):
        testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

elif Testsample == "validation_sample":
    testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2_tep = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    testset2 = np.row_stack((testset2,testset2_tep[-int(testset2_tep.shape[0]*0.3):,:]))
    testset2 = np.row_stack((testset2,np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")))
    testset1 = testset1[-testset2.shape[0]:,:]
    if (EcalBDTCut):
        testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
        testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
    if (TrdProtonHeliumCut):
        testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

else:
        print("Please choose a kind of validation sample correctly!")
        os._exit(0)


######################### Test Target Sample  ##########################
trainingset = np.r_[correct,confused]
testset = np.r_[testset1,testset2]
if Testsample == "electron" or "B1119_ele" : 
    testsetMvA = testset[:,ElectronCCMVABDT] ## get the Electron Charge Confusion BDT results
else:
    testsetMvA = testset[:,MvAresult] ## get the Proton Charge Confusion MvA results
trainingset = np.expand_dims(trainingset, axis=2)
testset = np.expand_dims(testset, axis=2)
sign_train = trainingset[:,-1]/np.abs(trainingset[:,-1])
sign_test = testset[:,-1]/np.abs(testset[:,-1])
trainingset = trainingset[:,inputlist]
testset = testset[:,inputlist]
for sign_ite in sign_independence:
    testset[:,sign_ite] = testset[:,sign_ite] * sign_test
    trainingset[:,sign_ite] = trainingset[:,sign_ite] * sign_train

##############  model prediction   #################
import CNN_models
model = CNN_models.VGG16(features)
if Energybin == "80.5_93.0" or Energybin == "93.0_108.0" or Energybin == "108.0_125.0" or Energybin == "125.0_147.0" or Energybin == "147_175" or Energybin == "175_211" or Energybin == "211_259" or Energybin == "259_330" or Energybin == "330_525":
    model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/VGG16_'+Energybin+'_'+arguments.inputlists+'.h5')
else:
    model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/VGG16_'+'38.9_80.5'+'_'+arguments.inputlists+'.h5')
y_pred = model.predict(testset)
y_pred_train = model.predict(trainingset)


'''
########################## Prediction ##################################
protonpred = 0
CCprotonpred = 0
for protonpredcount in range(testset1.shape[0]):
    if np.argmax(y_pred[protonpredcount])==1:
        protonpred+=1
for CCprotonpredcount in range(testset1.shape[0],y_pred.shape[0]):
    if np.argmax(y_pred[CCprotonpredcount])==0:
        CCprotonpred+=1
print("ChargeConrrect Prediction on test sample:"+str(protonpred)+" out of "+str(testset1.shape[0]))
print("ChargeConfused Prediction on test sample:"+str(CCprotonpred)+" out of "+str(testset2.shape[0]))
'''
######################### fraction plot ################################
# Neural Network application on training sample
correct_accumulate_count = 0
fractionoprotonset_train = np.array([])
for correct_iterate in range(binnumber):
    correct_accumulate_count = correct_accumulate_count+np.histogram(y_pred_train[0:correct.shape[0],1],binnumber,range=(0,1))[0][binnumber-correct_iterate-1]
    fractionoproton_train = correct_accumulate_count/(correct.shape[0])
    fractionoprotonset_train = np.append(fractionoprotonset_train,fractionoproton_train)
CC_accumulate_count = 0
rejectionCCprotonset_train = np.array([])
for CCiterate in range(binnumber):
    CC_accumulate_count = CC_accumulate_count+np.histogram(y_pred_train[correct.shape[0]:y_pred_train.shape[0],1],binnumber,range=(0,1))[0][binnumber-CCiterate-1]
    rejectionCCproton_train = 1.0/(CC_accumulate_count/(confused.shape[0]))
    rejectionCCprotonset_train = np.append(rejectionCCprotonset_train,rejectionCCproton_train)
# Neural Network application on test sample  
protonaccumulatecount=0
fractionoprotonset=np.array([])
for protoniterate in range(binnumber):
    protonaccumulatecount = protonaccumulatecount+np.histogram(y_pred[0:testset1.shape[0],1],binnumber,range=(0,1))[0][binnumber-protoniterate-1]
    fractionoproton = protonaccumulatecount/(testset1.shape[0])
    fractionoprotonset = np.append(fractionoprotonset,fractionoproton)
CCprotonaccumulatecount = 0
rejectionCCprotonset = np.array([])
for CCprotoniterate in range(binnumber):
    CCprotonaccumulatecount = CCprotonaccumulatecount+np.histogram(y_pred[testset1.shape[0]:y_pred.shape[0],1],binnumber,range=(0,1))[0][binnumber-CCprotoniterate-1]
    rejectionCCproton = 1.0/(CCprotonaccumulatecount/(testset2.shape[0]))
    if CCprotonaccumulatecount == 0:
        print('CCprotonaccumulatecount='+str(CCprotonaccumulatecount))
        print('rejectionCCproton='+str(rejectionCCproton))
    rejectionCCprotonset = np.append(rejectionCCprotonset,rejectionCCproton)
print(fractionoprotonset)
print(rejectionCCprotonset)
# TMVA application on test sample (Proton from Andi, Electron from Niko)
if Testsample == "electron" or Testsample == "B1119_ele":
    maxvalue = max(max(testsetMvA[0:testset1.shape[0]]),max(testsetMvA[testset1.shape[0]:y_pred.shape[0]]))
    minvalue = min(min(testsetMvA[0:testset1.shape[0]]),min(testsetMvA[testset1.shape[0]:y_pred.shape[0]]))
    electronaccumulatecountMvA=0
    fractionoelectronsetMvA=np.array([])
    for electroniterateMvA in range(binnumber):
        electronaccumulatecountMvA=electronaccumulatecountMvA+np.histogram(testsetMvA[0:testset1.shape[0]],binnumber,range=(minvalue,maxvalue))[0][binnumber-electroniterateMvA-1]
        fractionoelectronMvA=electronaccumulatecountMvA/(testset1.shape[0])
        fractionoelectronsetMvA=np.append(fractionoelectronsetMvA,fractionoelectronMvA)
    CCelectronaccumulatecountMvA=0
    rejectionCCelectronsetMvA=np.array([])
    for CCelectroniterateMvA in range(binnumber):
        CCelectronaccumulatecountMvA=CCelectronaccumulatecountMvA+np.histogram(testsetMvA[testset1.shape[0]:y_pred.shape[0]],binnumber,range=(minvalue,maxvalue))[0][binnumber-CCelectroniterateMvA-1]
        rejectionCCelectronMvA=1.0/(CCelectronaccumulatecountMvA/(testset2.shape[0]))
        rejectionCCelectronsetMvA=np.append(rejectionCCelectronsetMvA,rejectionCCelectronMvA)
else:
    protonaccumulatecountMvA=0
    fractionoprotonsetMvA=np.array([])
    for protoniterateMvA in range(binnumber):
        protonaccumulatecountMvA=protonaccumulatecountMvA+np.histogram(testsetMvA[0:testset1.shape[0]],binnumber,range=(0,1))[0][binnumber-protoniterateMvA-1]
        fractionoprotonMvA=protonaccumulatecountMvA/(testset1.shape[0])
        fractionoprotonsetMvA=np.append(fractionoprotonsetMvA,fractionoprotonMvA)
    CCprotonaccumulatecountMvA=0
    rejectionCCprotonsetMvA=np.array([])
    for CCprotoniterateMvA in range(binnumber):
        CCprotonaccumulatecountMvA=CCprotonaccumulatecountMvA+np.histogram(testsetMvA[testset1.shape[0]:y_pred.shape[0]],binnumber,range=(0,1))[0][binnumber-CCprotoniterateMvA-1]
        rejectionCCprotonMvA=1.0/(CCprotonaccumulatecountMvA/(testset2.shape[0]))
        rejectionCCprotonsetMvA=np.append(rejectionCCprotonsetMvA,rejectionCCprotonMvA)
#
plt.figure(figsize=(18,9))
#plt.ylim(1, rejectionCCprotonsetMvA[0]+2)
plt.plot(fractionoprotonset, rejectionCCprotonset, 'g-',lw=3, label='Neurual Network on validation sample')
plt.plot(fractionoprotonset_train, rejectionCCprotonset_train, 'r-',lw=3, label='Neurual Network on training sample')
if Testsample == "electron" or Testsample == "B1119_ele":
    plt.plot(fractionoelectronsetMvA, rejectionCCelectronsetMvA, 'b-', lw=3,label='Charge Confusion Electron BDT')
    plt.xlabel('Efficiency of Electron',fontsize=22)
    plt.ylabel('Rejection for CCElectron',fontsize=22)
else:
    plt.plot(fractionoprotonsetMvA, rejectionCCprotonsetMvA, 'b-', lw=3,label='Charge Confusion Proton TMVA')
    plt.xlabel('Efficiency of Proton',fontsize=22)
    plt.ylabel('Rejection for CCproton',fontsize=22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.legend( loc='best',fontsize=20)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/rejection_'+Energybin+'_'+arguments.inputlists+'.png')




########################## NN Prediction ##################################
#
plt.figure(figsize=(18,18))
plt.hist(y_pred[0:testset1.shape[0],1],bins=binnumber,range=(0,1),density=True,log=True,alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.xlabel('Estimator $_{CC}$',fontsize=35)
plt.ylabel('Count',fontsize=35)
plt.legend(loc='upper center',fontsize=35)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/ML_'+Energybin+'ChargeCorrect.png')
#
plt.figure(figsize=(18,18))
plt.hist(y_pred[testset1.shape[0]:y_pred.shape[0],1],bins=binnumber,range=(0,1),density=True,log=True,alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.xlabel('Estimator $_{CC}$',fontsize=35)
plt.ylabel('Count',fontsize=35)
plt.legend(loc='upper center',fontsize=35)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/ML_'+Energybin+'ChargeConfused.png')

plt.figure(figsize=(18,18))
plt.hist(y_pred[0:testset1.shape[0],1],bins=binnumber,range=(0,1),density=True,log=True,alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(y_pred[testset1.shape[0]:y_pred.shape[0],1],bins=binnumber,range=(0,1),density=True,log=True,alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.xlabel('Estimator $_{CC}$',fontsize=35)
plt.ylabel('Count',fontsize=35)
plt.legend(loc='upper center',fontsize=35)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/ML_'+Energybin+'.png')

########################  MvA Prediction ##############################
plt.figure(figsize=(18,18))
plt.hist(testsetMvA[0:testset1.shape[0]],binnumber,range=(-1,1),density=True,log=True, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black'  )
plt.hist(testsetMvA[testset1.shape[0]:y_pred.shape[0]],binnumber,range=(-1,1),density=True,log=True, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black'  )
plt.xlabel('Estimator $_{CC}$',fontsize=30)
plt.ylabel('Count',fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(loc='upper center',fontsize=30)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/CCMVA_'+Energybin+'.png')



###################### ROOT PLOT #################################################################################################################
ML_p = TH1D("ML_p","", binnumber,0,1)
ML_n = TH1D("ML_n","", binnumber,0,1)
CCMVA_p = TH1D("CCMVA_p","", binnumber,-1,1)
CCMVA_n = TH1D("CCMVA_n","", binnumber,-1,1)

ML_p.SetFillColor(6)
ML_p.SetFillStyle(3004)
ML_p.SetLineColor(6)
ML_n.SetFillColor(4)
ML_n.SetFillStyle(3005)
ML_n.SetLineColor(4)
CCMVA_p.SetFillColor(6)
CCMVA_p.SetFillStyle(3004)
CCMVA_p.SetLineColor(6)
CCMVA_n.SetFillColor(4)
CCMVA_n.SetFillStyle(3005)
CCMVA_n.SetLineColor(4)


fill_hist(ML_p, y_pred[0:testset1.shape[0],1])
fill_hist(CCMVA_p, testsetMvA[0:testset1.shape[0]])
fill_hist(ML_n, y_pred[testset1.shape[0]:y_pred.shape[0],1])
fill_hist(CCMVA_n, testsetMvA[testset1.shape[0]:y_pred.shape[0]])


scale = 150/ML_p.Integral()
ML_p.Scale(scale)
scale = 150/ML_n.Integral()
ML_n.Scale(scale)
scale = 150/CCMVA_p.Integral()
CCMVA_p.Scale(scale)
scale = 150/CCMVA_n.Integral()
CCMVA_n.Scale(scale)


c1 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
ML_p.Draw("HIST")
ML_n.Draw('HIST same')
leg =TLegend(.4,.7,.6,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
ML_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c1.Update()
c1.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/ML_seperation_"+Energybin+".pdf")

c2 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
CCMVA_p.Draw("HIST")
CCMVA_n.Draw('HIST same')
leg =TLegend(.7,.7,.9,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
CCMVA_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c2.Update()
c2.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/CCMVA_seperation_"+Energybin+".pdf")

c3 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
ML_p.Draw("HIST")
ML_n.Draw('HIST same')
leg =TLegend(.4,.7,.6,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
ML_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c3.SetLogy()
c3.Update()
c3.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/ML_seperation_log_"+Energybin+".pdf")

c4 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
CCMVA_p.Draw("HIST")
CCMVA_n.Draw('HIST same')
leg =TLegend(.7,.7,.9,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
CCMVA_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c4.SetLogy()
c4.Update()
c4.SaveAs("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/plot/CCMVA_seperation_log_"+Energybin+".pdf")



