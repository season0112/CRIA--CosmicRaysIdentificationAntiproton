#!/usr/bin/env python
########################################################################
#                                                                      #
#              Pronton Charge Confusion with VGG16                     #
#                          03.05.2019                                  #
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

parser = argparse.ArgumentParser()
parser.add_argument('--inputlists', choices=["mva_variables", "mva_variables+rigiditywiouthit"], help='inputlists')
parser.add_argument('--validationsample', choices=["proton", "antiproton", "electron", "negative_events","validation_sample"], help='validation_data_choice')
parser.add_argument('--energyrange',  choices=["147-1000GeV", "147-1000GeV_v2", "16.6-147GeV_v2"], help='energy_range.')
parser.add_argument('--energybin',  choices=["16.6_38.9","16.6_80.5","80.5_93.0","93.0_108.0","108.0_125.0","125.0_147.0","147_175", "175_211", "211_259", "259_330", "330_525",], help='energybin')
parser.add_argument('--epochsnumber', type=int, help='epochsnumber for training. Default:22')
parser.add_argument('--batchsize', type=int, help='batchsize for training. Default:1000')
parser.add_argument('--learningrate', type=float, help='learningrate for training. Default:0.00005')
parser.add_argument('--optimizer',  choices=["Adam", "RMSprop", "SGD", ], help='optimizer. Default:Adam')
parser.add_argument('--binnumber',  type=int, help='binnumber for training. Default:80')

parser.add_argument('--loadweights', dest='loadweights', action='store_true',help='loadweights. Default:False')
parser.add_argument('--no-loadweights', dest='loadweights', action='store_false',help='loadweights. Default:False')
parser.set_defaults(loadweights=False)

parser.add_argument('--ecalBDTCut', dest='ecalBDTCut', action='store_true',help='ecalBDTCut. Default:False')
parser.add_argument('--no-ecalBDTCut', dest='ecalBDTCut', action='store_false',help='ecalBDTCut. Default:False')
parser.set_defaults(ecalBDTCut=False)

parser.add_argument('--trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_true',help='trdProtonHeliumCut. Default:False')
parser.add_argument('--no-trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_false',help='trdProtonHeliumCut. Default:False')
parser.set_defaults(trdProtonHeliumCut=False)

arguments = parser.parse_args()

########################## Free Parameters #############################
if (arguments.validationsample):
    Validationsample = arguments.validationsample
else:
    print("You need to choose a testdataset! Potential choises are: proton, antiproton, electron")
    os._exit(0)

if (arguments.energyrange):
    Energyrange = arguments.energyrange
else:
    print("You need to choose a energyrange! Potential choises are: 147-1000GeV, 147-1000GeV_v2, 16.6-147GeV_v2")
    os._exit(0)

if (arguments.energybin):
    Energybin = arguments.energybin
    if Energybin == "16.6_80.5":
        binnings = np.array([16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5])     
    elif Energybin == "16.6_38.9":
        binnings = np.array([16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9])
else:
    print("You need to choose a energybin! Potential choises are: 147_175, 175_211, 211_250, 250_330，330_525")
    os._exit(0)

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

if (arguments.epochsnumber):
    epochsnumber = arguments.epochsnumber
else:
    epochsnumber = 22

if (arguments.learningrate):
    learningrate = arguments.learningrate
else:
    learningrate = 0.00005

if (arguments.optimizer):
    opt = arguments.optimizer
else:
    opt = "Adam"

if opt == "RMSprop":
    opt = RMSprop(lr=learningrate, rho=0.9, epsilon=1e-6)
elif opt == "SDG":
    decay_rate = learning_rate / epochsnumber
    opt = SGD(lr=learningrate, momentum=0.8, decay=decay_rate, nesterov=False)
elif opt == "Adam":
    opt = Adam(lr=learningrate, beta_1=0.9, beta_2=0.999, epsilon=1e-08)

if (arguments.batchsize):
    batchsize = arguments.batchsize
else:
    batchsize = 1000

if (arguments.binnumber):
    binnumber = arguments.binnumber
else:
    binnumber = 80

Loadweights = arguments.loadweights

EcalBDTCut = arguments.ecalBDTCut

TrdProtonHeliumCut = arguments.trdProtonHeliumCut

########################################################################
EcalBDTCutvalue_select_proton = 0.0  ## -1 denote proton, 1 denote electron
EcalBDTCutvalue_select_electron = 0.0  ## -1 denote proton, 1 denote electron
TrdProtonHeliumCutValue = 0.3

########################## Properties Label ############################
MvAresult = 37
TrdLogLikelihoodRatioProtonHeliumTracker = 38
EcalBDT_EnergyD = 39

################### load data and features #############################
if Energybin == "16.6_80.5" or Energybin == "16.6_38.9":
    for binleft in range(binnings.shape[0]-1):    
        positiveprotonB1042_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
        CCprotonB1042_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
        CCprotonB1042_tem = CCprotonB1042_tem[0:int(CCprotonB1042_tem.shape[0]*0.7),:]
        positiveprotonB1042_tem = positiveprotonB1042_tem[0:CCprotonB1042_tem.shape[0],:]
        electronpositive_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy')
        electronnegative_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy')
        electronnegative_tem = electronnegative_tem[0:electronpositive_tem.shape[0],:]
        if binleft == 0:
            positiveprotonB1042 = positiveprotonB1042_tem
            CCprotonB1042 = CCprotonB1042_tem
            electronpositive = electronnegative_tem
            electronnegative = electronnegative_tem
        else:
            positiveprotonB1042 = np.row_stack((positiveprotonB1042_tem, positiveprotonB1042))
            CCprotonB1042 = np.row_stack((CCprotonB1042_tem, CCprotonB1042))
            electronpositive = np.row_stack((electronpositive_tem, electronpositive))
            electronnegative = np.row_stack((electronnegative_tem, electronnegative))
else:
    positiveprotonB1042 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
    CCprotonB1042 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
    CCprotonB1042 = CCprotonB1042[0:int(CCprotonB1042.shape[0]*0.7),:]
    positiveprotonB1042 = positiveprotonB1042[0:CCprotonB1042.shape[0],:]
    if Energybin == "211_259" or  Energybin == "259_330" or Energybin == "330_525":
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
    elif Energybin == "80.5_93.0" or "93.0_108.0" or "108.0_125.0" or "125.0_147.0" or "147_175":
        electronpositive = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy')
        electronnegative = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy')
        electronnegative = electronnegative[0:electronpositive.shape[0],:]
    else:
        print("Energybin error!")
##
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

########################## Target Sample ############################
trainingset = np.r_[correct,confused]
y_proton = [1]*correct.shape[0]
y_CCproton = [0]*confused.shape[0]
y_target = np.array(y_proton+y_CCproton)
y_target = np_utils.to_categorical(y_target, num_classes=2)

######################### Test events ##################################
if Validationsample == "proton":
    if Energybin == "16.6_80.5" or Energybin == "16.6_38.9":
        for binleft in range(binnings.shape[0]-1):
            testset1_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tem = np.row_stack((testset2_tem,np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")))
            testset1_tem = testset1_tem[-testset2_tem.shape[0]:,:]
            if binleft == 0:
                testset1 = testset1_tem
                testset2 = testset2_tem
            else:
                testset1 = np.row_stack((testset1_tem, testset1))
                testset2 = np.row_stack((testset2_tem, testset2))
    else:
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.row_stack((testset2,np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")))
        testset1 = testset1[-testset2.shape[0]:,:]

elif Validationsample == "antiproton":
    if Energybin == "16.6_80.5" or Energybin == "16.6_38.9":
        for binleft in range(binnings.shape[0]-1):
            testset1_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset1_tem = testset1_tem[-testset2_tem.shape[0]:,:]
            if binleft == 0:
                testset1 = testset1_tem
                testset2 = testset2_tem
            else:
                testset1 = np.row_stack((testset1_tem, testset1))
                testset2 = np.row_stack((testset2_tem, testset2))
    else:
            testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset1 = testset1[-testset2.shape[0]:,:]

elif Validationsample == "electron":
    if Energybin == "16.6_80.5" or Energybin == "16.6_38.9":
        for binleft in range(binnings.shape[0]-1):
            testset1_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset1_tem = testset1_tem[-testset2_tem.shape[0]:,:]  
            if binleft == 0:
                testset1 = testset1_tem
                testset2 = testset2_tem
            else:
                testset1 = np.row_stack((testset1_tem, testset1))
                testset2 = np.row_stack((testset2_tem, testset2))
    else:
        if Energybin == "211_259" or  Energybin == "259_330" or Energybin == "330_525":  
            testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset1 = testset1[-testset2.shape[0]:,:]
        elif Energybin == "175_211":
            testset1_1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset2_1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.2002000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset1_1 = testset1_1[np.where(testset1_1[:,-1]<-200)[0],:]
            testset2_1 = testset2_1[np.where(testset2_1[:,-1]>200)[0],:]
            testset1_1 = testset1_1[-testset2_1.shape[0]:,:]
            testset1_2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset2_2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset1_2 = testset1_2[np.where(testset1_2[:,-1]>-200)[0],:]
            testset2_2 = testset2_2[np.where(testset2_2[:,-1]<200)[0],:]
            testset1_2 = testset1_2[-testset2_2.shape[0]:,:]
            testset1 = np.row_stack((testset1_1,testset1_2))
            testset2 = np.row_stack((testset2_1,testset2_2))
        elif Energybin == "80.5_93.0" or "93.0_108.0" or "108.0_125.0" or "125.0_147.0" or "147_175":    
            testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1091_el.pl1.0_25200_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
            testset1 = testset1[-testset2.shape[0]:,:]
        else:
            print("Energy Range Error!")
            os._exit(0)  

elif Validationsample == "negative_events":
    if Energybin == "16.6_80.5" or Energybin == "16.6_38.9":
        for binleft in range(binnings.shape[0]-1):
            testset1_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset1_tem = testset1_tem[-testset2_tem.shape[0]:,:]
            if binleft == 0:
                testset1 = testset1_tem
                testset2 = testset2_tem
            else:
                testset1 = np.row_stack((testset1_tem, testset1))
                testset2 = np.row_stack((testset2_tem, testset2))        
    else:
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset1 = testset1[-testset2.shape[0]:,:]

elif Validationsample == "validation_sample":
    if Energybin == "16.6_80.5" or Energybin == "16.6_38.9":
        for binleft in range(binnings.shape[0]-1):
            testset1_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tep_tem = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/negative_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1")
            testset2_tem = np.row_stack((testset2_tem,testset2_tep_tem[-int(testset2_tep_tem.shape[0]*0.3):,:]))
#            testset2_tem = np.row_stack((testset2_tem,np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/positive/'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'GeV/pattern0/positive_'+str(binnings[binleft])+'_'+str(binnings[binleft+1])+'_pattern_0.npy',encoding="latin1"))) #### zero event is small bin, therefore not use them.
            testset1_tem = testset1_tem[-testset2_tem.shape[0]:,:]
            if binleft == 0:
                testset1 = testset1_tem
                testset2 = testset2_tem
            else:
                testset1 = np.row_stack((testset1_tem, testset1))
                testset2 = np.row_stack((testset2_tem, testset2))
    else:
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2_tep = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/'+Energyrange+'/transferdata/negative/'+Energybin+'GeV/pattern0/negative_'+Energybin+'_pattern_0.npy',encoding="latin1")
        testset2 = np.row_stack((testset2,testset2_tep[-int(testset2_tep.shape[0]*0.3):,:]))
        testset2 = np.row_stack((testset2,np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/'+Energyrange+'/transferdata/positive/'+Energybin+'GeV/pattern0/positive_'+Energybin+'_pattern_0.npy',encoding="latin1")))
        testset1 = testset1[-testset2.shape[0]:,:]
else:
        print("Please choose a kind of validation sample correctly!")
        os._exit(0)
##
if (EcalBDTCut):
    testset1 = testset1[np.where((testset1[:,EcalBDT_EnergyD]>-2) & (testset1[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
    testset2 = testset2[np.where((testset2[:,EcalBDT_EnergyD]>-2) & (testset2[:,EcalBDT_EnergyD] < EcalBDTCutvalue_select_proton))[0],:]
if (TrdProtonHeliumCut):
    testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
    testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

######################### Test Target Sample  ##########################
testset = np.r_[testset1,testset2]
testsetMvA = testset[:,MvAresult] ## get the CCProton MvA results
testset1_pre = [1]*len(testset1)
testset2_pre = [0]*len(testset2)
testset_pre = np.array(testset1_pre+testset2_pre)
testset_pre = np_utils.to_categorical(testset_pre, num_classes=2)


################  add a axis for convolutional 1D ######################
trainingset = np.expand_dims(trainingset, axis=2)
testset = np.expand_dims(testset, axis=2)

################  which features are used ###############################
sign_train = trainingset[:,-1]/np.abs(trainingset[:,-1])
sign_test = testset[:,-1]/np.abs(testset[:,-1])
trainingset = trainingset[:,inputlist]
testset = testset[:,inputlist]
for sign_ite in sign_independence:
    trainingset[:,sign_ite] = trainingset[:,sign_ite] * sign_train
    testset[:,sign_ite] = testset[:,sign_ite] * sign_test

#######################################################################
trainingset_for_predict = trainingset

####################### training data shuffling ########################
indices = np.arange(len(trainingset))
np.random.shuffle(indices)
trainingset = trainingset[indices]
y_target = y_target[indices]


##################### model and compile #################################
import CNN_models
model = CNN_models.VGG16(features)
if (Loadweights): 
    model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/VGG16_'+Energybin+'_'+arguments.inputlists+'.h5')
model.compile(optimizer=opt,
              loss='categorical_crossentropy',
              metrics=['accuracy'])
model.summary()

######################### Train the model, #############################
history = model.fit(trainingset, y_target, 
                    epochs=epochsnumber, 
                    batch_size=batchsize,
                    validation_data=(testset,testset_pre))

########################## Prediction ##################################
# on validation sample
y_pred = model.predict(testset)  
protonpred = 0
CCprotonpred = 0
for protonpredcount in range(testset1.shape[0]):
    if np.argmax(y_pred[protonpredcount])==1:
        protonpred+=1
for CCprotonpredcount in range(testset1.shape[0],y_pred.shape[0]):
    if np.argmax(y_pred[CCprotonpredcount])==0:
        CCprotonpred+=1
print("ChargeConrrect Prediction on validation sample:"+str(protonpred)+" out of "+str(testset1.shape[0]))
print("ChargeConfused Prediction on validation sample:"+str(CCprotonpred)+" out of "+str(testset2.shape[0]))
# on training sample
y_pred_on_train = model.predict(trainingset_for_predict)
protonpred_on_train = 0
CCprotonpred_on_train = 0
for protonpredcount_on_train in range(correct.shape[0]):
    if np.argmax(y_pred_on_train[protonpredcount_on_train])==1:
        protonpred_on_train+=1
for CCprotonpredcount_on_train in range(correct.shape[0],y_pred_on_train.shape[0]):
    if np.argmax(y_pred_on_train[CCprotonpredcount_on_train])==0:
        CCprotonpred_on_train+=1
print("ChargeConrrect Prediction on training sample:"+str(protonpred_on_train)+" out of "+str(correct.shape[0]))
print("ChargeConfused Prediction on training sample:"+str(CCprotonpred_on_train)+" out of "+str(confused.shape[0]))

######################### model weights ################################
model.save_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/VGG16_'+Energybin+'_'+arguments.inputlists+'.h5')

################## Plot Loss and Accuracy Prediction ###################
plt.figure(figsize=(18,18))
plt.subplot(221)
loss = history.history['loss']
val_loss = history.history['val_loss']
epochs = range(1, len(loss) + 1)
plt.plot(epochs, loss, 'bo', label='Training loss')
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and validation loss',fontsize=26)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Epochs',fontsize=22)
plt.ylabel('Loss',fontsize=22)
plt.legend(loc='best',fontsize=20)
plt.subplot(222)
acc = history.history['acc']
val_acc = history.history['val_acc']
plt.plot(epochs, acc, 'bo', label='Training acc')
plt.plot(epochs, val_acc, 'b', label='Validation acc')
plt.title('Training and validation accuracy',fontsize=26)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Epochs',fontsize=22)
plt.ylabel('Accuracy',fontsize=22)
plt.legend(loc='best',fontsize=20)
plt.subplot(223)
plt.hist(y_pred[0:testset1.shape[0],1],bins=binnumber,range=(0,1),log=True ) 
plt.title("Proton",fontsize=26) 
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Estimator $_{CC}$',fontsize=22)
plt.ylabel('Count',fontsize=22)
plt.subplot(224)
plt.hist(y_pred[testset1.shape[0]:y_pred.shape[0],1],bins=binnumber,range=(0,1),log=True ) 
plt.title("CCProton",fontsize=26) 
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Estimator $_{CC}$',fontsize=22)
plt.ylabel('Count',fontsize=22)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/plot/ML_'+Energybin+'.png')

######################## Plot MvA results ##############################
protonMvApred = 0
CCprotonMvApred = 0
for protonMvApredcount in range(testset1.shape[0]):
    if testsetMvA[protonMvApredcount]>0.5:
        protonMvApred+=1
for CCprotonMvApredcount in range(testset1.shape[0],y_pred.shape[0]):
    if testsetMvA[CCprotonMvApredcount]<0.5:
        CCprotonMvApred+=1
plt.figure(figsize=(18,18))
plt.hist(testsetMvA[0:testset1.shape[0]],binnumber,range=(0,1),log=True, alpha=0.5,label='Proton',facecolor='blue',edgecolor='black'  )
plt.hist(testsetMvA[testset1.shape[0]:y_pred.shape[0]],binnumber,range=(0,1),log=True, alpha=0.5,label='CCProton',facecolor='green',edgecolor='black'  ) 
plt.xlabel('Estimator $_{CC}$',fontsize=30)
plt.ylabel('Count',fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(loc='upper center',fontsize=30)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/plot/CCMVA_'+Energybin+'.png')

######################### fraction plot ################################
#
protonaccumulatecount=0
fractionoprotonset=np.array([])
for protoniterate in range(binnumber):
    protonaccumulatecount=protonaccumulatecount+np.histogram(y_pred[0:testset1.shape[0],1],binnumber,range=(0,1))[0][binnumber-protoniterate-1]
    fractionoproton=protonaccumulatecount/(testset1.shape[0])
    fractionoprotonset=np.append(fractionoprotonset,fractionoproton)
CCprotonaccumulatecount=0
rejectionCCprotonset=np.array([])
for CCprotoniterate in range(binnumber):
    CCprotonaccumulatecount=CCprotonaccumulatecount+np.histogram(y_pred[testset1.shape[0]:y_pred.shape[0],1],binnumber,range=(0,1))[0][binnumber-CCprotoniterate-1]
    rejectionCCproton=1.0/(CCprotonaccumulatecount/(testset2.shape[0]))
    rejectionCCprotonset=np.append(rejectionCCprotonset,rejectionCCproton)
#
protonaccumulatecount_on_train = 0
fractionoprotonset_on_train = np.array([])
for protoniterate_on_train in range(binnumber):
    protonaccumulatecount_on_train = protonaccumulatecount_on_train +np.histogram(y_pred_on_train[0:correct.shape[0],1],binnumber,range=(0,1))[0][binnumber-protoniterate_on_train-1]
    fractionoproton_on_train = protonaccumulatecount_on_train/(correct.shape[0])
    fractionoprotonset_on_train = np.append(fractionoprotonset_on_train, fractionoproton_on_train)
CCprotonaccumulatecount_on_train = 0
rejectionCCprotonset_on_train = np.array([])
for CCprotoniterate_on_train in range(binnumber):
    CCprotonaccumulatecount_on_train = CCprotonaccumulatecount_on_train + np.histogram(y_pred_on_train[correct.shape[0]:y_pred_on_train.shape[0],1],binnumber,range=(0,1))[0][binnumber-CCprotoniterate_on_train-1]
    rejectionCCproton_on_train = 1.0/(CCprotonaccumulatecount_on_train/(confused.shape[0]))
    rejectionCCprotonset_on_train = np.append(rejectionCCprotonset_on_train,rejectionCCproton_on_train)
#
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
plt.plot(fractionoprotonsetMvA, rejectionCCprotonsetMvA, 'b-', lw=3,label='MVA')
plt.plot(fractionoprotonset_on_train, rejectionCCprotonset_on_train, 'r-',lw=3, label='Neurual Network on training sample')
plt.xlabel('Efficiency of Proton',fontsize=22)
plt.ylabel('Rejection for CCproton',fontsize=22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.legend( loc='best',fontsize=20)
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/plot/rejection_'+Energybin+'.png')

