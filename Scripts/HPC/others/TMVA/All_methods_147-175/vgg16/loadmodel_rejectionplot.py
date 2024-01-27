#!/usr/bin/env python
########################################################################
#                                                                      #
#       Toy Model for Pronton Charge Confusion with Dense Layer        #
#                      version 1.0   07.03.2018                        #
#                                                                      # 
########################################################################

########################## packages ####################################
from __future__ import division
import numpy as np
import math
import json
import collections
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LogNorm
from keras.layers import Activation, Dropout, Flatten, Dense  
from keras.layers import Dense, Conv2D, MaxPooling2D  
from keras.models import Sequential, Model, model_from_json   
from keras.optimizers import SGD, RMSprop, Adam
from keras.preprocessing.image import ImageDataGenerator   
from keras.utils import np_utils
from keras import initializers
#from plt import ImageFont
#from plt import Image
#from plt import ImageDraw
plt.switch_backend('agg')
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist


#########################################################
MLPBNN_signal = np.load('../mva_MLPBNN_signal_array.npy')
MLPBNN_background = np.load('../mva_MLPBNN_background_array.npy')

SVM_signal = np.load('../mva_SVM_signal_array.npy')
SVM_background = np.load('../mva_SVM_background_array.npy')

LikelihoodKDE_signal = np.load('../mva_LikelihoodKDE_signal_array.npy')
LikelihoodKDE_background = np.load('../mva_LikelihoodKDE_background_array.npy')

Fisher_signal = np.load('../mva_Fisher_signal_array.npy')
Fisher_background = np.load('../mva_Fisher_background_array.npy')

Likelihood_signal = np.load('../mva_Likelihood_signal_array.npy')
Likelihood_background = np.load('../mva_Likelihood_background_array.npy')

BDTG_signal = np.load('../mva_BDTG_signal_array.npy')
BDTG_background = np.load('../mva_BDTG_background_array.npy')

# rerange
Fisher_signal[np.where(Fisher_signal<-1)[0]] = -1
Fisher_background[np.where(Fisher_background<-1)[0]] = -1
Fisher_signal[np.where(Fisher_signal>1)[0]] = 1
Fisher_background[np.where(Fisher_background>1)[0]] = 1

Likelihood_signal[np.where(Likelihood_signal<-1)[0]] = -1
Likelihood_background[np.where(Likelihood_background<-1)[0]] = -1
Likelihood_signal[np.where(Likelihood_signal>1)[0]] = 1
Likelihood_background[np.where(Likelihood_background>1)[0]] = 1

########################## Free Parameters #############################
binnumber=150
binnumber_KDE=20000
binnumber_SVM=2000
binnumber_NN=500
EcalBDTCut = False
EcalBDTCutvalue_select_proton = 0.0  ## -1 denote proton, 1 denote electron
EcalBDTCutvalue_select_electron = 0.0  ## -1 denote proton, 1 denote electron
TrdProtonHeliumCut = False
TrdProtonHeliumCutValue = 0.3
Testsample = "proton"

########################## Properties Label ############################
MvAresult=16
TrdLogLikelihoodRatioProtonHeliumTracker=17

######################### Test events ##################################
if Testsample is "proton":
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/147-1000GeV/transferdata/positive/147_175GeV/pattern0/positive_147_175_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/147-1000GeV/transferdata/negative/147_175GeV/pattern0/negative_147_175_pattern_0.npy',encoding="latin1")
        testset2 = np.row_stack((testset2,np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/147-1000GeV/transferdata/negative/147_175GeV/pattern0/negative_147_175_pattern_0.npy',encoding="latin1")))
        testset1 = testset1[-testset2.shape[0]:,:]
        if (EcalBDTCut):
                testset1 = testset1[np.where((testset1[:,-5]>-2) & (testset1[:,-5] < EcalBDTCutvalue_select_proton))[0],:]
                testset2 = testset2[np.where((testset2[:,-5]>-2) & (testset2[:,-5] < EcalBDTCutvalue_select_proton))[0],:]
        if (TrdProtonHeliumCut):
                testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
                testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
        

elif Testsample is "antiproton":
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/negative/147_175GeV/pattern0/negative_147_175_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/positive/147_175GeV/pattern0/positive_147_175_pattern_0.npy',encoding="latin1")
#        testset1 = testset1[-testset2.shape[0]:,:]
        if (EcalBDTCut):
                testset1 = testset1[np.where((testset1[:,-5]>-2) & (testset1[:,-5] < EcalBDTCutvalue_select_proton))[0],:]
                testset2 = testset2[np.where((testset2[:,-5]>-2) & (testset2[:,-5] < EcalBDTCutvalue_select_proton))[0],:]
        if (TrdProtonHeliumCut):
                testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
                testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]

elif Testsample is "electron":
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1091_el.pl1.0_25200_7.6_all/147-1000GeV/transferdata/negative/147_175GeV/pattern0/negative_147_175_pattern_0.npy',encoding="latin1")
        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1091_el.pl1.0_25200_7.6_all/147-1000GeV/transferdata/positive/147_175GeV/pattern0/positive_147_175_pattern_0.npy',encoding="latin1")
#        testset1 = testset1[-testset2.shape[0]:,:]
        if (EcalBDTCut):
                testset1 = testset1[np.where((testset1[:,-5]>-2) & (testset1[:,-5] < EcalBDTCutvalue_select_electron))[0],:]
                testset2 = testset2[np.where((testset2[:,-5]>-2) & (testset2[:,-5] < EcalBDTCutvalue_select_electron))[0],:]
        if (TrdProtonHeliumCut):
                testset1 = testset1[np.where(testset1[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
                testset2 = testset2[np.where(testset2[:,TrdLogLikelihoodRatioProtonHeliumTracker] < TrdProtonHeliumCutValue)[0],:]
else:
        print("Please notify the kind of testsample correctly!")
        testset1 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/negative/147_175GeV/pattern0/negative_147_175_pattern_0.npy',encoding="latin1")

        testset2 = np.load('/hpcwork/jara0052/sichen/analysis_6.2/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/147-1000GeV/transferdata/negative/147_175GeV/pattern0/negative_147_175_pattern_0.npy',encoding="latin1")
#        testset1 = testset1[-testset2.shape[0]:,:]
##################################################################################
testset = np.r_[testset1,testset2]
testsetMvA = testset[:,MvAresult] ## get the CCProton MvA results
testset = np.expand_dims(testset, axis=2)
inputlist = np.arange(0,16)
testset = testset[:,inputlist,:]
testset[:,4] = np.fabs(testset[:,4])
testset[:,5] = np.fabs(testset[:,5])
testset[:,6] = np.fabs(testset[:,6])
features=16

##############  model prediction   #################
import CNN_models
model = CNN_models.VGG16(features)
model.load_weights('VGG16.h5')
y_pred = model.predict(testset)

######################### fraction plot ################################
protonaccumulatecount=0
fractionoprotonset=np.array([])
for protoniterate in range(binnumber_NN):
        protonaccumulatecount=protonaccumulatecount+np.histogram(y_pred[0:testset1.shape[0],1],binnumber_NN,range=(0,1))[0][binnumber_NN-protoniterate-1]
        fractionoproton=float(protonaccumulatecount)/float((testset1.shape[0]))
        fractionoprotonset=np.append(fractionoprotonset,fractionoproton)
CCprotonaccumulatecount=0
rejectionCCprotonset=np.array([])
for CCprotoniterate in range(binnumber_NN):
        CCprotonaccumulatecount=CCprotonaccumulatecount+np.histogram(y_pred[testset1.shape[0]:y_pred.shape[0],1],binnumber_NN,range=(0,1))[0][binnumber_NN-CCprotoniterate-1]
        rejectionCCproton=1.0/(float(CCprotonaccumulatecount)/float((testset2.shape[0])))
        rejectionCCprotonset=np.append(rejectionCCprotonset,rejectionCCproton)
protonaccumulatecountMvA=0
fractionoprotonsetMvA=np.array([])
for protoniterateMvA in range(binnumber_NN):
        protonaccumulatecountMvA=protonaccumulatecountMvA+np.histogram(testsetMvA[0:testset1.shape[0]],binnumber_NN,range=(0,1))[0][binnumber_NN-protoniterateMvA-1]
        fractionoprotonMvA=float(protonaccumulatecountMvA)/float((testset1.shape[0]))
        fractionoprotonsetMvA=np.append(fractionoprotonsetMvA,fractionoprotonMvA)
CCprotonaccumulatecountMvA=0
rejectionCCprotonsetMvA=np.array([])
for CCprotoniterateMvA in range(binnumber_NN):
        CCprotonaccumulatecountMvA=CCprotonaccumulatecountMvA+np.histogram(testsetMvA[testset1.shape[0]:y_pred.shape[0]],binnumber_NN,range=(0,1))[0][binnumber_NN-CCprotoniterateMvA-1]
        rejectionCCprotonMvA=1.0/(float(CCprotonaccumulatecountMvA)/float((testset2.shape[0])))
        rejectionCCprotonsetMvA=np.append(rejectionCCprotonsetMvA,rejectionCCprotonMvA)
#MLPBNN
protonaccumulatecount_MLPBNN=0
fractionoprotonset_MLPBNN=np.array([])
for protoniterate_MLPBNN in range(binnumber):
        protonaccumulatecount_MLPBNN=protonaccumulatecount_MLPBNN+np.histogram(MLPBNN_signal,binnumber,range=(0,1))[0][binnumber-protoniterate_MLPBNN-1]
        fractionoproton_MLPBNN=float(protonaccumulatecount_MLPBNN)/float((MLPBNN_signal.shape[0]))
        fractionoprotonset_MLPBNN=np.append(fractionoprotonset_MLPBNN,fractionoproton_MLPBNN)
CCprotonaccumulatecount_MLPBNN=0
rejectionCCprotonset_MLPBNN=np.array([])
for CCprotoniterate_MLPBNN in range(binnumber):
        CCprotonaccumulatecount_MLPBNN=CCprotonaccumulatecount_MLPBNN+np.histogram(MLPBNN_background,binnumber,range=(0,1))[0][binnumber-CCprotoniterate_MLPBNN-1]
        rejectionCCproton_MLPBNN=1.0/(float(CCprotonaccumulatecount_MLPBNN)/float((MLPBNN_background.shape[0])))
        rejectionCCprotonset_MLPBNN=np.append(rejectionCCprotonset_MLPBNN,rejectionCCproton_MLPBNN)
#SVM
protonaccumulatecount_SVM=0
fractionoprotonset_SVM=np.array([])
for protoniterate_SVM in range(binnumber):
        protonaccumulatecount_SVM=protonaccumulatecount_SVM+np.histogram(SVM_signal,binnumber,range=(0,1))[0][binnumber-protoniterate_SVM-1]
        fractionoproton_SVM=float(protonaccumulatecount_SVM)/float((SVM_signal.shape[0]))
        fractionoprotonset_SVM=np.append(fractionoprotonset_SVM,fractionoproton_SVM)
CCprotonaccumulatecount_SVM=0
rejectionCCprotonset_SVM=np.array([])
for CCprotoniterate_SVM in range(binnumber):
        CCprotonaccumulatecount_SVM=CCprotonaccumulatecount_SVM+np.histogram(SVM_background,binnumber,range=(0,1))[0][binnumber-CCprotoniterate_SVM-1]
        if CCprotonaccumulatecount_SVM==0:
            CCprotonaccumulatecount_SVM=1
        rejectionCCproton_SVM=1.0/(float(CCprotonaccumulatecount_SVM)/float((SVM_background.shape[0])))
        rejectionCCprotonset_SVM=np.append(rejectionCCprotonset_SVM,rejectionCCproton_SVM)
#LikelihoodKDE
protonaccumulatecount_LikelihoodKDE=0
fractionoprotonset_LikelihoodKDE=np.array([])
for protoniterate_LikelihoodKDE in range(binnumber_KDE):
        protonaccumulatecount_LikelihoodKDE=protonaccumulatecount_LikelihoodKDE+np.histogram(LikelihoodKDE_signal,binnumber_KDE,range=(0,1))[0][binnumber_KDE-protoniterate_LikelihoodKDE-1]
        fractionoproton_LikelihoodKDE=float(protonaccumulatecount_LikelihoodKDE)/float((LikelihoodKDE_signal.shape[0]))
        fractionoprotonset_LikelihoodKDE=np.append(fractionoprotonset_LikelihoodKDE,fractionoproton_LikelihoodKDE)
CCprotonaccumulatecount_LikelihoodKDE=0
rejectionCCprotonset_LikelihoodKDE=np.array([])
for CCprotoniterate_LikelihoodKDE in range(binnumber_KDE):
        CCprotonaccumulatecount_LikelihoodKDE=CCprotonaccumulatecount_LikelihoodKDE+np.histogram(LikelihoodKDE_background,binnumber_KDE,range=(0,1))[0][binnumber_KDE-CCprotoniterate_LikelihoodKDE-1]
        if CCprotonaccumulatecount_LikelihoodKDE==0:
            CCprotonaccumulatecount_LikelihoodKDE=1
        rejectionCCproton_LikelihoodKDE=1.0/(float(CCprotonaccumulatecount_LikelihoodKDE)/float((LikelihoodKDE_background.shape[0])))
        rejectionCCprotonset_LikelihoodKDE=np.append(rejectionCCprotonset_LikelihoodKDE,rejectionCCproton_LikelihoodKDE)
#
#BDT * Attention: range from -1 to 1
protonaccumulatecount_BDTG=0
fractionoprotonset_BDTG=np.array([])
for protoniterate_BDTG in range(binnumber):
        protonaccumulatecount_BDTG=protonaccumulatecount_BDTG+np.histogram(BDTG_signal,binnumber,range=(-1,1))[0][binnumber-protoniterate_BDTG-1]
        fractionoproton_BDTG=float(protonaccumulatecount_BDTG)/float((BDTG_signal.shape[0]))
        fractionoprotonset_BDTG=np.append(fractionoprotonset_BDTG,fractionoproton_BDTG)
CCprotonaccumulatecount_BDTG=0
rejectionCCprotonset_BDTG=np.array([])
for CCprotoniterate_BDTG in range(binnumber):
        CCprotonaccumulatecount_BDTG=CCprotonaccumulatecount_BDTG+np.histogram(BDTG_background,binnumber,range=(-1,1))[0][binnumber-CCprotoniterate_BDTG-1]
        rejectionCCproton_BDTG=1.0/(float(CCprotonaccumulatecount_BDTG)/float((BDTG_background.shape[0])))
        rejectionCCprotonset_BDTG=np.append(rejectionCCprotonset_BDTG,rejectionCCproton_BDTG)
#
#Fisher * Attention: range from -1 to 1
protonaccumulatecount_Fisher=0
fractionoprotonset_Fisher=np.array([])
for protoniterate_Fisher in range(binnumber):
        protonaccumulatecount_Fisher=protonaccumulatecount_Fisher+np.histogram(Fisher_signal,binnumber,range=(-1,1))[0][binnumber-protoniterate_Fisher-1]
        fractionoproton_Fisher=float(protonaccumulatecount_Fisher)/float((Fisher_signal.shape[0]))
        fractionoprotonset_Fisher=np.append(fractionoprotonset_Fisher,fractionoproton_Fisher)
CCprotonaccumulatecount_Fisher=0
rejectionCCprotonset_Fisher=np.array([])
for CCprotoniterate_Fisher in range(binnumber):
        CCprotonaccumulatecount_Fisher=CCprotonaccumulatecount_Fisher+np.histogram(Fisher_background,binnumber,range=(-1,1))[0][binnumber-CCprotoniterate_Fisher-1]
        rejectionCCproton_Fisher=1.0/(float(CCprotonaccumulatecount_Fisher)/float((Fisher_background.shape[0])))
        rejectionCCprotonset_Fisher=np.append(rejectionCCprotonset_Fisher,rejectionCCproton_Fisher)
#
#Likelihood * Attention: range from -1 to 1
protonaccumulatecount_Likelihood=0
fractionoprotonset_Likelihood=np.array([])
for protoniterate_Likelihood in range(binnumber):
        protonaccumulatecount_Likelihood=protonaccumulatecount_Likelihood+np.histogram(Likelihood_signal,binnumber,range=(-1,1))[0][binnumber-protoniterate_Likelihood-1]
        fractionoproton_Likelihood=float(protonaccumulatecount_Likelihood)/float((Likelihood_signal.shape[0]))
        fractionoprotonset_Likelihood=np.append(fractionoprotonset_Likelihood,fractionoproton_Likelihood)
CCprotonaccumulatecount_Likelihood=0
rejectionCCprotonset_Likelihood=np.array([])
for CCprotoniterate_Likelihood in range(binnumber):
        CCprotonaccumulatecount_Likelihood=CCprotonaccumulatecount_Likelihood+np.histogram(Likelihood_background,binnumber,range=(-1,1))[0][binnumber-CCprotoniterate_Likelihood-1]
        rejectionCCproton_Likelihood=1.0/(float(CCprotonaccumulatecount_Likelihood)/float((Likelihood_background.shape[0])))
        rejectionCCprotonset_Likelihood=np.append(rejectionCCprotonset_Likelihood,rejectionCCproton_Likelihood)
#
plt.figure(figsize=(18,9))
#plt.ylim(1, rejectionCCprotonsetMvA[0]+2)
plt.plot(fractionoprotonset, rejectionCCprotonset, 'g-',lw=3, label='Neural Network')
#plt.plot(fractionoprotonsetMvA, rejectionCCprotonsetMvA, 'b-', lw=3,label='MVA')
plt.plot(fractionoprotonset_BDTG, rejectionCCprotonset_BDTG, 'm-', lw=3,label='BDTG')
plt.plot(fractionoprotonset_MLPBNN, rejectionCCprotonset_MLPBNN, 'r-', lw=3,label='MLPBNN')
plt.plot(fractionoprotonset_SVM[-135:], rejectionCCprotonset_SVM[-135:], 'y-', lw=3,label='SVM')
plt.plot(fractionoprotonset_LikelihoodKDE, rejectionCCprotonset_LikelihoodKDE, 'c-', lw=3,label='LikelihoodKDE')
plt.plot(fractionoprotonset_Fisher, rejectionCCprotonset_Fisher, 'k-', lw=3,label='Fisher')
plt.plot(fractionoprotonset_Likelihood, rejectionCCprotonset_Likelihood, 'b-', lw=3,label='Likelihood')
#plt.yscale('log')
plt.xlim(left=0.4)
#plt.ylim(0.0,20000.0)
plt.xlabel('Efficiency of Proton',fontsize=22)
plt.ylabel('Rejection for CCproton',fontsize=22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.legend( loc='best',fontsize=20)
plt.grid(True)
plt.savefig('fractionfigure.png')

########################## NN Prediction ##################################
plt.figure(figsize=(18,18))
plt.hist(y_pred[0:testset1.shape[0],1],bins=binnumber,range=(0,1),alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black' )
plt.hist(y_pred[testset1.shape[0]:y_pred.shape[0],1],bins=binnumber,range=(0,1),alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black' )
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Estimator $_{CC}$',fontsize=22)
plt.ylabel('Count',fontsize=22)
plt.legend(loc='upper center',fontsize=30)
plt.savefig('ML_test.png')

########################  MvA Prediction ##############################
plt.figure(figsize=(18,18))
plt.hist(testsetMvA[0:testset1.shape[0]],binnumber,range=(0,1),normed=True,log=True, alpha=0.5,label='ChargeCorrect',facecolor='blue',edgecolor='black'  )
plt.hist(testsetMvA[testset1.shape[0]:y_pred.shape[0]],binnumber,range=(0,1),normed=True,log=True, alpha=0.5,label='ChargeConfused',facecolor='green',edgecolor='black'  )
plt.xlabel('Estimator $_{CC}$',fontsize=30)
plt.ylabel('Count',fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.legend(loc='upper center',fontsize=30)
plt.savefig('CCMVA_test.png')
'''
###################### ROOT PLOT #################################################################################################################
ML_p = TH1D("ML_p","", 100,0,1)
ML_n = TH1D("ML_n","", 100,0,1)
CCMVA_p = TH1D("CCMVA_p","", 100,0,1)
CCMVA_n = TH1D("CCMVA_n","", 100,0,1)

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
ML_p.Draw()
ML_n.Draw('same')
leg =TLegend(.4,.7,.6,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
ML_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c1.Update()
c1.SaveAs("ML_seperation.pdf")

c2 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
CCMVA_p.Draw()
CCMVA_n.Draw('same')
leg =TLegend(.4,.7,.6,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
CCMVA_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c2.Update()
c2.SaveAs("CCMVA_seperation.pdf")

c3 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
ML_p.Draw()
ML_n.Draw('same')
leg =TLegend(.4,.7,.6,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
ML_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c3.SetLogy()
c3.Update()
c3.SaveAs("ML_seperation_log.pdf")

c4 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
CCMVA_p.Draw()
CCMVA_n.Draw('same')
leg =TLegend(.4,.7,.6,.9,)
leg.SetFillColor(0)
leg.AddEntry(ML_p,"Charge Correct")
leg.AddEntry(ML_n,"Charge Confused")
leg.Draw()
CCMVA_p.GetXaxis().SetTitle("Charge Confusion Estimator");
c4.SetLogy()
c4.Update()
c4.SaveAs("CCMVA_seperation_log.pdf")
'''




