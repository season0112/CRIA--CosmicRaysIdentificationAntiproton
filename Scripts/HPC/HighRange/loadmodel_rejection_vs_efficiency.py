#!/usr/bin/env python
########################################################################
#                                                                      #
#                   Pronton Charge Confusion with VGG16                #
#                                  03.05.2019                          #
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
from root_numpy import fill_hist, root2array, tree2array

parser = argparse.ArgumentParser()
parser.add_argument('--testsample', choices=["proton", "antiproton", "electron", "negative_events","validation_sample","B1119_ele"], help='validation_data_choice')
parser.add_argument('--binnumber',  type=int, help='binnumber for training. Default:80')
parser.add_argument('--ecalBDTCut', dest='ecalBDTCut', action='store_true',help='ecalBDTCut. Default:False')
parser.add_argument('--no-ecalBDTCut', dest='ecalBDTCut', action='store_false',help='ecalBDTCut. Default:False')
parser.set_defaults(ecalBDTCut=False)
parser.add_argument('--trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_true',help='trdProtonHeliumCut. Default:False')
parser.add_argument('--no-trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_false',help='trdProtonHeliumCut. Default:False')
parser.set_defaults(trdProtonHeliumCut=False)
arguments = parser.parse_args()

NN_rejection_all_eff = np.array([])
MVA_rejection_all_eff = np.array([])

########################## Free Parameters #############################
if (arguments.testsample):
    Testsample = arguments.testsample
else:
    print("You need to choose a testdataset! Potential choises are: proton, antiproton, electron")
    os._exit(0)

if (arguments.binnumber):
    binnumber = arguments.binnumber
else:
    binnumber = 80

EcalBDTCut = arguments.ecalBDTCut
TrdProtonHeliumCut = arguments.trdProtonHeliumCut
highpath = os.getenv('HPCHIGHENERGYDATADIR')

########################################################################
EcalBDTCutvalue_select_proton = 0.0  ## -1 denote proton, 1 denote electron
EcalBDTCutvalue_select_electron = 0.0  ## -1 denote proton, 1 denote electron
TrdProtonHeliumCutValue = 0.3

########################## Properties Label ############################
ElectronCCMVABDT = 36
MvAresult = 37
TrdLogLikelihoodRatioProtonHeliumTracker = 38
EcalBDT_EnergyD = 39

Efficiencyall = [0.8,0.9]
for Efficiency in Efficiencyall:
    NN_rejection = np.array([])
    MVA_rejection = np.array([])
    for Energybin in ["38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175", "175_211", "211_259", "259_330", "330_525",]:
            ######################### Test Sample #################################
            if Testsample == "electron":
                testset1 = root2array(highpath + '/B1091_el.pl1.0_25_2000_7.6_all_Tree_' + Energybin + '.root')
                testset2 = XXX
                testset1 = testset1[-testset2.shape[0]:]
                if (EcalBDTCut):
                    testset1 = testset1[np.where((testset1['EcalBDT_EnergyD']>-2) & (testset1['EcalBDT_EnergyD'] > EcalBDTCutvalue_select_electron))[0]]
                    testset2 = testset2[np.where((testset2['EcalBDT_EnergyD']>-2) & (testset2['EcalBDT_EnergyD'] > EcalBDTCutvalue_select_electron))[0]]
                if (TrdProtonHeliumCut):
                    testset1 = testset1[np.where(testset1['TrdLogLikelihoodRatioProtonHeliumTracker'] < TrdProtonHeliumCutValue)[0]]
                    testset2 = testset2[np.where(testset2['TrdLogLikelihoodRatioProtonHeliumTracker'] < TrdProtonHeliumCutValue)[0]]
                testset1 = testset1[np.where(testset1['ElectronCCMVABDT'] > -2.0)[0]]
                testset2 = testset2[np.where(testset2['ElectronCCMVABDT'] > -2.0)[0]]

            elif Testsample == "proton":
                testset1 = root2array(highpath + '/B1042_antipr.pl1.1800_7.6_all_Tree_' + Energybin + '.root')
                #testset1 = root2array(highpath + '/B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_' + Energybin + '.root')
                testset2 = root2array(highpath + '/B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_' + Energybin + '.root')
                testset1 = testset1[-testset2.shape[0]:]
                if (EcalBDTCut):
                    testset1 = testset1[np.where((testset1['EcalBDT_EnergyD']>-2) & (testset1['EcalBDT_EnergyD'] < EcalBDTCutvalue_select_proton))[0]]
                    testset2 = testset2[np.where((testset2['EcalBDT_EnergyD']>-2) & (testset2['EcalBDT_EnergyD'] < EcalBDTCutvalue_select_proton))[0]]
                if (TrdProtonHeliumCut):
                    testset1 = testset1[np.where(testset1['TrdLogLikelihoodRatioProtonHeliumTracker'] < TrdProtonHeliumCutValue)[0] ]
                    testset2 = testset2[np.where(testset2['TrdLogLikelihoodRatioProtonHeliumTracker'] < TrdProtonHeliumCutValue)[0] ]
                testset1 = testset1[np.where(testset1['ProtonCCMVABDT'] > -2.0)[0]]
                testset2 = testset2[np.where(testset2['ProtonCCMVABDT'] > -2.0)[0]]

            else:
                print("Please choose a kind of validation sample correctly!")
                os._exit(0)


            ######################### Test Target Sample  ##########################
            testset = np.append(testset1,testset2)
            if Testsample == "electron" or Testsample == "B1119_ele" : 
                testsetMvA = testset['ElectronCCMVABDT'] ## get the Electron Charge Confusion BDT results
            else:
                testsetMvA = testset['ProtonCCMVABDT'] ## get the Proton Charge Confusion MvA results

            sign_test = testset['Rigidity'] / np.abs(testset['Rigidity'])
            testset = np.array([ testset['TrdLogLikelihoodRatioElectronProtonTracker'], testset['RigidityAsymmetry'], testset['RigidityAsymmetryL9'], testset['Chi2TrackerYAsymmetry'], testset['InnerMaxSpanRigidityMatching']*sign_test, testset['L1L9RigidityMatching']*sign_test, testset['L24L58RigidityMatching']*sign_test, testset['Log10Chi2TrackerXInner'], testset['Log10Chi2TrackerYInner'], testset['Log10Chi2TrackerX'], testset['Log10Chi2TrackerY'], testset['TrackerL58L24ChargeAsymmetry'], testset['TrackerL9Charge'], testset['TrackerL78Charge'], testset['UpperTofCharge'], testset['LowerTofCharge']] )
            testset = np.transpose(testset)
            testset = np.expand_dims(testset, axis=2)


            ##############  model prediction   #################
            import CNN_models
            model = CNN_models.VGG16(16)
            if Energybin == "80.5_93.0" or Energybin == "93.0_108.0" or Energybin == "108.0_125.0" or Energybin == "125.0_147.0" or Energybin == "147_175" or Energybin == "175_211" or Energybin == "211_259" or Energybin == "259_330" or Energybin == "330_525":
                model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/mva_variables/VGG16_'+Energybin+'_mva_variables'+'.h5')
            else:
                model.load_weights('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/weight/mva_variables/VGG16_'+'38.9_80.5'+'_mva_variables'+'.h5')

            y_pred = model.predict(testset)


            ######################### fraction plot ################################
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
                if CCprotonaccumulatecount == 0 :
                    CCprotonaccumulatecount = 0.0000000000000001
                rejectionCCproton = 1.0/(CCprotonaccumulatecount/(testset2.shape[0]))
                rejectionCCprotonset = np.append(rejectionCCprotonset,rejectionCCproton)
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
                    if CCelectronaccumulatecountMvA == 0 :
                        CCelectronaccumulatecountMvA = 0.0000000000000001
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

            index_NN = np.argwhere( fractionoprotonset == (min(fractionoprotonset-Efficiency, key=lambda x: abs(x)) + Efficiency) )
            if Testsample == "electron" or Testsample == "B1119_ele":
                index_MVA = np.argwhere( fractionoelectronsetMvA == (min(fractionoelectronsetMvA-Efficiency, key=lambda x: abs(x)) + Efficiency) )
            elif Testsample == "proton":
                index_MVA = np.argwhere( fractionoprotonsetMvA == (min(fractionoprotonsetMvA-Efficiency, key=lambda x: abs(x)) + Efficiency) )

            NN_rejection = np.append(NN_rejection,rejectionCCprotonset[index_NN[0,0]])
            if Testsample == "electron" or Testsample == "B1119_ele":
                MVA_rejection = np.append(MVA_rejection,rejectionCCelectronsetMvA[index_MVA[0,0]])
            elif Testsample == "proton":
                MVA_rejection = np.append(MVA_rejection,rejectionCCprotonsetMvA[index_MVA[0,0]])  


    if NN_rejection_all_eff.shape[0] == 0 :
        NN_rejection_all_eff = NN_rejection
        MVA_rejection_all_eff = MVA_rejection
    else:
        NN_rejection_all_eff = np.row_stack((NN_rejection_all_eff, NN_rejection))
        MVA_rejection_all_eff = np.row_stack((MVA_rejection_all_eff, MVA_rejection))

NNbinnings = np.array([38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93, 108, 125, 147, 175, 211, 250, 330, 525])
NNpoint=np.array( np.zeros(NNbinnings.shape[0]-1) )
for i in range(NNpoint.shape[0]):
    NNpoint[i]=(NNbinnings[i]+NNbinnings[i+1])/2



##################################################################################
###################################### Plot ######################################
##################################################################################
plt.figure(figsize=(18,9))
for i in range(NN_rejection_all_eff.shape[0]):
    plt.plot(NNpoint[-15:], NN_rejection_all_eff[i,-15:], 'o-',lw=3, label='Neural Network (eff=' + str(Efficiencyall[i])+')')
    if Testsample == "electron" or Testsample == "B1119_ele":
        plt.plot(NNpoint[-15:], MVA_rejection_all_eff[i,-15:], 'o-',lw=3, label='Electron CC BDT (eff='+ str(Efficiencyall[i])+')')
    elif Testsample == "proton":
        plt.plot(NNpoint[-15:], MVA_rejection_all_eff[i,-15:], 'o-',lw=3, label='Proton CC MLPBNN (eff='+ str(Efficiencyall[i])+')')

plt.legend(loc='best',fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.xlabel('Rigidity (GV)',fontsize=35)
plt.ylabel('Rejection',fontsize=35)
if Testsample == "electron" or Testsample == "B1119_ele":
    plt.savefig(highpath + '/RejectionPlots/Electron_Compare.png')
elif Testsample == "proton":
    plt.savefig(highpath + '/RejectionPlots/Proton_Compare.png')

'''
plt.figure(figsize=(18,9))
#plt.ylim(1, rejectionCCprotonsetMvA[0]+2)
plt.plot(fractionoprotonset, rejectionCCprotonset, 'g-',lw=3, label='Neurual Network on validation sample')
if Testsample == "electron" or "B1119_ele":
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
plt.savefig('plot/rejection_' + Energybin + '_mva_variables' + '.png')
'''

