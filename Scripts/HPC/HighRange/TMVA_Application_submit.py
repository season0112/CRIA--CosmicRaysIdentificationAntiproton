#!/usr/bin/env python
from __future__ import division
import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from mpl_toolkits.mplot3d import Axes3D
from root_numpy import root2array, tree2array, fill_hist, hist2array
from ROOT import TFile, TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
import binning
plt.switch_backend('agg')
import CC_Application_tool 


def FurtherCut(CutArray, ArraySpice, EcalBDTCut, TrdProtonHeliumCut, PhysicstriggerCut):
    # ArraySpice: ProtonLike(Correct and Confused), ElectronLike
    if (EcalBDTCut):
        if ArraySpice == 'ProtonLike':
            CutArray = CutArray[np.where((CutArray['EcalBDT_EnergyD']>-2) & (CutArray['EcalBDT_EnergyD']<EcalBDTCutvalue_select_proton))[0]]
        elif ArraySpice == 'ElectronLike':
            CutArray = CutArray[np.where((CutArray['EcalBDT_EnergyD']>-2) & (CutArray['EcalBDT_EnergyD']>EcalBDTCutvalue_select_electron))[0]] # select electron from EcalBDT
    if (TrdProtonHeliumCut):
        CutArray = CutArray[np.where((CutArray['TrdLogLikelihoodRatioProtonHeliumTracker'] < TrdProtonHeliumCutValue))[0]] # to be optimized.
    if (PhysicstriggerCut):
        CutArray = CutArray[np.where((CutArray["TriggerFlags"]&0x3e)>0)[0]]
    return CutArray 


def RawArrray_To_ReducedArray(RawArray, RawArraySign):
    reduced_dtype = np.dtype( [ ('TrdLogLikelihoodRatioElectronProtonTracker','<f4'), ('RigidityAsymmetry','<f4'), ('RigidityAsymmetryL9','<f4'), ('Chi2TrackerYAsymmetry','<f4'), ('InnerMaxSpanRigidityMatching','<f4'), ('L1L9RigidityMatching','<f4'), ('L24L58RigidityMatching','<f4'), ('Log10Chi2TrackerXInner','<f4'), ('Log10Chi2TrackerYInner','<f4'), ('Log10Chi2TrackerX','<f4'), ('Log10Chi2TrackerY','<f4'), ('TrackerL58L24ChargeAsymmetry','<f4'), ('TrackerL9Charge','<f4'), ('TrackerL78Charge','<f4'), ('UpperTofCharge','<f4'), ('LowerTofCharge','<f4') ] )
    ReducedArray = np.transpose( [RawArray['TrdLogLikelihoodRatioElectronProtonTracker'], RawArray['RigidityAsymmetry'], RawArray['RigidityAsymmetryL9'], RawArray['Chi2TrackerYAsymmetry'], RawArray['InnerMaxSpanRigidityMatching']*RawArraySign, RawArray['L1L9RigidityMatching']*RawArraySign, RawArray['L24L58RigidityMatching']*RawArraySign, RawArray['Log10Chi2TrackerXInner'], RawArray['Log10Chi2TrackerYInner'], RawArray['Log10Chi2TrackerX'], RawArray['Log10Chi2TrackerY'], RawArray['TrackerL58L24ChargeAsymmetry'], RawArray['TrackerL9Charge'], RawArray['TrackerL78Charge'], RawArray['UpperTofCharge'], RawArray['LowerTofCharge']] )
    ReducedArray = np.asarray(ReducedArray, order='C')
    ReducedArray.dtype = reduced_dtype
    return ReducedArray, RawArray['Weight'], RawArray['MCPrimaryMomentum']


def CutOnPatterns(PorcessedArray, Pattern):
    PorcessedArray = PorcessedArray[np.where(PorcessedArray['Pattern'] == int(Pattern))[0]]
    return PorcessedArray 


def Plot(InputArray, RigidityBin, pdfName, Pattern):
    plt.figure()
    plt.hist(InputArray, 50, facecolor='g', alpha=0.5)
    plt.yscale('log')
    plt.xlabel('CC estimator')
    plt.ylabel('Counts')
    plt.title('')
    plt.savefig( "/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/CC_Estimator_Check/" + pdfName + "_" + str(RigidityBin) + "_Patttern_" + str(Pattern) + ".pdf" )
    plt.close()


def main():

    if RigidityBinSection == "All":
        RigidityBinList = ["14.1_15.3", "15.3_16.6", "16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175", "175_211", "211_259", "211_250", "250_330", "259_450", "330_525"] 

    elif RigidityBinSection == "14.1_15.3":
        RigidityBinList = ["14.1_15.3"]
    elif RigidityBinSection == "15.3_16.6":
        RigidityBinList = ["15.3_16.6"]
    elif RigidityBinSection == "16.6_18":
        RigidityBinList = ["16.6_18"]
    elif RigidityBinSection == "18_19.5":
        RigidityBinList = ["18_19.5"]
    elif RigidityBinSection == "19.5_21.1":
        RigidityBinList = ["19.5_21.1"]
    elif RigidityBinSection == "21.1_22.8":
        RigidityBinList = ["21.1_22.8"]
    elif RigidityBinSection == "22.8_24.7":
        RigidityBinList = ["22.8_24.7"]
    elif RigidityBinSection == "24.7_26.7":
        RigidityBinList = ["24.7_26.7"]
    elif RigidityBinSection == "26.7_28.8":
        RigidityBinList = ["26.7_28.8"]
    elif RigidityBinSection == "28.8_31.1":
        RigidityBinList = ["28.8_31.1"]
    elif RigidityBinSection == "31.1_33.5":
        RigidityBinList = ["31.1_33.5"]
    elif RigidityBinSection == "33.5_36.1":
        RigidityBinList = ["33.5_36.1"]
    elif RigidityBinSection == "36.1_38.9":
        RigidityBinList = ["36.1_38.9"]
    elif RigidityBinSection == "38.9_41.9":
        RigidityBinList = ["38.9_41.9"]
    elif RigidityBinSection == "41.9_45.1":
        RigidityBinList = ["41.9_45.1"]
    elif RigidityBinSection == "45.1_48.5":
        RigidityBinList = ["45.1_48.5"]
    elif RigidityBinSection == "48.5_52.2":
        RigidityBinList = ["48.5_52.2"]
    elif RigidityBinSection == "52.2_56.1":
        RigidityBinList = ["52.2_56.1"]
    elif RigidityBinSection == "56.1_60.3":
        RigidityBinList = ["56.1_60.3"]
    elif RigidityBinSection == "60.3_64.8":
        RigidityBinList = ["60.3_64.8"]
    elif RigidityBinSection == "64.8_69.7":
        RigidityBinList = ["64.8_69.7"]
    elif RigidityBinSection == "69.7_74.9":
        RigidityBinList = ["69.7_74.9"]
    elif RigidityBinSection == "74.9_80.5":
        RigidityBinList = ["74.9_80.5"]
    elif RigidityBinSection == "80.5_93":
        RigidityBinList = ["80.5_93"]
    elif RigidityBinSection == "93_108":
        RigidityBinList = ["93_108"]
    elif RigidityBinSection == "108_125":
        RigidityBinList = ["108_125"]
    elif RigidityBinSection == "125_147":
        RigidityBinList = ["125_147"]
    elif RigidityBinSection == "147_175":
        RigidityBinList = ["147_175"]
    elif RigidityBinSection == "175_211":
        RigidityBinList = ["175_211"]
    elif RigidityBinSection == "211_259":
        RigidityBinList = ["211_259"]
    elif RigidityBinSection == "211_250":
        RigidityBinList = ["211_250"]
    elif RigidityBinSection == "250_330":
        RigidityBinList = ["250_330"]
    elif RigidityBinSection == "259_450":
        RigidityBinList = ["259_450"]
    elif RigidityBinSection == "330_525":
        RigidityBinList = ["330_525"]
    else:
        print("Wrong Choice of RigidityBinSection ! ")


    for RigidityBin in RigidityBinList:
        print("Now is: " + str(RigidityBin))

        #### Load data
        ## 1. MC data as choices of templates
        #Correct_MC_name  = 'B1042_antipr.pl1.1800_7.6_all_Tree'                                        # Default
        #Correct_MC_name  = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive'                      # Large (For CCLevel)
        #Confused_MC_name = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative'                     # Default
        #Electron_MC_name = 'B1091_el.pl1.0_25_2000_7.6_all_Tree'                                      # Default

        ## For CCLevel
        #Correct_MC_name  = 'B1042_pr.pl1.1800_7.6_all_Tree_positive'                                  # Test for CCLevel
        #Confused_MC_name = 'B1042_pr.pl1.1800_7.6_all_Tree_negative'                                  # Test for CCLevel
        #Correct_MC_name  = 'B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_positive'           # Test for CCLevel
        #Confused_MC_name = 'B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_negative'           # Test for CCLevel
        #Correct_MC_name  = 'B1220_pr.pl1ph.021000_7.8_all_Tree_positive'                               # Test for CCLevel
        #Confused_MC_name = 'B1220_pr.pl1ph.021000_7.8_all_Tree_negative'                               # Test for CCLevel
        #Correct_MC_name  = 'B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_Tree_positive'                  # Test for CCLevel
        #Confused_MC_name = 'B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_Tree_negative'                  # Test for CCLevel

        #Correct_MC = root2array(datapath + '/' + Correct_MC_name + '_' + RigidityBin + '.root', 'ExampleAnalysisTree')
        #Confused_MC = root2array(datapath + '/' + Confused_MC_name + '_' + RigidityBin + '.root', 'ExampleAnalysisTree')
        #Electron_MC = root2array(datapath + '/' + Electron_MC_name + '_' + RigidityBin + '.root', 'ExampleAnalysisTree')

        ## 2. ISS data 
        ISSpositive = root2array(datapath + '/B1130_pass7_7.8_all_Tree_positive' + suffix + RigidityBin + '.root', 'ExampleAnalysisTree')
        ISSnegative = root2array(datapath + '/B1130_pass7_7.8_all_Tree_negative' + suffix + RigidityBin + '.root', 'ExampleAnalysisTree')
        Electron_Template_Data = root2array(datapath + '/B1130_pass7_7.8_all_Tree_negative' + suffix + RigidityBin + '.root', 'ExampleAnalysisTree', selection="EcalBDT_EnergyD>-999 && EcalBDT_EnergyD>0" ) 
        
        #### cut on EcalBDT, TrdProtonHeliumCut, and Physicstrigger
        if Pattern == 0:
            EcalBDTCut = True
        else:
            EcalBDTCut = False

        ISSpositive            = FurtherCut(ISSpositive           , "ProtonLike"  , EcalBDTCut, TrdProtonHeliumCut, PhysicstriggerCut)
        ISSnegative            = FurtherCut(ISSnegative           , "ProtonLike"  , EcalBDTCut, TrdProtonHeliumCut, PhysicstriggerCut)
        #Correct_MC             = FurtherCut(Correct_MC            , "ProtonLike"  , EcalBDTCut, TrdProtonHeliumCut, PhysicstriggerCut)
        #Confused_MC            = FurtherCut(Confused_MC           , "ProtonLike"  , EcalBDTCut, TrdProtonHeliumCut, PhysicstriggerCut)
        #Electron_MC            = FurtherCut(Electron_MC           , "ElectronLike", EcalBDTCut, TrdProtonHeliumCut, PhysicstriggerCut)
        Electron_Template_Data = FurtherCut(Electron_Template_Data, "ElectronLike", EcalBDTCut, TrdProtonHeliumCut, PhysicstriggerCut)

        #### Cut on Tracker Patterns
        ISSpositive            = CutOnPatterns(ISSpositive           , Pattern)
        ISSnegative            = CutOnPatterns(ISSnegative           , Pattern)
        #Correct_MC             = CutOnPatterns(Correct_MC            , Pattern)
        #Confused_MC            = CutOnPatterns(Confused_MC           , Pattern)
        #Electron_MC            = CutOnPatterns(Electron_MC           , Pattern)
        Electron_Template_Data = CutOnPatterns(Electron_Template_Data, Pattern)


        #### Get Rigidity Sign to correct some MVA varaibles.
        sign_ISSPositive               = ISSpositive['Rigidity']/np.abs(ISSpositive['Rigidity'])
        sign_ISSNegative               = ISSnegative['Rigidity']/np.abs(ISSnegative['Rigidity'])
        #sign_Correct_MC                = Correct_MC['Rigidity']/np.abs(Correct_MC['Rigidity'])
        #sign_Confused_MC               = Confused_MC['Rigidity']/np.abs(Confused_MC['Rigidity'])
        #sign_Electron_MC               = Electron_MC['Rigidity']/np.abs(Electron_MC['Rigidity'])
        sign_Electron_Template_Data    = Electron_Template_Data['Rigidity']/np.abs(Electron_Template_Data['Rigidity'])


        #### From all variables to only MVA variables
        ISSpositive_simplified           , ISSpositive_Weight           , ISSpositive_MCPrimaryMomentum            = RawArrray_To_ReducedArray(ISSpositive, sign_ISSPositive)
        ISSnegative_simplified           , ISSnegative_Weight           , ISSnegative_MCPrimaryMomentum            = RawArrray_To_ReducedArray(ISSnegative, sign_ISSNegative)
        #Correct_MC_simplified            , Correct_MC_Weight            , Correct_MC_MCPrimaryMomentum             = RawArrray_To_ReducedArray(Correct_MC , sign_Correct_MC)
        #Confused_MC_simplified           , Confused_MC_Weight           , Confused_MC_MCPrimaryMomentum            = RawArrray_To_ReducedArray(Confused_MC, sign_Confused_MC)
        #Electron_MC_simplified           , Electron_MC_Weight           , Electron_MC_MCPrimaryMomentum            = RawArrray_To_ReducedArray(Electron_MC, sign_Electron_MC)
        Electron_Template_Data_simplified, Electron_Template_Data_Weight, Electron_Template_Data_MCPrimaryMomentum = RawArrray_To_ReducedArray(Electron_Template_Data, sign_Electron_Template_Data)

        #### Release some memories in case of out-of-memory error.
        del sign_ISSPositive
        del sign_ISSNegative
        #del sign_Correct_MC
        #del sign_Confused_MC
        #del sign_Electron_MC
        del sign_Electron_Template_Data

        del ISSpositive
        del ISSnegative
        #del Correct_MC
        #del Confused_MC
        #del Electron_MC
        del Electron_Template_Data

        CC_Application_tool.TMVAEvaluate(ISSpositive_simplified           , ISSpositive_Weight           , ISSpositive_MCPrimaryMomentum           , Pattern, RigidityBin, "ISS_positive" + suffix                                         , datapath)
        CC_Application_tool.TMVAEvaluate(ISSnegative_simplified           , ISSnegative_Weight           , ISSnegative_MCPrimaryMomentum           , Pattern, RigidityBin, "ISS_negative" + suffix                                         , datapath)
        #CC_Application_tool.TMVAEvaluate(Correct_MC_simplified            , Correct_MC_Weight            , Correct_MC_MCPrimaryMomentum            , Pattern, RigidityBin, "ChargeCorrectProtonTemplate_MC_" + Correct_MC_name + str("_")  , datapath)
        #CC_Application_tool.TMVAEvaluate(Confused_MC_simplified           , Confused_MC_Weight           , Confused_MC_MCPrimaryMomentum           , Pattern, RigidityBin, "ChargeConfusedProtomTemplate_MC_" + Confused_MC_name + str("_"), datapath)
        #CC_Application_tool.TMVAEvaluate(Electron_MC_simplified           , Electron_MC_Weight           , Electron_MC_MCPrimaryMomentum           , Pattern, RigidityBin, "ElectronTemplate_MC_" + Electron_MC_name + str("_")            , datapath)
        CC_Application_tool.TMVAEvaluate(Electron_Template_Data_simplified, Electron_Template_Data_Weight, Electron_Template_Data_MCPrimaryMomentum, Pattern, RigidityBin, "ElectronTemplate_Data_" + suffix                               , datapath)

        '''
        #### Release some memories in case of out-of-memory error again.
        #del ISSpositive_simplified 
        #del ISSnegative_simplified
        del Correct_MC_simplified
        del Confused_MC_simplified
        #del Electron_MC_simplified
        #del Electron_Template_Data_simplified

        #del ISSpositive_Weight
        #del ISSnegative_Weight           
        del Correct_MC_Weight
        del Confused_MC_Weight
        #del Electron_MC_Weight
        #del Electron_Template_Data_Weight
        '''

        '''
        #### Plot response
        Proton   = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/ChargeCorrectProtonTemplate_MC_B1042_antipr.pl1.1800_7.6_all_Tree_" + str(RigidityBin) + "_Pattern_" + str(Pattern) + ".npy")
        CCProton = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_" + str(RigidityBin) + "_Pattern_" + str(Pattern) + ".npy")
        Electron = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/ElectronTemplate_MC_B1091_el.pl1.0_25_2000_7.6_all_Tree_" + str(RigidityBin) + "_Pattern_" + str(Pattern) + ".npy")

        Plot(Proton[:,0]  , RigidityBin, "Proton"  , Pattern)
        Plot(CCProton[:,0], RigidityBin, "CCProton", Pattern)
        Plot(Electron[:,0], RigidityBin, "Electron", Pattern)
        '''


if __name__ == '__main__':
 
    #### Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--issversion', help='ISS data version, pass7.8 or published2016')
    parser.add_argument('--pattern', help='Which tracker patterns you choose')

    parser.add_argument('--rigidityBinSection', default="All", help='Since some rigidity bins are large, so this help to train seperately.')

    parser.add_argument('--trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_true',help='trdProtonHeliumCut. Default:False')
    parser.add_argument('--no-trdProtonHeliumCut', dest='trdProtonHeliumCut', action='store_false',help='trdProtonHeliumCut. Default:False')
    parser.set_defaults(trdProtonHeliumCut=False)
    parser.add_argument('--physicstriggerCut', dest='physicstriggerCut', action='store_true',help='physicstriggerCut. Default:False')
    parser.add_argument('--no-physicstriggerCut', dest='physicstriggerCut', action='store_false',help='physicstriggerCut. Default:False')
    parser.set_defaults(physicstriggerCut=False)

    arguments = parser.parse_args()

    RigidityBinSection = arguments.rigidityBinSection

    if (arguments.issversion):
        ISSversion = arguments.issversion
    else:
        print("You need to choose a ISS data cersion !")
        os._exit(0)

    if (arguments.pattern):
        Pattern = arguments.pattern
    else:
        print("You need to choose a tracker pattern !")
        os._exit(0)

    if (arguments.trdProtonHeliumCut):
        TrdProtonHeliumCut = arguments.trdProtonHeliumCut
    else:
        TrdProtonHeliumCut = False

    if (arguments.physicstriggerCut):
        PhysicstriggerCut = arguments.physicstriggerCut
    else:
        PhysicstriggerCut = False

    #### Free Parameters
    EcalBDTCutvalue_select_proton = 0.0  ## -1 denote proton, 1 denote electron
    EcalBDTCutvalue_select_electron = 0.0  ## -1 denote proton, 1 denote electron
    TrdProtonHeliumCutValue = 0.3
    datapath = os.getenv('HPCHIGHENERGYDATADIR')


    if ISSversion == "pass7.8":
        suffix = "_"
    elif ISSversion == "published2016":
        suffix = "_May2015_"
    elif ISSversion == "PhyRep2021":
        suffix = "_Nov2017_"

    main()


