#!/usr/bin/env python
#To do: GeneralBinNumber, Last2BinNumber

#GeneralBinNumber = "CCN20TRDN12"
#Last2BinNumber = "CCN9TRDN11"
#Last2BinNumber = "CCN20TRDN12"

from __future__ import division
import numpy as np
import math
import json
import collections
import matplotlib.pyplot as plt
import argparse
import os
import binning
from root_numpy import root2array, tree2array
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import uncertainties
from uncertainties import unumpy
from matplotlib.ticker import ScalarFormatter
import CCLevel_function
import PythonPlotDefaultParameters

def main():
    ## Here no CCcut applied???
    #cclevelMCLow , cclevelErrorMCLow, cclevelMCLow_allPattern, cclevelErrorMCLow_allPattern = CCLevel_function.LoadCCLevelMCLow(lowpath  , binnameLow)    // Allpatterrns and single pattern
    #cclevelMCHigh, cclevelErrorMCHigh                                                       = CCLevel_function.LoadCCLevelMCHigh(highpath, binname)       // Allpatterrns

    f_CCLevel_ISSOverMC = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + "newall.root", "RECREATE")

    CCNumber, TRDNumber, CCcut_TF_all = CCLevel_function.FitSettingUsage(pattern)

    #CCcut_TF_all = ["0.20", "0.65", "0.70", "0.80"]
    #CCcut_TF_all = ["0.00", "0.20", "0.40", "0.65", "0.70", "0.80", "0.90"]
    CCcut_TF_all  = ["0.00", "0.20", "0.40", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95"]
    
    for CCcut_TF in CCcut_TF_all:

        #CCcut_TF = "0.65"
        #CCcut_TF = "0.00"
        #CCNumber = 20
        #TRDNumber = 20
        print("Fit Result from this setting: ")
        print("CCNumber:" + str(CCNumber))
        print("TRDNumber:" + str(TRDNumber))
        print("CCcut_TF:" + str(CCcut_TF))


        #### Load Fit Result
        FitResult_negative = TFile.Open(workpath + "/templatefit/negative/FitResult/FitResult_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + "_CCFree.root", "READ")
        FitResult_positive = TFile.Open(workpath + "/templatefit/positive/FitResult/FitResult_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + "_CCFree.root", "READ")
        # load Fit result absolute numbers
        antiproton_number = np.asarray(FitResult_negative.Get("AntiprotonNumber_cccut_" + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        ccproton_number   = np.asarray(FitResult_negative.Get("CCprotonNumber_cccut_"   + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        electron_number   = np.asarray(FitResult_negative.Get("ElectronNumber_cccut_"   + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        # load Fit result error absolute numbers
        delta_antiproton = np.asarray(FitResult_negative.Get("AntiprotonNumberError_cccut_" + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        delta_ccproton   = np.asarray(FitResult_negative.Get("CCprotonNumberError_cccut_"   + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        # load Fit result relative number
        antiproton_number_per = np.asarray(FitResult_negative.Get("AntiprotonNumber_Relative_cccut_" + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        ccproton_number_per   = np.asarray(FitResult_negative.Get("CCprotonNumber_Relative_cccut_" + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        electron_number_per   = np.asarray(FitResult_negative.Get("ElectronNumber_Relative_cccut_" + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        # load Proton number before template fit
        f_RawProton   = TFile.Open(workpath + "/templatefit/positive/ProtonRawNumber/ProtonRawNumber_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")
        proton_number = np.asarray(f_RawProton.Get("ProtonRawNumber_CCcut_" + CCcut_TF))              
        # load Proton number from template fit
        proton_number_TF = np.asarray(FitResult_positive.Get("AntiprotonNumber_cccut_" + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))
        # load Proton number error from template fit
        p_number_error_TF = np.asarray(FitResult_positive.Get("AntiprotonNumberError_cccut_" + str(CCcut_TF) + "_CCN_" + str(CCNumber) + "_" + "TRDN_" + str(TRDNumber) + "_ISSVersion" + suffix))



        #### Calculate and Plot CCLevel, Signal to Background and CCcut Efficiency.
        ## Create object to hold values in RIGIDITY and SIG EFF.
        CCcutValue_AllRigidity_AllEff                      = np.array([])
        CCProtonNumber_AllRigidity_AllEff                  = np.array([])
        ProtonNumber_AllRigidity_AllEff                    = np.array([])

        CCLevel_uncertainty_AllRigidity_AllEff             = np.array([])
        CCLevelMC_uncertainty_AllRigidity_AllEff           = np.array([])
        CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff  = np.array([])
        CCLevelRatio_uncertainty_AllRigidity_AllEff        = np.array([])
        CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff = np.array([])
        CCLevelRatio_AboveCCcutValue_uncertainty_AllRigidity_AllEff          = np.array([])
        CCLevelRatio_AboveCCcutValue_uncertainty_Reweight_AllRigidity_AllEff = np.array([])

        CCPToProtonRatioMC_AllRigidity_AllEff              = np.array([])
        CCPToProtonRatio_AllRigidity_AllEff                = np.array([])
        CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff = np.array([])  

        SignalToBackgroundRatio_AllRigidity_AllEff         = np.array([])
        NNbinningCenterNow_AllRigidity_AllEff              = np.array([])

        binname_showpoints = binname[-ShowPoints:] # Only show last few points

        for signal_efficiency in signal_efficiency_value:
            print('*************************************************************')
            print("Signal Efficiency is " + str(signal_efficiency))

            ## Create object to hold values in RIGIDITY. 
            CCcutValue_AllRigidity                                             = np.array([])
            CCProtonNumber_AllRigidity                                         = np.array([])
            ProtonNumber_AllRigidity                                           = np.array([])
            CCLevel_uncertainty_AllRigidity                                    = unumpy.uarray(np.array([]), np.array([]))
            CCLevelMC_uncertainty_AllRigidity                                  = unumpy.uarray(np.array([]), np.array([]))
            CCLevelMC_uncertainty_Reweight_AllRigidity                         = unumpy.uarray(np.array([]), np.array([]))
            CCLevel_uncertainty_ISSOverMC_AllRigidity                          = unumpy.uarray(np.array([]), np.array([]))
            CCLevel_uncertainty_ISSOverMC_Reweight_AllRigidity                 = unumpy.uarray(np.array([]), np.array([]))
            CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_AllRigidity          = unumpy.uarray(np.array([]), np.array([]))
            CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_Reweight_AllRigidity = unumpy.uarray(np.array([]), np.array([]))
            CCPToProtonRatioMC_AllRigidity                                     = np.array([])
            CCPToProtonRatio_AllRigidity                                       = np.array([])
            CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity                 = np.array([])
            SignalToBackground_AllRigidity                                     = np.array([])
            NNbinningCenterNow_AllRigidity                                     = np.array([])

            ## Loop on Rigidity bins
            if binmergenumber == 2:
                MergeEndIndex = 20
                FluctuateBinStart = 999
                FluctuateBinEnd   = 999 
                BinList = np.append( np.arange(0, MergeEndIndex, 2), np.arange(MergeEndIndex,32,1))
            elif binmergenumber == 5:
                MergeEndIndex = 25
                ## Option1: No merge for Fluctuated range
                #BinList = np.append( np.arange(0, MergeEndIndex, 5), np.arange(MergeEndIndex, 32, 1))
                ## Option2: With merge for Fluctuated range
                FluctuateBinStart = 25
                FluctuateBinEnd = 26
                BinList = np.arange(0, MergeEndIndex, 5)
                BinList = np.append( BinList, np.arange(FluctuateBinStart, FluctuateBinEnd, 2))
                BinList = np.append( BinList, np.arange(MergeEndIndex+(FluctuateBinEnd-FluctuateBinStart+1),32,1))

            for index in BinList:
                print("")

                BinningNow_pi, NNbinningCenterNow = CCLevel_function.MergeRigidityBins(index, MergeEndIndex, binmergenumber, binname_showpoints, NNbinningCenter, FluctuateBinStart)

                ## Load MC
                #MCName = "B1042_pr.pl1.flux.l1a9.2016000_7.6_all"
                MCName = "B1042_pr.pl1.1800_7.6_all"
                #MCName = "B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all"
                MCccp, MCproton = CCLevel_function.LoadMC_CCPAndProton(index, MergeEndIndex, binmergenumber, workpath, MCName, BinningNow_pi, pattern, NNsuffix, FluctuateBinStart)
                # Fix on TRDLikelihood (TH1D has TRD low/high limit, for example 0.0-1.6...)
                MCccp    = MCccp[np.where( (MCccp[:,1]<1.6) & (MCccp[:,1]>0.0) )[0]]
                MCproton = MCproton[np.where( (MCproton[:,1]<1.6) & (MCproton[:,1]>0.0) )[0]]
                # Fix for 'psa' prescaling for positive rigidity events.
                if MCName == 'B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all':
                    MCproton = np.row_stack((MCproton, MCproton, MCproton, MCproton, MCproton))

                ## Load DATA
                ISSp, ISSn = CCLevel_function.LoadISS_PositiveAndNegative(index, MergeEndIndex, binmergenumber, workpath, suffix, BinningNow_pi, pattern, NNsuffix, FluctuateBinStart)

                ## (1). Calculate cc_cut_value according to different signal efficiency in ISS positive data
                CCcutValue = CCLevel_function.Calculate_CC_cutvalue_fromSigEff(ISSp, signal_efficiency, CCcut_TF)
                ## (2). Calcualte CClevel 
                # (Depends on CCcut_TF only, not using CCcutValue or SigEff)
                ## (2.1) Calculate ISS CClevel
                (cclevel_CCcutTF_uncertainty        , ccproton_number_CCcutTF_uncertainty        , proton_number_CCcutTF_uncertainty, 
                 cclevel_AboveCCcutValue_uncertainty, ccproton_number_AboveCCcutValue_uncertainty, proton_number_AboveCCcutValue_uncertainty) = CCLevel_function.CalculateISSCCLevel(MCccp, CCcutValue, CCcut_TF, ISSn, ISSp, ccproton_number_per, antiproton_number, index, ShowPoints, ccproton_number, electron_number, delta_ccproton, proton_number_TF, p_number_error_TF, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd)
                ## (2.2) Calcualte MC CClevel
                (cclevelMC_CCcutTF_uncertainty        , cclevelMC_CCcutTF_Reweighting_uncertainty        , 
                 MCccp_CCcutTF_uncertainty            , MCproton_CCcutTF_uncertainty                     , MCccp_CCcutTF_Reweighting_uncertainty, MCproton_CCcutTF_Reweighting_uncertainty, 
                 cclevelMC_AboveCCcutValue_uncertainty, cclevelMC_AboveCCcutValue_Reweighting_uncertainty) = CCLevel_function.CalculateMCCCLevel(MCName, MCccp, CCcutValue, MCproton, CCcut_TF, pattern)
                ## (2.3) Calcualte ISS/MC CC level 
                CCLevelRatio_CCcutTF_uncertainty                  = cclevel_CCcutTF_uncertainty / cclevelMC_CCcutTF_uncertainty
                CCLevelRatio_CCcutTF_Reweight_uncertainty         = cclevel_CCcutTF_uncertainty / cclevelMC_CCcutTF_Reweighting_uncertainty
                print(cclevel_AboveCCcutValue_uncertainty)
                print(cclevelMC_AboveCCcutValue_uncertainty)
                CCLevelRatio_AboveCCcutValue_uncertainty          = cclevel_AboveCCcutValue_uncertainty / cclevelMC_AboveCCcutValue_uncertainty
                CCLevelRatio_AboveCCcutValue_Reweight_uncertainty = cclevel_AboveCCcutValue_uncertainty / cclevelMC_AboveCCcutValue_Reweighting_uncertainty 
                print("\033[34mCalculateCCLevel(ISSOverMC): \033[0m")
                print("CCLevelRatio (No Reweighting):" + str(CCLevelRatio_CCcutTF_uncertainty.n))
                print("CCLevelRatio (With Reweighting):" + str(CCLevelRatio_CCcutTF_Reweight_uncertainty.n))
                ## (3). CCProton/ChargeCorrectProton 
                # (Depends on CCcut_TF only, not using CCcutValue or SigEff)
                ccratio_uncertainty, ccratioMC_uncertainty, ccratioRatio_uncertainty = CCLevel_function.Calculate_CCProtonOverProton(ccproton_number_CCcutTF_uncertainty, proton_number_CCcutTF_uncertainty, MCccp_CCcutTF_Reweighting_uncertainty, MCproton_CCcutTF_Reweighting_uncertainty) 
                '''
                ## (4). Calcualte signal to background ratio
                # (Depends on CCcut_TF only, not using CCcutValue or SigEff)
                #signal_number = signal_efficiency * antiproton_number[index+(32-ShowPoints)] 
                signal_number = antiproton_number[index+(32-ShowPoints)]
                signal_to_background = signal_number / ccproton_number_CCcutTF_uncertainty.n
                '''
                ## (5). collect values in different Rigidity bins        
                CCcutValue_AllRigidity                                             = np.append(CCcutValue_AllRigidity                     , CCcutValue)
                CCProtonNumber_AllRigidity                                         = np.append(CCProtonNumber_AllRigidity                 , ccproton_number_CCcutTF_uncertainty.n)
                ProtonNumber_AllRigidity                                           = np.append(ProtonNumber_AllRigidity                   , proton_number_CCcutTF_uncertainty.n)
                # CCLevel
                CCLevel_uncertainty_AllRigidity                                    = np.append(CCLevel_uncertainty_AllRigidity            , cclevel_CCcutTF_uncertainty)
                CCLevelMC_uncertainty_AllRigidity                                  = np.append(CCLevelMC_uncertainty_AllRigidity          , cclevelMC_CCcutTF_uncertainty)
                CCLevelMC_uncertainty_Reweight_AllRigidity                         = np.append(CCLevelMC_uncertainty_Reweight_AllRigidity , cclevelMC_CCcutTF_Reweighting_uncertainty)
                CCLevel_uncertainty_ISSOverMC_AllRigidity                          = np.append(CCLevel_uncertainty_ISSOverMC_AllRigidity                , CCLevelRatio_CCcutTF_uncertainty) 
                CCLevel_uncertainty_ISSOverMC_Reweight_AllRigidity                 = np.append(CCLevel_uncertainty_ISSOverMC_Reweight_AllRigidity       , CCLevelRatio_CCcutTF_Reweight_uncertainty)
                CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_AllRigidity          = np.append(CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_AllRigidity, CCLevelRatio_AboveCCcutValue_uncertainty) 
                CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_Reweight_AllRigidity = np.append(CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_Reweight_AllRigidity, CCLevelRatio_AboveCCcutValue_Reweight_uncertainty)
                # CCPToProtonRatio                
                CCPToProtonRatio_AllRigidity                       = np.append(CCPToProtonRatio_AllRigidity                      , ccratio_uncertainty.n)
                CCPToProtonRatioMC_AllRigidity                     = np.append(CCPToProtonRatioMC_AllRigidity                    , ccratioMC_uncertainty.n)
                CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity = np.append(CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity, ccratioRatio_uncertainty)
                #SignalToBackground_AllRigidity                     = np.append(SignalToBackground_AllRigidity             , signal_to_background)
                NNbinningCenterNow_AllRigidity                     = np.append( NNbinningCenterNow_AllRigidity                    , NNbinningCenterNow)

            ## collect values in all Signal Efficiency
            CCcutValue_AllRigidity_AllEff                      = CCLevel_function.VStack(CCcutValue_AllRigidity_AllEff                     , CCcutValue_AllRigidity)
            CCProtonNumber_AllRigidity_AllEff                  = CCLevel_function.VStack(CCProtonNumber_AllRigidity_AllEff                 , CCProtonNumber_AllRigidity)
            ProtonNumber_AllRigidity_AllEff                    = CCLevel_function.VStack(ProtonNumber_AllRigidity_AllEff                   , ProtonNumber_AllRigidity)
            # CCLevel
            CCLevel_uncertainty_AllRigidity_AllEff             = np.append(CCLevel_uncertainty_AllRigidity_AllEff                          , CCLevel_uncertainty_AllRigidity)
            CCLevelMC_uncertainty_AllRigidity_AllEff           = np.append(CCLevelMC_uncertainty_AllRigidity_AllEff                        , CCLevelMC_uncertainty_AllRigidity)
            CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff  = np.append(CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff               , CCLevelMC_uncertainty_Reweight_AllRigidity) 
            CCLevelRatio_uncertainty_AllRigidity_AllEff                          = np.append(CCLevelRatio_uncertainty_AllRigidity_AllEff                         , CCLevel_uncertainty_ISSOverMC_AllRigidity)
            CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff                 = np.append(CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff                , CCLevel_uncertainty_ISSOverMC_Reweight_AllRigidity)
            CCLevelRatio_AboveCCcutValue_uncertainty_AllRigidity_AllEff          = np.append(CCLevelRatio_AboveCCcutValue_uncertainty_AllRigidity_AllEff         , CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_AllRigidity)
            CCLevelRatio_AboveCCcutValue_uncertainty_Reweight_AllRigidity_AllEff = np.append(CCLevelRatio_AboveCCcutValue_uncertainty_Reweight_AllRigidity_AllEff, CCLevel_AboveCCcutValue_uncertainty_ISSOverMC_Reweight_AllRigidity)
            # CCPToProtonRatio
            CCPToProtonRatioMC_AllRigidity_AllEff                     = CCLevel_function.VStack(CCPToProtonRatioMC_AllRigidity_AllEff             , CCPToProtonRatioMC_AllRigidity)
            CCPToProtonRatio_AllRigidity_AllEff                       = CCLevel_function.VStack(CCPToProtonRatio_AllRigidity_AllEff               , CCPToProtonRatio_AllRigidity)
            CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff = np.append(CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff       , CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity) 
            #SignalToBackgroundRatio_AllRigidity_AllEff                = CCLevel_function.VStack(SignalToBackgroundRatio_AllRigidity_AllEff        , SignalToBackground_AllRigidity)
            NNbinningCenterNow_AllRigidity_AllEff                      = np.append(NNbinningCenterNow_AllRigidity_AllEff                           , NNbinningCenterNow_AllRigidity)

        '''
        #### Plot
        ## (1). Plot Signal To Background_Ratio
        #CCLevel_function.Plot_Signal_To_Background_Ratio(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, signal_efficiency_value, pattern, NNsuffix)
        #CCLevel_function.Plot_CCcut_Efficiency          (resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, CCcutValue_AllRigidity_AllEff, signal_efficiency_value, pattern, NNsuffix)
        ## (2). Plot CCLevel
        #CCLevel_function.Plot_CCLevel                   (resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, unumpy.nominal_values(CCLevel_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCLevel_uncertainty_AllRigidity_AllEff), signal_efficiency_value, binningCenter_Low, cclevelMCLow, cclevelErrorMCLow, cclevelMCLow_allPattern, cclevelErrorMCLow_allPattern, cclevelMCHigh, cclevelErrorMCHigh, unumpy.nominal_values(CCLevelMC_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCLevelMC_uncertainty_AllRigidity_AllEff), ShowPoints, pattern, NNsuffix)
        CCLevel_function.Plot_CCLevel                   (resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, unumpy.nominal_values(CCLevel_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCLevel_uncertainty_AllRigidity_AllEff), signal_efficiency_value, binningCenter_Low, unumpy.nominal_values(CCLevelMC_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCLevelMC_uncertainty_AllRigidity_AllEff), unumpy.nominal_values(CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff), unumpy.std_devs(CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff), ShowPoints, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd)
        CCLevel_function.Plot_CCLevel_ISSOverMC         (resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, unumpy.nominal_values(CCLevelRatio_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCLevelRatio_uncertainty_AllRigidity_AllEff), unumpy.nominal_values(CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff), unumpy.std_devs(CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff), unumpy.nominal_values(CCLevelRatio_AboveCCcutValue_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCLevelRatio_AboveCCcutValue_uncertainty_AllRigidity_AllEff), unumpy.nominal_values(CCLevelRatio_AboveCCcutValue_uncertainty_Reweight_AllRigidity_AllEff), unumpy.std_devs(CCLevelRatio_AboveCCcutValue_uncertainty_Reweight_AllRigidity_AllEff), signal_efficiency_value, ShowPoints, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd)
        ## (3). Plot CCProton and Proton numbers (Depends on CCcut_TF only, not using CCcutValue or SigEff)
        #CCLevel_function.Plot_CCprotonNumbers           (resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, CCProtonNumber_AllRigidity_AllEff, signal_efficiency_value, ShowPoints, pattern, NNsuffix)
        #CCLevel_function.Plot_ProtonNumbers             (resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, ProtonNumber_AllRigidity_AllEff, ShowPoints, signal_efficiency_value, pattern, NNsuffix)
        ## (4). Plot CCProtonToProtonRatio (Depends on CCcut_TF only, not using CCcutValue or SigEff)
        CCLevel_function.Plot_CCProton_Proton_Ratio     (resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, unumpy.nominal_values(CCPToProtonRatio_AllRigidity_AllEff), unumpy.std_devs(CCPToProtonRatio_AllRigidity_AllEff), unumpy.nominal_values(CCPToProtonRatioMC_AllRigidity_AllEff), unumpy.std_devs(CCPToProtonRatioMC_AllRigidity_AllEff), signal_efficiency_value, ShowPoints, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd)
        CCLevel_function.Plot_CCProton_Proton_Ratio_ISSOverMC(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, SignalToBackgroundRatio_AllRigidity_AllEff, NNbinningCenter, unumpy.nominal_values(CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff), ShowPoints, signal_efficiency_value, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd)
        ## (5). Plot Antiproton to Proton Ratio (Depends on CCcut_TF only, not using CCcutValue or SigEff)
        #AntiprotonToProtonRatio = antiproton_number/proton_number_TF
        #CCLevel_function.PlotAntiprotonToProtonRatio(resultpath, ISSversion, Binningversion, NNbinningCenter, AntiprotonToProtonRatio, RigidityBinCenter_450version, CCcut_TF, suffix, pattern, NNsuffix)
        ## (6). Plot Statistical Error 
        #CCLevel_function.Plot_Sta_Error(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, NNbinningCenter, RigidityBinCenter_450version, delta_antiproton, proton_number_TF, AntiprotonToProtonRatio, pattern, NNsuffix)
        ## (7). Plot Proton and Antiproton Number (Only Full Span: 175-211, 211-259, 259-450 or 259-330 and 330-525)
        #CCLevel_function.Plot_Antiproton_Number(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, NNbinningCenter, antiproton_number, pattern, NNsuffix)
        #CCLevel_function.Plot_Proton_number(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, NNbinningCenter, proton_number_TF, pattern, NNsuffix)
        '''
        
        #### Save CCProton_To_Proton_Ratio
        print("signal_efficiency_value:") 
        print(signal_efficiency_value)
        print("CCLevel_uncertainty_AllRigidity_AllEff:")
        print(CCLevel_uncertainty_AllRigidity_AllEff.shape)
        print("CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff:")
        print(CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff.shape)
        f_CCLevel_ISSOverMC.cd()
        ''' 
        f_CCLevel_ISSOverMC.WriteObject(binmergenumber                            , "binmergenumber")
        f_CCLevel_ISSOverMC.WriteObject(NNbinningCenter                           , "NNbinningCenter")
        f_CCLevel_ISSOverMC.WriteObject(MergeEndIndex                             , "MergeEndIndex")
        f_CCLevel_ISSOverMC.WriteObject(FluctuateBinStart                         , "FluctuateBinStart")
        f_CCLevel_ISSOverMC.WriteObject(FluctuateBinEnd                           , "FluctuateBinEnd")

        unumpy.nominal_values(CCLevel_uncertainty_AllRigidity_AllEff), unumpy.std_devs(CCLevel_uncertainty_AllRigidity_AllEff),
        f_CCLevel_ISSOverMC.WriteObject(cclevel_all_diffeff                       , "cclevel_all_diffeff")
        f_CCLevel_ISSOverMC.WriteObject(cclevelerror_all_diffeff                  , "cclevelerror_all_diffeff")

        unumpy.nominal_values(CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff), unumpy.std_devs(CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff),
        f_CCLevel_ISSOverMC.WriteObject(CCLevelMC_Reweight_AllRigidity_AllEff     , "CCLevelMC_Reweight_AllRigidity_AllEff")
        f_CCLevel_ISSOverMC.WriteObject(CCLevelMCError_Reweight_AllRigidity_AllEff, "CCLevelMCError_Reweight_AllRigidity_AllEff")
        '''
        if len(signal_efficiency_value) > 1: 
            for i in range(CCLevel_uncertainty_AllRigidity_AllEff.shape[0]):
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(NNbinningCenterNow_AllRigidity_AllEff[i,:].tolist())                                            , "RigidityBinningCenter")
                # CCProton Over Proton
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCPToProtonRatio_AllRigidity_AllEff[i,:]).tolist()), "CCProtonOverProtonISS_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCPToProtonRatio_AllRigidity_AllEff[i,:]).tolist()), "CCProtonOverProtonISS_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCPToProtonRatioMC_AllRigidity_AllEff[i,:]).tolist()), "CCProtonOverProtonMC_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCPToProtonRatioMC_AllRigidity_AllEff[i,:]).tolist()), "CCProtonOverProtonMC_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff[i,:]).tolist()), "CCProtonOverProtonRatio_ISSOverMC_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff[i,:]).tolist()), "CCProtonOverProtonRatio_ISSOverMC_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                # CC Level Ratio
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCLevel_uncertainty_AllRigidity_AllEff[i,:]).tolist())                   , "CCLevelISS_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCLevel_uncertainty_AllRigidity_AllEff[i,:]).tolist())                   , "CCLevelISS_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )

                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff[i,:]).tolist())        , "CCLevelMC_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff[i,:]).tolist())        , "CCLevellMC_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )

                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff[i,:]).tolist())         , "CLevelRatio_" + "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
                f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff[i,:]).tolist())         , "CLevelRatioError_" + "Signal_Efficiency_" + str(signal_efficiency_value[i]) + "_CCcut_TF_" + str(CCcut_TF) )
        else:
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(NNbinningCenterNow_AllRigidity_AllEff.tolist())                                            , "RigidityBinningCenter")
            # CCProton Over Proton
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCPToProtonRatio_AllRigidity_AllEff).tolist()), "CCProtonOverProtonISS_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCPToProtonRatio_AllRigidity_AllEff).tolist()), "CCProtonOverProtonISS_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCPToProtonRatioMC_AllRigidity_AllEff).tolist()), "CCProtonOverProtonMC_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCPToProtonRatioMC_AllRigidity_AllEff).tolist()), "CCProtonOverProtonMC_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff).tolist()), "CCProtonOverProtonRatio_ISSOverMC_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCPToProtonRatio_ISSOverMC_uncertainty_AllRigidity_AllEff).tolist()), "CCProtonOverProtonRatio_ISSOverMC_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            # CC Level Ratio
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCLevel_uncertainty_AllRigidity_AllEff).tolist())                   , "CCLevelISS_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCLevel_uncertainty_AllRigidity_AllEff).tolist())                   , "CCLevelISS_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )         
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff).tolist())                   , "CCLevelMC_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCLevelMC_uncertainty_Reweight_AllRigidity_AllEff).tolist())                   , "CCLevellMC_Error_"+ "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.nominal_values( CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff).tolist())     , "CLevelRatio_" + "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) ) 
            f_CCLevel_ISSOverMC.WriteObject( CCLevel_function.ListToVector(unumpy.std_devs      ( CCLevelRatio_uncertainty_Reweight_AllRigidity_AllEff).tolist())     , "CLevelRatioError_" + "Signal_Efficiency_" + str(signal_efficiency_value[0]) + "_CCcut_TF_" + str(CCcut_TF) )


if __name__ == '__main__':

    #### Parser Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--issversion', help='ISS data version, pass7.8 or published2016')
    parser.add_argument('--binningversion', help='binningversion, 450version or 525version')
    parser.add_argument('--pattern', help='which tracker pattern you choose')
    parser.add_argument('--ifVGGNN', help='if you want to turn to VGGNN, please turn to Yes.')
    arguments = parser.parse_args()

    if (arguments.issversion):
        ISSversion = arguments.issversion
    else:
        print("You need to choose a ISS data version!")
        os._exit(0)

    if ISSversion == "pass7.8":
        suffix = ""
    elif ISSversion == "published2016":
        suffix = "_May2015"
    elif ISSversion == "PhyRep2021":
        suffix = "_Nov2017"

    if (arguments.binningversion):
        Binningversion = arguments.binningversion
    else:
        print("You need to choose a binning version!")
        os._exit(0)

    if (arguments.pattern):
        pattern = arguments.pattern
    else:
        print("You need to choose a tracker pattern!")
        os._exit(0)

    if (arguments.ifVGGNN):
        IfVGGNN = arguments.ifVGGNN
    else:
        print("You need to choose which CC estimator you want to use!")
        os._exit(0)
    if IfVGGNN == "No":
        NNsuffix = ""
    elif IfVGGNN == "Yes":
        NNsuffix = "_VGG16NN"

    ####  free parameters
    workpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis/data'
    resultpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis'
    highpath = os.getenv('HPCHIGHENERGYDATADIR')
    lowpath = os.getenv('HPCLOWENERGYDIR') + '/totalall/'

    #signal_efficiency_value = [0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0] # for signal efiiciency cut
    #signal_efficiency_value = [1.0]
    #signal_efficiency_value = [0.6, 0.7, 0.9, 1.0]
    signal_efficiency_value = [1.0]

    #ShowPoints = 32 # total: 32 points
    #ShowPoints = 8
    #ShowPoints = 12
    ShowPoints = 32

    binmergenumber = 5

    #### Define Binning
    if Binningversion == "525version":
        NNbinningEdges  = binning.Newbinnings_525_zhili[26:59] # 14.1-525:published2016binnings[26:59]
        resultdir       = "results_525version"
        binname         = binning.bins_525_zhili
        NNbinningCenter = binning.Newbinnings_525_center[26:]
    elif Binningversion == "450version":
        NNbinningEdges  = binning.published2016binnings[26:58] #  14.1-450:published2016binnings[26:58]
        resultdir       = "results_450version"
        binname         = binning.bins_450
        NNbinningCenter = binning.published2016binnings_center[26:]

    RigidityBinCenter_450version = binning.published2016binnings_center


    #### A Naive fix to add MC cclevel data in low rigidity range. (But no CCcut?)
    binnameLow        = binning.binslow
    binningCenter_Low = binning.Newbinnings_525_center
    binningCenter_Low = np.insert(binningCenter_Low, 0, 0.9) # Naive Fix to add bincenter for 0.8-1.0 GV
    binningCenter_Low = binningCenter_Low[0:30]              # 0.8-18.0 GV


    main()


