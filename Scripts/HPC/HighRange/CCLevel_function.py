#!/usr/bin/env python
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
from uncertainties import ufloat
from matplotlib.ticker import ScalarFormatter
import ROOT
import PythonPlotDefaultParameters
import uproot


def FitSettingUsage(pattern):
    if pattern == "0":
        CCNumber = 9
        TRDNumber = 11
        CCcut_TF_all = ["0.65"]
    if pattern == "1":
        CCNumber = 9
        TRDNumber = 11
        CCcut_TF_all = ["0.65"]
    elif pattern == "2":
        CCNumber = 20
        TRDNumber = 12
        CCcut_TF_all = ["0.65"]
    elif pattern == "3":
        CCNumber =  20
        TRDNumber = 20
        CCcut_TF_all = ["0.20"]
    elif pattern == "4":
        CCNumber =  20
        TRDNumber = 20
        CCcut_TF_all = ["0.20"]
    elif pattern == "5":
        CCNumber =  20
        TRDNumber = 20
        CCcut_TF_all = ["0.20"]
    elif pattern == "-1":
        CCNumber =  20
        TRDNumber = 20
        CCcut_TF_all = ["0.20"]
    return CCNumber, TRDNumber, CCcut_TF_all


def LoadCCLevelMCLow(lowpath, binnameLow):
    print('\n')
    print("Calculating CCLevel (MC) in Low Rigidity Now.")
    cclevelMCLow                 = []
    cclevelErrorMCLow            = []
    cclevelMCLow_allPattern      = []
    cclevelErrorMCLow_allPattern = []
    for LowIndex in range(binnameLow.shape[0]):
        MCccpLow_allPattern     = root2array(lowpath + "B1042_pr.pl1.1800_7.6_all_Tree_negative_" + str(binnameLow[LowIndex]) + ".root", "AntiprotonLowEnergyTree").shape[0]
        MCprotonLow_allPattern  = root2array(lowpath + "B1042_pr.pl1.1800_7.6_all_Tree_positive_" + str(binnameLow[LowIndex]) + ".root", "AntiprotonLowEnergyTree").shape[0]
        MCccpLow                = root2array(lowpath + "B1042_pr.pl1.1800_7.6_all_Tree_negative_" + str(binnameLow[LowIndex]) + ".root", "AntiprotonLowEnergyTree", selection="Pattern==0").shape[0]
        MCprotonLow             = root2array(lowpath + "B1042_pr.pl1.1800_7.6_all_Tree_positive_" + str(binnameLow[LowIndex]) + ".root", "AntiprotonLowEnergyTree", selection="Pattern==0").shape[0]

        ## Fix for statistic too low
        if MCccpLow==0:
            MCccpLow = 1
        #print("Rigidity Binning: " + str(binnameLow[LowIndex]))
        #print("MCccpLow_allPattern numbers: " + str(MCccpLow_allPattern))
        #print("MCccpLow numbers:" + str(MCccpLow))

        MCccpLow_allPattern_uncertainty       = unumpy.uarray(MCccpLow_allPattern, np.sqrt(MCccpLow_allPattern))
        MCprotonLow_allPattern_uncertainty    = unumpy.uarray(MCprotonLow_allPattern, np.sqrt(MCprotonLow_allPattern))
        cclevelMCLow_allPattern_uncertainty   = MCccpLow_allPattern_uncertainty / (MCccpLow_allPattern_uncertainty + MCprotonLow_allPattern_uncertainty)
        cclevelMCLow_allPattern.append(cclevelMCLow_allPattern_uncertainty.n)
        cclevelErrorMCLow_allPattern.append(cclevelMCLow_allPattern_uncertainty.s)

        MCccpLow_uncertainty       = unumpy.uarray(MCccpLow, np.sqrt(MCccpLow))
        MCprotonLow_uncertainty    = unumpy.uarray(MCprotonLow, np.sqrt(MCprotonLow))
        cclevelMCLow_uncertainty   = MCccpLow_uncertainty / (MCccpLow_uncertainty + MCprotonLow_uncertainty)
        cclevelMCLow.append(cclevelMCLow_uncertainty.n)
        cclevelErrorMCLow.append(cclevelMCLow_uncertainty.s)

    cclevelMCLow      = np.array(cclevelMCLow)
    cclevelErrorMCLow = np.array(cclevelErrorMCLow)
    cclevelMCLow_allPattern      = np.array(cclevelMCLow_allPattern)
    cclevelErrorMCLow_allPattern = np.array(cclevelErrorMCLow_allPattern)
    return cclevelMCLow, cclevelErrorMCLow, cclevelMCLow_allPattern, cclevelErrorMCLow_allPattern


def LoadCCLevelMCHigh(highpath, binname):
    print('\n')
    print("Calculate CCLevel (MC) in High Rigidity: ")
    #### A Naive fix to add MC cclevel data in high rigidity range.
    cclevelMCHigh                 = []
    cclevelErrorMCHigh            = []

    for HighIndex in range(binname.shape[0]):
        MCccpHigh                = root2array(highpath + "/B1042_pr.pl1.1800_7.6_all_Tree_negative_" + str(binname[HighIndex]) + ".root", "ExampleAnalysisTree").shape[0]
        MCprotonLowHigh          = root2array(highpath + "/B1042_pr.pl1.1800_7.6_all_Tree_positive_" + str(binname[HighIndex]) + ".root", "ExampleAnalysisTree").shape[0]
        #MCccpHigh                = root2array(highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_" + str(binname[HighIndex]) + ".root", "ExampleAnalysisTree").shape[0]
        #MCprotonLowHigh          = root2array(highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_" + str(binname[HighIndex]) + ".root", "ExampleAnalysisTree").shape[0]

        ## Fix for statistic too low
        if MCccpHigh == 0:
            MCccpHigh = 1

        MCccpHigh_uncertainty       = unumpy.uarray(MCccpHigh, np.sqrt(MCccpHigh))
        MCprotonHigh_uncertainty    = unumpy.uarray(MCprotonLowHigh, np.sqrt(MCprotonLowHigh))
        cclevelMCHigh_uncertainty   = MCccpHigh_uncertainty / (MCccpHigh_uncertainty + MCprotonHigh_uncertainty)
        cclevelMCHigh.append(cclevelMCHigh_uncertainty.n)
        cclevelErrorMCHigh.append(cclevelMCHigh_uncertainty.s)
    cclevelMCHigh      = np.array(cclevelMCHigh)
    cclevelErrorMCHigh = np.array(cclevelErrorMCHigh)
    return cclevelMCHigh, cclevelErrorMCHigh


def MergeRigidityBins(index, MergeEndIndex, binmergenumber, binname_showpoints, NNbinningCenter, FluctuateBinStart): 
    if index < MergeEndIndex:
        if binmergenumber == 2:
            BinningNow_p1 = binname_showpoints[index]
            BinningNow_p2 = binname_showpoints[index+1]
            BinningNow_pi = [BinningNow_p1, BinningNow_p2]
            NNbinningCenterNow = (NNbinningCenter[index] + NNbinningCenter[index+1]) / 2
            print("\033[31mRigidiy Bin (part1) :\033[0m" + str(BinningNow_p1))
            print("\033[31mRigidiy Bin (part2) :\033[0m" + str(BinningNow_p2))
            print("\033[31mIndex (part1) is:\033[0m" + str(index) )
            print("\033[31mIndex (part2) is:\033[0m" + str(index+1) )
            print("\033[31m31mRigidiy Binning Center is:\033[0m" + str(NNbinningCenterNow) )
        elif binmergenumber == 5:
            BinningNow_p1 = binname_showpoints[index]
            BinningNow_p2 = binname_showpoints[index+1]
            BinningNow_p3 = binname_showpoints[index+2]
            BinningNow_p4 = binname_showpoints[index+3]
            BinningNow_p5 = binname_showpoints[index+4]
            BinningNow_pi = [BinningNow_p1, BinningNow_p2, BinningNow_p3, BinningNow_p4, BinningNow_p5]
            NNbinningCenterNow = ( NNbinningCenter[index] + NNbinningCenter[index+1] + NNbinningCenter[index+2] + NNbinningCenter[index+3] + NNbinningCenter[index+4] ) /5
            print("\033[31mRigidiy Bin (part1) :\033[0m" + str(BinningNow_p1))
            print("\033[31mRigidiy Bin (part2) :\033[0m" + str(BinningNow_p2))
            print("\033[31mRigidiy Bin (part3) :\033[0m" + str(BinningNow_p3))
            print("\033[31mRigidiy Bin (part4) :\033[0m" + str(BinningNow_p4))
            print("\033[31mRigidiy Bin (part5) :\033[0m" + str(BinningNow_p5))
            print("\033[31mIndex (part1) is:\033[0m" + str(index) )
            print("\033[31mIndex (part2) is:\033[0m" + str(index+1) )
            print("\033[31mIndex (part3) is:\033[0m" + str(index+2) )
            print("\033[31mIndex (part4) is:\033[0m" + str(index+3) )
            print("\033[31mIndex (part5) is:\033[0m" + str(index+4) )
            print("\033[31m31mRigidiy Binning Center is:\033[0m" + str(NNbinningCenterNow) )
    elif index == FluctuateBinStart:
            BinningNow_p1 = binname_showpoints[index]
            BinningNow_p2 = binname_showpoints[index+1]
            BinningNow_pi = [BinningNow_p1, BinningNow_p2]
            NNbinningCenterNow = (NNbinningCenter[index] + NNbinningCenter[index+1]) / 2
            print("\033[31mRigidiy Bin (part1) :\033[0m" + str(BinningNow_p1))
            print("\033[31mRigidiy Bin (part2) :\033[0m" + str(BinningNow_p2))
            print("\033[31mIndex (part1) is:\033[0m" + str(index) )
            print("\033[31mIndex (part2) is:\033[0m" + str(index+1) )
            print("\033[31m31mRigidiy Binning Center is:\033[0m" + str(NNbinningCenterNow) )
    else:
        BinningNow    = binname_showpoints[index]
        BinningNow_pi = [BinningNow]
        NNbinningCenterNow = NNbinningCenter[index]
        print("\033[31mRigidiy Bin:\033[0m" + str(BinningNow))
        print("\033[31mIndex is:\033[0m" + str(index) )
        print("\033[31m31mRigidiy Binning Center is:\033[0m" + str(NNbinningCenterNow) )

    return BinningNow_pi, NNbinningCenterNow


def LoadMC_CCPAndProton(index, MergeEndIndex, binmergenumber, workpath, MCName, BinningNow_pi, pattern, NNsuffix, FluctuateBinStart):

    if index < MergeEndIndex:
        if binmergenumber == 2:
            MCccp_p1    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p1 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp_p2    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p2 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp       = np.row_stack((MCccp_p1   , MCccp_p2))
            MCproton    = np.row_stack((MCproton_p1, MCproton_p2))
        elif binmergenumber == 5:
            MCccp_p1    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p1 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp_p2    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p2 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp_p3    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[2]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p3 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[2]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp_p4    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[3]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p4 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[3]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp_p5    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[4]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p5 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[4]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp       = np.row_stack((MCccp_p1   , MCccp_p2   , MCccp_p3   , MCccp_p4   , MCccp_p5))
            MCproton    = np.row_stack((MCproton_p1, MCproton_p2, MCproton_p3, MCproton_p4, MCproton_p5))
    elif index == FluctuateBinStart:
            MCccp_p1    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p1 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp_p2    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCproton_p2 = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            MCccp       = np.row_stack((MCccp_p1   , MCccp_p2))
            MCproton    = np.row_stack((MCproton_p1, MCproton_p2))
    else:
        MCccp    = np.load(workpath + "/ChargeConfusedProtomTemplate_MC_" + MCName + "_Tree_negative_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
        MCproton = np.load(workpath + "/ChargeCorrectProtonTemplate_MC_"  + MCName + "_Tree_positive_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")

    return MCccp, MCproton


def LoadISS_PositiveAndNegative(index, MergeEndIndex, binmergenumber, workpath, suffix, BinningNow_pi, pattern, NNsuffix, FluctuateBinStart):
    if index < MergeEndIndex:
        if binmergenumber == 2:
            ISSp_p1  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p1  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp_p2  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p2  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp     = np.row_stack((ISSp_p1, ISSp_p2))
            ISSn     = np.row_stack((ISSn_p1, ISSn_p2))
        elif binmergenumber == 5:
            ISSp_p1  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p1  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp_p2  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p2  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp_p3  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[2]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p3  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[2]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp_p4  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[3]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p4  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[3]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp_p5  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[4]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p5  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[4]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp     = np.row_stack((ISSp_p1, ISSp_p2, ISSp_p3, ISSp_p4, ISSp_p5))
            ISSn     = np.row_stack((ISSn_p1, ISSn_p2, ISSn_p3, ISSn_p4, ISSn_p5))
    elif index == FluctuateBinStart:
            ISSp_p1  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p1  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp_p2  = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSn_p2  = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[1]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
            ISSp     = np.row_stack((ISSp_p1, ISSp_p2))
            ISSn     = np.row_stack((ISSn_p1, ISSn_p2))
    else:
        ISSp     = np.load(workpath + "/ISS_positive" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")
        ISSn     = np.load(workpath + "/ISS_negative" + suffix +"_" + str(BinningNow_pi[0]) + "_Pattern_" + str(pattern) + NNsuffix + ".npy")

    return ISSp, ISSn



#### Calculate
def Calculate_CC_cutvalue_fromSigEff(ISSp, signal_efficiency, CCcut_TF):
    print("\033[34mCalculate_CC_cutvalue_fromSigEff: \033[0m")
    numbers_corresSigff = np.int(ISSp.shape[0]*signal_efficiency)
    ISSp_sorted         = sorted(ISSp[:,0])
    cc_cut_value        = ISSp_sorted[-numbers_corresSigff]
    print("CCcut value: " + str(cc_cut_value))
    print("CCcut_TF: " + str(CCcut_TF))
    return cc_cut_value


def CalculateISSCCLevel(MCccp, cc_cut_value, CCcut_TF, ISSn, ISSp, ccproton_number_per, antiproton_number, index, ShowPoints, ccproton_number, electron_number, delta_ccproton, proton_number_TF, p_number_error_TF, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd):
    print("\033[34mCalculateISSCCLevel: \033[0m")

    #### Compare CCproton with "CCcut Value from SigEff" and "CCcut used for tempalte fit"
    MCccp_abovecut    = MCccp[np.where(MCccp[:,0]>float(cc_cut_value))[0]].shape[0]
    MCccp_aboveTFCut  = MCccp[np.where(MCccp[:,0]>float(CCcut_TF))[0]].shape[0]
    #Percentage_MCCCcut_Value = MCccp_abovecut / MCccp[:,0].shape[0]
    #Percentage_MCCCcut_TF = MCccp_aboveTFCut / MCccp[:,0].shape[0]
    #print("ISS: CCproton_abovecut/CCproton_aboveTFCut: "+str(MCccp_abovecut/MCccp_aboveTFCut))
    
    # Simple Fix for low statistics MC
    if MCccp_abovecut == 0:
        MCccp_abovecut = 1
    if MCccp_aboveTFCut == 0:
        MCccp_aboveTFCut =1


    #### Calculate CCProton number number 
    #(1). CCcutTF
    if index< MergeEndIndex:
        if binmergenumber == 2:
            background_number_CCcutTF_p1        = ccproton_number[index+(32-ShowPoints)]
            background_number_error_CCcutTF_p1  = delta_ccproton[index+(32-ShowPoints)]
            background_number_CCcutTF_p1_uncertainty = ufloat(background_number_CCcutTF_p1, background_number_error_CCcutTF_p1)

            background_number_CCcutTF_p2        = ccproton_number[index+(32-ShowPoints)+1]
            background_number_error_CCcutTF_p2  = delta_ccproton[index+(32-ShowPoints)+1]
            background_number_CCcutTF_p2_uncertainty = ufloat(background_number_CCcutTF_p2, background_number_error_CCcutTF_p2)

            background_number_CCcutTF           = (background_number_CCcutTF_p1_uncertainty+background_number_CCcutTF_p2_uncertainty).n
            background_number_error_CCcutTF     = (background_number_CCcutTF_p1_uncertainty+background_number_CCcutTF_p2_uncertainty).s
            ccproton_number_CCcutTF_uncertainty = ufloat(background_number_CCcutTF_p1, background_number_error_CCcutTF_p1) + ufloat(background_number_CCcutTF_p2, background_number_error_CCcutTF_p2)
        elif binmergenumber == 5:
            background_number_CCcutTF_p1        = ccproton_number[index+(32-ShowPoints)]
            background_number_error_CCcutTF_p1  = delta_ccproton[index+(32-ShowPoints)]
            background_number_CCcutTF_p1_uncertainty = ufloat(background_number_CCcutTF_p1, background_number_error_CCcutTF_p1)

            background_number_CCcutTF_p2        = ccproton_number[index+(32-ShowPoints)+1]
            background_number_error_CCcutTF_p2  = delta_ccproton[index+(32-ShowPoints)+1]
            background_number_CCcutTF_p2_uncertainty = ufloat(background_number_CCcutTF_p2, background_number_error_CCcutTF_p2)

            background_number_CCcutTF_p3        = ccproton_number[index+(32-ShowPoints)+2]
            background_number_error_CCcutTF_p3  = delta_ccproton[index+(32-ShowPoints)+2]
            background_number_CCcutTF_p3_uncertainty = ufloat(background_number_CCcutTF_p3, background_number_error_CCcutTF_p3)

            background_number_CCcutTF_p4        = ccproton_number[index+(32-ShowPoints)+3]
            background_number_error_CCcutTF_p4  = delta_ccproton[index+(32-ShowPoints)+3]
            background_number_CCcutTF_p4_uncertainty = ufloat(background_number_CCcutTF_p4, background_number_error_CCcutTF_p4)

            background_number_CCcutTF_p5        = ccproton_number[index+(32-ShowPoints)+4]
            background_number_error_CCcutTF_p5  = delta_ccproton[index+(32-ShowPoints)+4]
            background_number_CCcutTF_p5_uncertainty = ufloat(background_number_CCcutTF_p5, background_number_error_CCcutTF_p5)

            background_number_CCcutTF           = (background_number_CCcutTF_p1_uncertainty + background_number_CCcutTF_p2_uncertainty + background_number_CCcutTF_p3_uncertainty + background_number_CCcutTF_p4_uncertainty + background_number_CCcutTF_p5_uncertainty).n
            background_number_error_CCcutTF     = (background_number_CCcutTF_p1_uncertainty + background_number_CCcutTF_p2_uncertainty + background_number_CCcutTF_p3_uncertainty + background_number_CCcutTF_p4_uncertainty + background_number_CCcutTF_p5_uncertainty).s
            ccproton_number_CCcutTF_uncertainty = ufloat(background_number_CCcutTF_p1, background_number_error_CCcutTF_p1) + ufloat(background_number_CCcutTF_p2, background_number_error_CCcutTF_p2) + ufloat(background_number_CCcutTF_p3, background_number_error_CCcutTF_p3) + ufloat(background_number_CCcutTF_p4, background_number_error_CCcutTF_p4) + ufloat(background_number_CCcutTF_p5, background_number_error_CCcutTF_p5)
    elif index == FluctuateBinStart:
            background_number_CCcutTF_p1        = ccproton_number[index+(32-ShowPoints)]
            background_number_error_CCcutTF_p1  = delta_ccproton[index+(32-ShowPoints)]
            background_number_CCcutTF_p1_uncertainty = ufloat(background_number_CCcutTF_p1, background_number_error_CCcutTF_p1)

            background_number_CCcutTF_p2        = ccproton_number[index+(32-ShowPoints)+1]
            background_number_error_CCcutTF_p2  = delta_ccproton[index+(32-ShowPoints)+1]
            background_number_CCcutTF_p2_uncertainty = ufloat(background_number_CCcutTF_p2, background_number_error_CCcutTF_p2)

            background_number_CCcutTF           = (background_number_CCcutTF_p1_uncertainty+background_number_CCcutTF_p2_uncertainty).n
            background_number_error_CCcutTF     = (background_number_CCcutTF_p1_uncertainty+background_number_CCcutTF_p2_uncertainty).s
            ccproton_number_CCcutTF_uncertainty = ufloat(background_number_CCcutTF_p1, background_number_error_CCcutTF_p1) + ufloat(background_number_CCcutTF_p2, background_number_error_CCcutTF_p2)
    else:
        background_number_CCcutTF       = ccproton_number[index+(32-ShowPoints)] 
        background_number_error_CCcutTF = delta_ccproton[index+(32-ShowPoints)]
        ccproton_number_CCcutTF_uncertainty = ufloat(background_number_CCcutTF, np.sqrt(background_number_CCcutTF))
        #ccproton_number_CCcutTF_uncertainty = ufloat(background_number_CCcutTF, background_number_error_CCcutTF)
    print("CCProton Error check: background_number_error: " + str(background_number_error_CCcutTF) + ", square root error is " + str(np.sqrt(background_number_CCcutTF)) )
    #print("ISS: ccproton_number_CCcutTF_uncertainty: " + str(ccproton_number_CCcutTF_uncertainty))
    #(2). CCcutValue
    if index<MergeEndIndex:
        if binmergenumber == 2:
            ccproton_number_AboveCCcutValue_p1             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)]
            ccproton_number_error_AboveCCcutValue_p1       = delta_ccproton[index+(32-ShowPoints)] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p1
            ccproton_number_AboveCCcutValue_p1_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p1, ccproton_number_error_AboveCCcutValue_p1) 

            ccproton_number_AboveCCcutValue_p2             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)+1]
            ccproton_number_error_AboveCCcutValue_p2       = delta_ccproton[index+(32-ShowPoints)+1] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p2
            ccproton_number_AboveCCcutValue_p2_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p2, ccproton_number_error_AboveCCcutValue_p2)

            ccproton_number_AboveCCcutValue       = (ccproton_number_AboveCCcutValue_p1_uncertainty + ccproton_number_AboveCCcutValue_p2_uncertainty).n
            ccproton_number_error_AboveCCcutValue = (ccproton_number_AboveCCcutValue_p1_uncertainty + ccproton_number_AboveCCcutValue_p2_uncertainty).s
            ccproton_number_AboveCCcutValue_uncertainty = ufloat(ccproton_number_AboveCCcutValue, ccproton_number_error_AboveCCcutValue)
        elif binmergenumber == 5:
            ccproton_number_AboveCCcutValue_p1             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)]
            ccproton_number_error_AboveCCcutValue_p1       = delta_ccproton[index+(32-ShowPoints)] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p1
            ccproton_number_AboveCCcutValue_p1_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p1, ccproton_number_error_AboveCCcutValue_p1)

            ccproton_number_AboveCCcutValue_p2             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)+1]
            ccproton_number_error_AboveCCcutValue_p2       = delta_ccproton[index+(32-ShowPoints)+1] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p2
            ccproton_number_AboveCCcutValue_p2_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p2, ccproton_number_error_AboveCCcutValue_p2)

            ccproton_number_AboveCCcutValue_p3             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)+2]
            ccproton_number_error_AboveCCcutValue_p3       = delta_ccproton[index+(32-ShowPoints)+2] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p3
            ccproton_number_AboveCCcutValue_p3_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p3, ccproton_number_error_AboveCCcutValue_p3)

            ccproton_number_AboveCCcutValue_p4             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)+3]
            ccproton_number_error_AboveCCcutValue_p4       = delta_ccproton[index+(32-ShowPoints)+3] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p4
            ccproton_number_AboveCCcutValue_p4_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p4, ccproton_number_error_AboveCCcutValue_p4)

            ccproton_number_AboveCCcutValue_p5             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)+4]
            ccproton_number_error_AboveCCcutValue_p5       = delta_ccproton[index+(32-ShowPoints)+4] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p5
            ccproton_number_AboveCCcutValue_p5_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p5, ccproton_number_error_AboveCCcutValue_p5)

            ccproton_number_AboveCCcutValue       = (ccproton_number_AboveCCcutValue_p1_uncertainty + ccproton_number_AboveCCcutValue_p2_uncertainty + ccproton_number_AboveCCcutValue_p3_uncertainty + ccproton_number_AboveCCcutValue_p4_uncertainty + ccproton_number_AboveCCcutValue_p5_uncertainty).n
            ccproton_number_error_AboveCCcutValue = (ccproton_number_AboveCCcutValue_p1_uncertainty + ccproton_number_AboveCCcutValue_p2_uncertainty + ccproton_number_AboveCCcutValue_p3_uncertainty + ccproton_number_AboveCCcutValue_p4_uncertainty + ccproton_number_AboveCCcutValue_p5_uncertainty).s
            ccproton_number_AboveCCcutValue_uncertainty = ufloat(ccproton_number_AboveCCcutValue, ccproton_number_error_AboveCCcutValue)
    elif index == FluctuateBinStart:
            ccproton_number_AboveCCcutValue_p1             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)]
            ccproton_number_error_AboveCCcutValue_p1       = delta_ccproton[index+(32-ShowPoints)] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p1
            ccproton_number_AboveCCcutValue_p1_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p1, ccproton_number_error_AboveCCcutValue_p1)

            ccproton_number_AboveCCcutValue_p2             = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)+1]
            ccproton_number_error_AboveCCcutValue_p2       = delta_ccproton[index+(32-ShowPoints)+1] / background_number_CCcutTF * ccproton_number_AboveCCcutValue_p2
            ccproton_number_AboveCCcutValue_p2_uncertainty = ufloat(ccproton_number_AboveCCcutValue_p2, ccproton_number_error_AboveCCcutValue_p2)

            ccproton_number_AboveCCcutValue       = (ccproton_number_AboveCCcutValue_p1_uncertainty + ccproton_number_AboveCCcutValue_p2_uncertainty).n
            ccproton_number_error_AboveCCcutValue = (ccproton_number_AboveCCcutValue_p1_uncertainty + ccproton_number_AboveCCcutValue_p2_uncertainty).s
            ccproton_number_AboveCCcutValue_uncertainty = ufloat(ccproton_number_AboveCCcutValue, ccproton_number_error_AboveCCcutValue)
    else:
        #ccproton_number_AboveCCcutValue = MCccp_abovecut/MCccp.shape[0] * ISSn.shape[0]                                                     # Set1: only works in high energy ragne, CCproton dominent.
        ccproton_number_AboveCCcutValue = MCccp_abovecut/MCccp_aboveTFCut * ccproton_number[index+(32-ShowPoints)]                          # Set2: CCproton_Number_fromTemplateFit
        #ccproton_number_AboveCCcutValue = MCccp_abovecut/MCccp.shape[0] * ISSn.shape[0] * ccproton_number_per[index+(32-ShowPoints)]         # Set3: ISSn.shape[0] * CC_Relative_Ratio
        #ccproton_number_AboveCCcutValue = MCccp_abovecut/MCccp.shape[0] * ISSn.shape[0] * (1-antiproton_number_per[index+(32-ShowPoints)])  # Set4.1: ISSn.shape[0] * (CC_Relative_Ratio + Electron_Relative_Ratio)
        #ccproton_number_AboveCCcutValue = MCccp_abovecut/MCccp.shape[0] * (antiproton_number[index+(32-ShowPoints)] + ccproton_number[index+(32-ShowPoints)] + electron_number[index+(32-ShowPoints)])  * (1-antiproton_number_per[index+(32-ShowPoints)])                                                                                              # Set4.2: TotalTFNumber * (CC_Relative_Ratio + Electron_Relative_Ratio)
        #ccproton_number_AboveCCcutValue = MCccp_abovecut/MCccp.shape[0] * ISSn.shape[0] * (1-electron_number_per[index+(32-ShowPoints)])    # Set5: ISSn.shape[0] * (CC_Relative_Ratio + Pbar_Relative_Ratio)
        #ccproton_number_error_AboveCCcutValue = np.sqrt(background_number_CCcutTF) / background_number_CCcutTF * ccproton_number_AboveCCcutValue
        ccproton_number_error_AboveCCcutValue = delta_ccproton[index+(32-ShowPoints)] / background_number_CCcutTF * ccproton_number_AboveCCcutValue
        ccproton_number_AboveCCcutValue_uncertainty = ufloat(ccproton_number_AboveCCcutValue, ccproton_number_error_AboveCCcutValue)

    #### Calcualte ISS Proton number
    #(1). CCcutTF
    if index<MergeEndIndex:
        if binmergenumber == 2:
            protonsignal_above_TF_p1             = proton_number_TF[index+(32-ShowPoints)]
            proton_TF_error_p1                   = p_number_error_TF[index+(32-ShowPoints)]
            protonsignal_above_TF_p1_uncertainty = ufloat(protonsignal_above_TF_p1, proton_TF_error_p1)

            protonsignal_above_TF_p2             = proton_number_TF[index+(32-ShowPoints)+1]
            proton_TF_error_p2                   = p_number_error_TF[index+(32-ShowPoints)+1]
            protonsignal_above_TF_p2_uncertainty = ufloat(protonsignal_above_TF_p2, proton_TF_error_p2)

            protonsignal_above_TF             = (protonsignal_above_TF_p1_uncertainty + protonsignal_above_TF_p2_uncertainty).n
            proton_TF_error                   = (protonsignal_above_TF_p1_uncertainty + protonsignal_above_TF_p2_uncertainty).s
            proton_number_CCcutTF_uncertainty = ufloat(protonsignal_above_TF, np.sqrt(protonsignal_above_TF))
        elif binmergenumber == 5:
            protonsignal_above_TF_p1             = proton_number_TF[index+(32-ShowPoints)]
            proton_TF_error_p1                   = p_number_error_TF[index+(32-ShowPoints)]
            protonsignal_above_TF_p1_uncertainty = ufloat(protonsignal_above_TF_p1, proton_TF_error_p1)

            protonsignal_above_TF_p2             = proton_number_TF[index+(32-ShowPoints)+1]
            proton_TF_error_p2                   = p_number_error_TF[index+(32-ShowPoints)+1]
            protonsignal_above_TF_p2_uncertainty = ufloat(protonsignal_above_TF_p2, proton_TF_error_p2)

            protonsignal_above_TF_p3             = proton_number_TF[index+(32-ShowPoints)+2]
            proton_TF_error_p3                   = p_number_error_TF[index+(32-ShowPoints)+2]
            protonsignal_above_TF_p3_uncertainty = ufloat(protonsignal_above_TF_p3, proton_TF_error_p3)

            protonsignal_above_TF_p4             = proton_number_TF[index+(32-ShowPoints)+3]
            proton_TF_error_p4                   = p_number_error_TF[index+(32-ShowPoints)+3]
            protonsignal_above_TF_p4_uncertainty = ufloat(protonsignal_above_TF_p4, proton_TF_error_p4)

            protonsignal_above_TF_p5             = proton_number_TF[index+(32-ShowPoints)+4]
            proton_TF_error_p5                   = p_number_error_TF[index+(32-ShowPoints)+4]
            protonsignal_above_TF_p5_uncertainty = ufloat(protonsignal_above_TF_p5, proton_TF_error_p5)

            protonsignal_above_TF             = (protonsignal_above_TF_p1_uncertainty + protonsignal_above_TF_p2_uncertainty + protonsignal_above_TF_p3_uncertainty + protonsignal_above_TF_p4_uncertainty + protonsignal_above_TF_p5_uncertainty).n
            proton_TF_error                   = (protonsignal_above_TF_p1_uncertainty + protonsignal_above_TF_p2_uncertainty + protonsignal_above_TF_p3_uncertainty + protonsignal_above_TF_p4_uncertainty + protonsignal_above_TF_p5_uncertainty).s
            proton_number_CCcutTF_uncertainty = ufloat(protonsignal_above_TF, np.sqrt(protonsignal_above_TF))
    elif index == FluctuateBinStart:
            protonsignal_above_TF_p1             = proton_number_TF[index+(32-ShowPoints)]
            proton_TF_error_p1                   = p_number_error_TF[index+(32-ShowPoints)]
            protonsignal_above_TF_p1_uncertainty = ufloat(protonsignal_above_TF_p1, proton_TF_error_p1)

            protonsignal_above_TF_p2             = proton_number_TF[index+(32-ShowPoints)+1]
            proton_TF_error_p2                   = p_number_error_TF[index+(32-ShowPoints)+1]
            protonsignal_above_TF_p2_uncertainty = ufloat(protonsignal_above_TF_p2, proton_TF_error_p2)

            protonsignal_above_TF             = (protonsignal_above_TF_p1_uncertainty + protonsignal_above_TF_p2_uncertainty).n
            proton_TF_error                   = (protonsignal_above_TF_p1_uncertainty + protonsignal_above_TF_p2_uncertainty).s
            proton_number_CCcutTF_uncertainty = ufloat(protonsignal_above_TF, np.sqrt(protonsignal_above_TF))
    else:
        protonsignal_above_TF             = proton_number_TF[index+(32-ShowPoints)]
        proton_TF_error                   = p_number_error_TF[index+(32-ShowPoints)]
        proton_number_CCcutTF_uncertainty = ufloat(protonsignal_above_TF, np.sqrt(protonsignal_above_TF))
    print("Proton Error check: proton_TF_error is " + str(proton_TF_error) + ", square root error is " + str(np.sqrt(protonsignal_above_TF)) )
    print("ISS: proton_number_CCcutTF_uncertainty: " + str(proton_number_CCcutTF_uncertainty))
    #(2). CCcutValue (it is using ISSp to calculate, first 20 points and last 12 points are using same way)
    proton_number_AboveCCcutValue_orig = ISSp[np.where(ISSp[:,0]>cc_cut_value)[0],0].shape[0]
    #Percentage_AboveCC     = np.where(ISSp[:,0]>cc_cut_value)[0].shape[0] / np.where(ISSp[:,0]>float(CCcut_TF))[0].shape[0]
    #protonsignal_above_TF   = protonsignal_above_TF * Percentage_AboveCC
    #protonsignal_above_TF   = protonsignal_above_TF * signal_efficiency
    #print("ISS: protonsignal_above_TF/proton_number_AboveCCcutValue_orig:" + str(protonsignal_above_TF/proton_number_AboveCCcutValue_orig))
    #proton_number_error_AboveCCcutValue = p_number_error_TF[index+(32-ShowPoints)] * (MCccp_abovecut/MCccp_aboveTFCut)
    proton_number_error_AboveCCcutValue = np.sqrt(proton_number_AboveCCcutValue_orig)
    proton_number_AboveCCcutValue_uncertainty = unumpy.uarray(proton_number_AboveCCcutValue_orig, proton_number_error_AboveCCcutValue)

    #### Calcualte ISS cclevel
    #(1). CCcutTF
    cclevel_uncertainty = ccproton_number_CCcutTF_uncertainty / (ccproton_number_CCcutTF_uncertainty + proton_number_CCcutTF_uncertainty)
    #print("ISS CCLevel Uncertainty:" + str(cclevel_uncertainty))
    #print("ISS: RelativeError:" + str(cclevel_uncertainty.s/cclevel_uncertainty.n))
    #(2). CCcutValue
    cclevel_AboveCCcutValue_uncertainty = ccproton_number_AboveCCcutValue_uncertainty / (ccproton_number_AboveCCcutValue_uncertainty + proton_number_AboveCCcutValue_uncertainty) 

    return cclevel_uncertainty, ccproton_number_CCcutTF_uncertainty, proton_number_CCcutTF_uncertainty, cclevel_AboveCCcutValue_uncertainty, ccproton_number_AboveCCcutValue_uncertainty, proton_number_AboveCCcutValue_uncertainty


def CalculateMCCCLevel(MCName, MCccp, cc_cut_value, MCproton, CCcut_TF, pattern):

    print("\033[34mCalculateMCCCLevel: \033[0m")
     
    #### Update Rewights for MC 
    ReweightingFile = uproot.open(os.getenv('HPCHIGHENERGYDATADIR') + "/ISS_anylsis/antiproton_to_proton_ratio_plot/" + "ReweightingFile_" + MCName + "_Pattern_" + str(pattern) + ".root")
    xaxis = ReweightingFile['ReweightingFactor_More'].axis().edges()
    yvalue = ReweightingFile['ReweightingFactor_More'].values()

    for index, value in enumerate(MCccp[:,3]):
        R_true = MCccp[:,3][index]
        for p in range(xaxis.shape[0]-1):
            if (xaxis[p]<R_true) and (R_true<xaxis[p+1]):
                MCccp[index,2] = yvalue[p]             
                break       
            if p == xaxis.shape[0]-1-1:
                print("A MCccp Anomaly R_true" + str(R_true) + " GV.") 

    for index, value in enumerate(MCproton[:,3]):
        R_true = MCproton[:,3][index]
        for p in range(xaxis.shape[0]-1):
            if (xaxis[p]<R_true) and (R_true<xaxis[p+1]):
                MCproton[index,2] = yvalue[p]
                break
            if p == xaxis.shape[0]-1-1:
                print("A MCproton Anomaly R_true" + str(R_true) + " GV.")

    #### Calculate CCProton number  
    #(1). CCcutTF
    MCccp_abovecut = MCccp[np.where(MCccp[:,0]>float(CCcut_TF))[0]]
    MCccp_abovecut_value = MCccp_abovecut.shape[0]
    # Simple Fix for low statistics MC
    if MCccp_abovecut_value == 0:
        MCccp_abovecut_value = 1
    MCccp_CCcutTF_uncertainty = ufloat(MCccp_abovecut_value, np.sqrt(MCccp_abovecut_value))
    #print("MC: MCccp_CCcutTF_uncertainty: " + str(MCccp_CCcutTF_uncertainty))
    #(2). CCcutValue
    MCccp_AboveCCcutValue = MCccp[np.where(MCccp[:,0]>cc_cut_value)[0]]
    MCccp_AboveCCcutValue_value = MCccp_AboveCCcutValue.shape[0]
    if MCccp_AboveCCcutValue_value == 0:
        MCccp_AboveCCcutValue_value = 1
    MCccp_AboveCCcutValue_uncertainty = ufloat(MCccp_AboveCCcutValue_value, np.sqrt(MCccp_AboveCCcutValue_value))

    #### Calcualte Proton number 
    #(1). CCcutTF
    MCproton_above = MCproton[np.where(MCproton[:,0]>float(CCcut_TF))[0]]
    MCproton_CCcutTF_uncertainty = ufloat(MCproton_above.shape[0], np.sqrt(MCproton_above.shape[0]))
    #print("MC: MCproton_CCcutTF_uncertainty: " + str(MCproton_CCcutTF_uncertainty))
    #(2). CCcutValue
    MCproton_AboveCCcutValue = MCproton[np.where(MCproton[:,0]>cc_cut_value)[0]]
    MCproton_AboveCCcutValue_uncertainty = ufloat(MCproton_AboveCCcutValue.shape[0], np.sqrt(MCproton_AboveCCcutValue.shape[0]))

    #### CCLevel WITHOUT Reweighting
    #(1). CCcutTF
    cclevelMC_uncertainty = MCccp_CCcutTF_uncertainty / (MCccp_CCcutTF_uncertainty + MCproton_CCcutTF_uncertainty)
    #print("MC: cclevelMC_uncertainty: " + str(cclevelMC_uncertainty))
    #print("MC: RelativeError:" + str(cclevelMC_uncertainty.s/cclevelMC_uncertainty.n))
    #(2). CCcutValue
    cclevelMC_AboveCCcutValue_uncertainty = MCccp_AboveCCcutValue_uncertainty / (MCccp_AboveCCcutValue_uncertainty + MCproton_AboveCCcutValue_uncertainty) 

    #### CCLevel WITH Reweighting
    #(1). CCcutTF
    MCccp_abovecut_weighted = sum(MCccp_abovecut[:,2])
    if MCccp_abovecut_weighted == 0:
        MCccp_abovecut_weighted = 0.0001
    MCccp_abovecut_Reweighting_uncertainty = unumpy.uarray(MCccp_abovecut_weighted, MCccp_abovecut_weighted * np.sqrt(MCccp_abovecut_value) / MCccp_abovecut_value ) # np.sqrt(sum(MCccp_abovecut_weighted)) too large.
    MCproton_above_Reweighting_uncertainty = unumpy.uarray(sum(MCproton_above[:,2]), sum(MCproton_above[:,2]) * np.sqrt(MCproton_above.shape[0]) / MCproton_above.shape[0] ) 
    cclevelMC_Reweighting_uncertainty = MCccp_abovecut_Reweighting_uncertainty / (MCccp_abovecut_Reweighting_uncertainty + MCproton_above_Reweighting_uncertainty) 
    #(2). CCcutValue
    MCccp_AboveCCcutValue_weighted = sum(MCccp_AboveCCcutValue[:,2])
    if MCccp_AboveCCcutValue_weighted == 0:
        MCccp_AboveCCcutValue_weighted = 0.0001
    MCccp_AboveCCcutValue_Reweighting_uncertainty = unumpy.uarray(MCccp_AboveCCcutValue_weighted, MCccp_AboveCCcutValue_weighted * np.sqrt(MCccp_AboveCCcutValue_value) / MCccp_AboveCCcutValue_value )
    MCproton_CCcutTF_Reweighting_uncertainty = unumpy.uarray(sum(MCproton_AboveCCcutValue[:,2]), sum(MCproton_AboveCCcutValue[:,2]) * np.sqrt(MCproton_AboveCCcutValue.shape[0]) / MCproton_AboveCCcutValue.shape[0] )
    cclevelMC_AboveCCcutValue_Reweighting_uncertainty = MCccp_AboveCCcutValue_Reweighting_uncertainty / (MCccp_AboveCCcutValue_Reweighting_uncertainty + MCproton_CCcutTF_Reweighting_uncertainty)

    return cclevelMC_uncertainty, cclevelMC_Reweighting_uncertainty, MCccp_CCcutTF_uncertainty, MCproton_CCcutTF_uncertainty, MCccp_abovecut_Reweighting_uncertainty, MCproton_above_Reweighting_uncertainty, cclevelMC_AboveCCcutValue_uncertainty, cclevelMC_AboveCCcutValue_Reweighting_uncertainty 


def Calculate_CCProtonOverProton(background_number_uncertainty, protonsignal_above_uncertainty, MCccp_abovecut_uncertainty, MCproton_above_uncertainty):
    ## (3.1) ISS
    CCProtonOverProton_uncertainty_ISS = background_number_uncertainty/protonsignal_above_uncertainty
    ## (3.2) MC
    CCProtonOverProton_uncertainty_MC  = MCccp_abovecut_uncertainty/MCproton_above_uncertainty
    ## (3.3) ISS/MC
    ccratioRatio_uncertainty = CCProtonOverProton_uncertainty_ISS/CCProtonOverProton_uncertainty_MC
    ccratioRatio = ccratioRatio_uncertainty.n
    #ccratioMC = CCProtonOverProton_uncertainty_MC.n
    #print("ISS/MC ccratioRatio: " + str(ccratioRatio) + ", ccratio="  + str(ccratio) + ", ccratioMC=" + str(ccratioMC))
    return CCProtonOverProton_uncertainty_ISS, CCProtonOverProton_uncertainty_MC, ccratioRatio_uncertainty 


#### Plot
def PlotAntiprotonToProtonRatio(resultpath, ISSversion, Binningversion, NNpoint, NN, RigidityBinCenter_450version, CCcut_TF, suffix, pattern, NNsuffix):
    plt.figure(figsize=(18,9))
    plt.xlim(0,550)
    #plt.ylim(0, max(NN)*1.5)
    plt.plot(NNpoint, NN , 'gs',lw=10, markersize=10, label='NeuralNetwork')
    #plt.plot(RigidityBinCenter_450version, binning.published2016 , 'bo',lw=3, label='2016antiprotonpaper')
    plt.errorbar(RigidityBinCenter_450version, binning.published2016, yerr=binning.publishederrorbar ,fmt='o',ecolor='r',color='b',elinewidth=2,capsize=4,label='Published (PRL 2016)')
    #plt.vlines(16.6, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    #plt.vlines(38.9, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    #plt.vlines(147, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    #plt.vlines(175, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    plt.xlabel('Rigidity(GV)',fontsize=30)
    plt.ylabel('Antiproton to proton flux ratio',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.yscale('log')
    plt.grid(True)
    plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'AntiprotonToProton_Ratio_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_'+ str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')

def Plot_Signal_To_Background_Ratio(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, signal_efficiency_value, pattern, NNsuffix):
    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1   
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    for i in range(EffNumber):
        plt.figure(figsize=(18,9))
        plt.yscale('log')
        plt.xlim(0,550)
        if EffNumber == 1:
            plt.errorbar(NNpoint[-8:], signal_to_background_all_diffeff[-8:], yerr=0, xerr=0, fmt='o', label='ISS CCcut=' + str(CCcut_TF) + ')' ) 
        else:
            plt.errorbar(NNpoint[-8:], signal_to_background_all_diffeff[i,-8:], yerr=0, xerr=0, fmt='o', label='ISS CCcut=' + str(CCcut_TF) + ')' ) # label='signal_efficiency='+str(signal_efficiency_value[i])
        plt.xlabel('Rigidity (GV)',fontsize=30)
        plt.ylabel('Signal to Background Ratio',fontsize=30)
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.grid(True, which='both', axis='y')
        plt.legend( loc='best',fontsize=20)
        plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'SignalToBackground_Ratio_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_'+ str(CCcut_TF) + suffix + '_' + str(Binningversion) + "_SignalEff_" + str(signal_efficiency_value[i])  + '.pdf')

def Plot_CCcut_Efficiency(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, cc_cut_value_all_diffeff, signal_efficiency_value, pattern, NNsuffix):
    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    plt.figure(figsize=(18,9))
    plt.xlim(0,550)
    if EffNumber == 1:
        plt.errorbar(NNpoint[-8:], cc_cut_value_all_diffeff[-8:], yerr=None, xerr=0, fmt='o', label='signal_efficiency='+str(signal_efficiency_value[0]))
    else:
        for i in range(signal_to_background_all_diffeff.shape[0]):
            plt.errorbar(NNpoint[-8:], cc_cut_value_all_diffeff[i,-8:], yerr=None, xerr=0, fmt='o', label='signal_efficiency='+str(signal_efficiency_value[i]))
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Charge Confusion Cut',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True)
    plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCcut_Efficiency_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '_' + ISSversion + '.pdf')


def Plot_CCLevel(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, cclevel_all_diffeff, cclevelerror_all_diffeff, signal_efficiency_value, binningCenter_Low, CCLevelMC_AllRigidity_AllEff, CCLevelMCError_AllRigidity_AllEff, CCLevelMC_Reweight_AllRigidity_AllEff, CCLevelMCError_Reweight_AllRigidity_AllEff, ShowPoints, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd):

    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    if binmergenumber == 2:
        NNpoint_Merged = np.append((NNpoint[0:MergeEndIndex:2] + NNpoint[1:MergeEndIndex:2])/2, NNpoint[MergeEndIndex:])
    elif binmergenumber == 5:
        ## Option1: No merge for Fluctuated range
        #NNpoint_Merged = np.append(NNpoint[2:MergeEndIndex:5]                                 , NNpoint[MergeEndIndex:])
        ## Option2: With merge for Fluctuated range
        NNpoint_Merged = NNpoint[2:MergeEndIndex:5]
        NNpoint_Merged = np.append( NNpoint_Merged, (NNpoint[FluctuateBinStart] + NNpoint[FluctuateBinEnd]) / 2 )
        NNpoint_Merged = np.append( NNpoint_Merged, NNpoint[FluctuateBinEnd+1:] )

    NNpoint_Merged                             = NNpoint_Merged[2:]
    cclevel_all_diffeff                        = cclevel_all_diffeff[2:]
    cclevelerror_all_diffeff                   = cclevelerror_all_diffeff[2:]
    CCLevelMC_Reweight_AllRigidity_AllEff      = CCLevelMC_Reweight_AllRigidity_AllEff[2:]
    CCLevelMCError_Reweight_AllRigidity_AllEff = CCLevelMCError_Reweight_AllRigidity_AllEff[2:]

    ## Plot
    fig, ax = plt.subplots(figsize=(18,9))
    plt.yscale('log')
    #plt.xlim(0,550)

    if EffNumber == 1:
        plt.errorbar(NNpoint_Merged, cclevel_all_diffeff                  , yerr=cclevelerror_all_diffeff                  , xerr=0, fmt='o', markersize=9, label='ISS') #label='ISS (CCcut=' + str(CCcut_TF) + ')'
        #plt.errorbar(NNpoint_Merged, CCLevelMC_AllRigidity_AllEff, yerr=CCLevelMCError_AllRigidity_AllEff, xerr=0, fmt='o', markersize=9, label='MC (No Reweighting) (CCcut=' + str(CCcut_TF) + ')')
        plt.errorbar(NNpoint_Merged, CCLevelMC_Reweight_AllRigidity_AllEff, yerr=CCLevelMCError_Reweight_AllRigidity_AllEff, xerr=0, fmt='o', markersize=9, label='MC')  #label='MC (With Reweighting) (CCcut=' + str(CCcut_TF) + ')'
    else:
        for i in range(signal_to_background_all_diffeff.shape[0]):
            plt.errorbar(NNpoint_Merged, cclevel_all_diffeff[i,:], yerr=cclevelerror_all_diffeff[i,:], xerr=0, fmt='o', markersize=9, label='ISS')    #label='ISS (Sig_eff='+str(signal_efficiency_value[i])+')'
            #plt.errorbar(NNpoint_Merged, CCLevelMC_AllRigidity_AllEff[i,:], yerr=CCLevelMCError_AllRigidity_AllEff[i,:], xerr=0, fmt='o', markersize=9, label='MC (No Reweighting) (CCcut=' + str(CCcut_TF) + ')' ) #label='MC (Sig_eff='+str(signal_efficiency_value[i])+")" 
            plt.errorbar(NNpoint_Merged, CCLevelMC_Reweight_AllRigidity_AllEff[i,:], yerr=CCLevelMCError_Reweight_AllRigidity_AllEff[i,:], xerr=0, fmt='o', markersize=9, label='MC') #label='MC (Sig_eff='+str(signal_efficiency_value[i])+")"   

    #plt.errorbar(binningCenter_Low, cclevelMCLow, yerr=cclevelErrorMCLow, xerr=0, fmt='o', markersize=9, label='MC (Low Range) (FullSpan)')
    #plt.errorbar(binningCenter_Low, cclevelMCLow_allPattern, yerr=cclevelErrorMCLow_allPattern, xerr=0, fmt='o', markersize=9, label='MC (Low Range) (All Pattern)')
    #plt.errorbar(NNpoint, cclevelMCHigh, yerr=cclevelErrorMCHigh, xerr=0, fmt='o', markersize=9, label='MC (High Range)')

    #plt.tick_params(axis='x', which='both', labelsize=20)
    #ax.xaxis.set_major_formatter(ScalarFormatter())
    #ax.xaxis.set_minor_formatter(ScalarFormatter())
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Charge Confusion Level',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    #plt.grid(True, which='both', axis='y')
    plt.legend( loc='best',fontsize=20)

    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCLevel(XLinear)_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '_' + ISSversion +'.pdf')
    plt.xscale('log')
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCLevel(XLog)_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '_' + ISSversion + '.pdf')


def Plot_CCLevel_ISSOverMC(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, CCLevelRatio_all_diffeff, CCLevelRatioError_AllRigidity_AllEff, CCLevelRatio_Reweight_AllRigidity_AllEff, CCLevelRatioError_Reweight_AllRigidity_AllEff, CCLevelRatio_AboveCCcutValue_AllRigidity_AllEff, CCLevelRatioError_AboveCCcutValue_AllRigidity_AllEff, CCLevelRatio_AboveCCcutValue_Reweight_AllRigidity_AllEff, CCLevelRatioError_AboveCCcutValue_Reweight_AllRigidity_AllEff, signal_efficiency_value, ShowPoints, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd):
    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    if binmergenumber == 2:
        NNpoint_Merged = np.append((NNpoint[0:MergeEndIndex:2] + NNpoint[1:MergeEndIndex:2])/2, NNpoint[MergeEndIndex:])
    elif binmergenumber == 5:
        ## Option1: No merge for Fluctuated range
        #NNpoint_Merged = np.append(NNpoint[2:MergeEndIndex:5]                                 , NNpoint[MergeEndIndex:])
        ## Option2: With merge for Fluctuated range
        NNpoint_Merged = NNpoint[2:MergeEndIndex:5]
        NNpoint_Merged = np.append( NNpoint_Merged, (NNpoint[FluctuateBinStart] + NNpoint[FluctuateBinEnd]) / 2 )
        NNpoint_Merged = np.append( NNpoint_Merged, NNpoint[FluctuateBinEnd+1:] )

    ## Plot
    plt.figure(figsize=(18,9))
    #plt.xlim(0,550)
    ax=plt.gca()

    if EffNumber == 1:
        #plt.errorbar(NNpoint_Merged, CCLevelRatio_all_diffeff                                , yerr=CCLevelRatioError_AllRigidity_AllEff         , xerr=0, fmt='o', markersize=9, label='ISS/MC (No Reweighting) (CCcut=' + str(CCcut_TF) + ')' ) # label='ISS/MC (Sig_eff='+str(signal_efficiency_value[0])+")"
        plt.errorbar(NNpoint_Merged, CCLevelRatio_Reweight_AllRigidity_AllEff                , yerr=CCLevelRatioError_Reweight_AllRigidity_AllEff, xerr=0, fmt='o', markersize=9, label='ISS/MC') # label='ISS/MC (Sig_eff='+str(signal_efficiency_value[0])+")"
        #plt.errorbar(NNpoint_Merged, CCLevelRatio_AboveCCcutValue_AllRigidity_AllEff         , yerr=CCLevelRatioError_AboveCCcutValue_AllRigidity_AllEff, xerr=0, fmt='o', markersize=9, label='ISS/MC (No Reweighting)  (Sig_eff=' + str(signal_efficiency_value[0]) + ')' )
        #plt.errorbar(NNpoint_Merged, CCLevelRatio_AboveCCcutValue_Reweight_AllRigidity_AllEff, yerr=CCLevelRatioError_AboveCCcutValue_Reweight_AllRigidity_AllEff, xerr=0, fmt='o', markersize=9, label='ISS/MC (With Reweighting)  (Sig_eff=' + str(signal_efficiency_value[0]) + ')' )

    else:
        for i in range(signal_to_background_all_diffeff.shape[0]):
            #plt.errorbar(NNpoint_Merged, CCLevelRatio_all_diffeff[i,:]                , yerr=CCLevelRatioError_AllRigidity_AllEff[i,:]         , xerr=0, fmt='o', markersize=9, label='ISS/MC (No Reweighting) (CCcut=' + str(CCcut_TF) + ')') #label='ISS/MC (Sig_eff='+str(signal_efficiency_value[i])+")"
            plt.errorbar(NNpoint_Merged, CCLevelRatio_Reweight_AllRigidity_AllEff[i,:], yerr=CCLevelRatioError_Reweight_AllRigidity_AllEff[i,:], xerr=0, fmt='o', markersize=9, label='ISS/MC (With Reweighting) (CCcut=' + str(CCcut_TF) + ')') #label='ISS/MC (Sig_eff='+str(signal_efficiency_value[i])+")"
            #plt.errorbar(NNpoint_Merged, CCLevelRatio_AboveCCcutValue_AllRigidity_AllEff[i,:], yerr=CCLevelRatioError_AboveCCcutValue_AllRigidity_AllEff[i,:], xerr=0, fmt='o', markersize=9, label='ISS/MC (No Reweighting)  (Sig_eff=' + str(signal_efficiency_value[i]) + ')' )
            #plt.errorbar(NNpoint_Merged, CCLevelRatio_AboveCCcutValue_Reweight_AllRigidity_AllEff[i,:], yerr=CCLevelRatioError_AboveCCcutValue_Reweight_AllRigidity_AllEff[i,:], xerr=0, fmt='o', markersize=9, label='ISS/MC (With Reweighting)  (Sig_eff=' + str(signal_efficiency_value[i]) + ')' )

    ylim = ax.get_ylim()
    minvalue = ylim[0]
    maxvalue = ylim[1]
    plt.vlines(175,  minvalue, maxvalue, linestyles = "dotted")

    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Charge Confusion Level (ISS/MC)',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True, which='both', axis='y')
    plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCLevel_ISSOverMC(XLinear)_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '_' + ISSversion + '.pdf')
    plt.xscale('log')
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCLevel_ISSOverMC(XLog)_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '_' + ISSversion + '.pdf')


def Plot_CCprotonNumbers(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, ccprotonnumber_all_diffeff, signal_efficiency_value, ShowPoints, pattern, NNsuffix):
    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    plt.figure(figsize=(18,9))
    ax=plt.gca()
    if EffNumber == 1:
        #### FIXME:First two points has limited statistics, temporary solution to plot.
        ccprotonnumber_all_diffeff[0] = ccprotonnumber_all_diffeff[2] * 0.8
        ccprotonnumber_all_diffeff[1] = ccprotonnumber_all_diffeff[2] * 0.6
        plt.errorbar(NNpoint[-ShowPoints:], ccprotonnumber_all_diffeff[-ShowPoints:] ,   yerr=np.sqrt(ccprotonnumber_all_diffeff[-ShowPoints:]), xerr=0, fmt='o', markersize=9, label='CCcut=' + str(CCcut_TF) ) # label='ISS (Sig_eff='+str(signal_efficiency_value[0])+")"
    else:
        for i in range(signal_to_background_all_diffeff.shape[0]):
            #### FIXME:First two points has limited statistics, temporary solution to plot.
            ccprotonnumber_all_diffeff[i,0] = ccprotonnumber_all_diffeff[i,2] * 0.8
            ccprotonnumber_all_diffeff[i,1] = ccprotonnumber_all_diffeff[i,2] * 0.6
            plt.errorbar(NNpoint[-ShowPoints:], ccprotonnumber_all_diffeff[i,-ShowPoints:] ,   yerr=np.sqrt(ccprotonnumber_all_diffeff[i,-ShowPoints:]), xerr=0, fmt='o', markersize=9, label='CCcut=' + str(CCcut_TF) )

    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Charge Cofunsed Proton Numbers',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)

    plt.xlim(10,550)
    ax.set_xscale('log')
    ax.set_xticks([20, 30, 50, 70, 100, 200, 400])
    ax.get_xaxis().set_major_formatter(ScalarFormatter())

    plt.yscale('log')
    #plt.grid(True, which='both', axis='y')
    plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCProtonNumbers(XLog)_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')

def Plot_ProtonNumbers(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, protonnumber_all_diffeff, ShowPoints, signal_efficiency_value, pattern, NNsuffix):
    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    plt.figure(figsize=(18,9))
    ax=plt.gca()
    if EffNumber == 1:
        plt.errorbar(NNpoint[-ShowPoints:], protonnumber_all_diffeff[-ShowPoints:] ,   yerr=np.sqrt(protonnumber_all_diffeff[-ShowPoints:]), xerr=0, fmt='o', markersize=9, label='CCcut=' + str(CCcut_TF)) # label='ISS (Sig_eff='+str(signal_efficiency_value[0])+")"
    else:
        for i in range(signal_to_background_all_diffeff.shape[0]):
            plt.errorbar(NNpoint[-ShowPoints:], protonnumber_all_diffeff[i,-ShowPoints:] ,   yerr=np.sqrt(protonnumber_all_diffeff[i,-ShowPoints:]), xerr=0, fmt='o', markersize=9, label='CCcut=' + str(CCcut_TF))

    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Proton Numbers',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)

    plt.xlim(10,550)
    ax.set_xscale('log')
    ax.set_xticks([20, 30, 50, 70, 100, 200, 400])
    ax.get_xaxis().set_major_formatter(ScalarFormatter())

    plt.yscale('log')
    #plt.grid(True, which='both', axis='y')
    plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'ProtonNumbers(XLog)_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')


def Plot_CCProton_Proton_Ratio(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, CCPToProtonRatio_ISS, CCPToProtonRatioError_ISS, CCPToProtonRatio_MC, CCPToProtonRatioError_MC, signal_efficiency_value, ShowPoints, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd):
    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    if binmergenumber == 2:
        NNpoint_Merged = np.append((NNpoint[0:MergeEndIndex:2] + NNpoint[1:MergeEndIndex:2])/2, NNpoint[MergeEndIndex:])
    elif binmergenumber == 5:
        ## Option1: No merge for Fluctuated range
        #NNpoint_Merged = np.append(NNpoint[2:MergeEndIndex:5]                                 , NNpoint[MergeEndIndex:])
        ## Option2: With merge for Fluctuated range
        NNpoint_Merged = NNpoint[2:MergeEndIndex:5]
        NNpoint_Merged = np.append( NNpoint_Merged, (NNpoint[FluctuateBinStart] + NNpoint[FluctuateBinEnd]) / 2 )
        NNpoint_Merged = np.append( NNpoint_Merged, NNpoint[FluctuateBinEnd+1:] )

    ## Plot
    plt.figure(figsize=(18,9))
    plt.yscale('log')
    plt.xlim(0,550)

    if EffNumber == 1:
        plt.errorbar(NNpoint_Merged, CCPToProtonRatio_ISS, yerr=CCPToProtonRatioError_ISS, xerr=0, fmt='o', markersize=9, label='ISS (Sig_eff='+str(signal_efficiency_value[0])+")") 
        plt.errorbar(NNpoint_Merged, CCPToProtonRatio_MC , yerr=CCPToProtonRatioError_MC , xerr=0, fmt='o', markersize=9, label='MC (Sig_eff='+str(signal_efficiency_value[0])+")")
    else:
        for i in range(signal_to_background_all_diffeff.shape[0]):
            plt.errorbar(NNpoint_Merged, CCPToProtonRatio_ISS[i,:], yerr=CCPToProtonRatioError_ISS[i,:], xerr=0, fmt='o', markersize=9, label='ISS (Sig_eff='+str(signal_efficiency_value[i])+")") 
            plt.errorbar(NNpoint_Merged, CCPToProtonRatio_MC[i,:] , yerr=CCPToProtonRatioError_MC[i,:] , xerr=0, fmt='o', markersize=9, label='MC (Sig_eff='+str(signal_efficiency_value[i])+")")
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('CCproton/Proton Ratio',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True, which='both', axis='y')
    plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCProtonToProtonRatio_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')


def Plot_CCProton_Proton_Ratio_ISSOverMC(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, signal_to_background_all_diffeff, NNpoint, ccratioRatio_all_diffeff, ccratioRatioError_all_diffeff, ShowPoints, signal_efficiency_value, pattern, NNsuffix, binmergenumber, MergeEndIndex, FluctuateBinStart, FluctuateBinEnd):
    if signal_to_background_all_diffeff.ndim == 1:
        EffNumber = 1
    else:
        EffNumber = signal_to_background_all_diffeff.shape[0]

    if binmergenumber == 2:
        NNpoint_Merged = np.append((NNpoint[0:MergeEndIndex:2] + NNpoint[1:MergeEndIndex:2])/2, NNpoint[MergeEndIndex:])
    elif binmergenumber == 5:
        ## Option1: No merge for Fluctuated range
        #NNpoint_Merged = np.append(NNpoint[2:MergeEndIndex:5]                                 , NNpoint[MergeEndIndex:])
        ## Option2: With merge for Fluctuated range
        NNpoint_Merged = NNpoint[2:MergeEndIndex:5]
        NNpoint_Merged = np.append( NNpoint_Merged, (NNpoint[FluctuateBinStart] + NNpoint[FluctuateBinEnd]) / 2 )
        NNpoint_Merged = np.append( NNpoint_Merged, NNpoint[FluctuateBinEnd+1:] )

    ## Plot
    plt.figure(figsize=(18,9))
    plt.xlim(0,550)
    plt.ylim(0,2)

    if EffNumber == 1:
        plt.errorbar(NNpoint_Merged, ccratioRatio_all_diffeff, yerr=ccratioRatioError_all_diffeff, xerr=0, fmt='o', markersize=9, label='ISS/MC (Sig_eff='+str(signal_efficiency_value[0])+")")
    else:
        for i in range(signal_to_background_all_diffeff.shape[0]):
            plt.errorbar(NNpoint_Merged, ccratioRatio_all_diffeff[i,:], yerr=ccratioRatioError_all_diffeff[i,:], xerr=0, fmt='o', markersize=9, label='ISS/MC (Sig_eff='+str(signal_efficiency_value[i])+")")
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('CCproton/Proton (ISS/MC)',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True, which='both', axis='y')
    plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'CCProtonToProtonRatio_ISSOverMC_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_' + str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')


def Plot_Sta_Error(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, NNpoint, RigidityBinCenter_450version, delta_antiproton, proton_number_TF, NN, pattern, NNsuffix):
    sta_error_2016_truevalue = (np.array([0.04, 0.04, 0.05, 0.05, 0.06, 0.06, 0.06, 0.07, 0.08, 0.08, 0.07, 0.08, 0.10, 0.12, 0.17, 0.26, 0.34, 0.34]) * 10**(-4)).astype(np.float)
    sta_error_2016_relative  = sta_error_2016_truevalue/binning.published2016[-sta_error_2016_truevalue.shape[0]:]
    NN_sta_error_truevalue   = delta_antiproton/proton_number_TF[-delta_antiproton.shape[0]:]
    NN_sta_error_relative    = NN_sta_error_truevalue/NN[-NN_sta_error_truevalue.shape[0]:]

    plt.figure(figsize=(18,9))
    plt.xlim(0,550)
    plt.plot(NNpoint[-NN_sta_error_relative.shape[0]:], NN_sta_error_relative*100, 'gs', markersize=10,label='NeuralNetwork')
    plt.plot(RigidityBinCenter_450version[-sta_error_2016_relative.shape[0]:], sta_error_2016_relative*100, 'bo',markersize=10, label='Published (PRL 2016)')
    #plt.vlines(16.6, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    #plt.vlines(38.9, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    #plt.vlines(147, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    #plt.vlines(175, 0, max(NN)*1.5, colors = "c", linestyles = "dashed")
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Relative Error [%]',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True)
    #plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'StatisticalError_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_'+ str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')

def Plot_Antiproton_Number(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, NNpoint, antiproton_number, pattern, NNsuffix):
    antiproton_number_fullspan = antiproton_number[-4:]

    plt.figure(figsize=(18,9))
    plt.xlim(0,550)
    plt.plot(NNpoint[-antiproton_number_fullspan.shape[0]:], antiproton_number_fullspan, 'gs', markersize=10,label='NeuralNetwork')
    #plt.vlines(16.6, 0, max(binning.antiproton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.vlines(38.9, 0, max(binning.antiproton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.vlines(147, 0, max(binning.antiproton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.vlines(175, 0, max(binning.antiproton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Antiproton Numbers',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True)
    plt.yscale('log')
    #plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'AntiprotonNumbers_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_'+ str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')

def Plot_Proton_number(resultpath, ISSversion, Binningversion, CCcut_TF, suffix, NNpoint, proton_number_fullspan, pattern, NNsuffix):
    proton_number_fullspan     = proton_number_TF[-4:]

    plt.figure(figsize=(18,9))
    plt.xlim(0,550)
    plt.plot(NNpoint[-proton_number_fullspan.shape[0]:], proton_number_fullspan, 'gs', markersize=10,label='NeuralNetwork')
    #plt.vlines(16.6, 0, max(binning.proton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.vlines(38.9, 0, max(binning.proton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.vlines(147, 0, max(binning.proton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.vlines(175, 0, max(binning.proton_number_published)*2, colors = "c", linestyles = "dashed")
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('Proton Numbers',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True)
    plt.yscale('log')
    #plt.legend( loc='best',fontsize=20)
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + 'ProtonNumbers_Pattern_' + pattern + NNsuffix + '_' + 'CCcut_TF_'+ str(CCcut_TF) + suffix + '_' + str(Binningversion) + '.pdf')

def ListToVector(rawlist):
    vector = ROOT.vector('double')(len(rawlist))
    vector.clear()
    for i in range(len(rawlist)):
        vector.insert(vector.begin()+i, rawlist[i])
    return vector


def VStack(array_all, array):
    if array_all.shape[0] == 0:
        array_all = array
    else:
        array_all = np.vstack((array_all, array ))
    return array_all









