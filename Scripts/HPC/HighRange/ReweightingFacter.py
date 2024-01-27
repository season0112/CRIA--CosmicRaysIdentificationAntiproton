#!/usr/bin/env python
import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TChain, gDirectory, std
from root_numpy import root2array, tree2array, fill_hist
import binning
import matplotlib.pyplot as plt
from ctypes import *
import uproot
import ReweightingFacter_function
import os
import argparse

def main():

    #### Get MC  R_true histogram/BinningWidth
    print("\033[34mLoadMCHist Now. \033[0m")
    h_MCPrimaryMomentum, h_MCPrimaryMomentum_more, MCWeight, MCPrimaryMomentum = ReweightingFacter_function.LoadMCHist(highpath, binningWidth, LowEdge, HighEdge, BinNumber, pattern, MCName)

    #### Reference
    print()
    print("\033[34mLoadProtonRPL2015 Now. \033[0m")
    h_ProtonPRL, h_ProtonPRL_PowerLawFit, h_ProtonPRL_AllRange = ReweightingFacter_function.LoadProtonRPL2015(resultpath, LowEdge, HighEdge, BinNumber)

    #### Ratio Calculation
    ReweightingFactor = TH1D("", "", 31, binning.Newbinnings_525_zhili[26:-1])
    #ReweightingFactor.Divide(h_MCPrimaryMomentum, h_ProtonPRL)
    ReweightingFactor.Divide(h_ProtonPRL, h_MCPrimaryMomentum)
    for i in range(ReweightingFactor.GetNbinsX()):
        ReweightingFactor.SetBinContent(i+1, ReweightingFactor.GetBinContent(i+1))
        ReweightingFactor.SetBinError(i+1, 0)

    ReweightingFactor_More = TH1D("", "", BinNumber, LowEdge, HighEdge)
    #ReweightingFactor_More.Divide(h_MCPrimaryMomentum_more, h_ProtonPRL_PowerLawFit)
    ReweightingFactor_More.Divide(h_ProtonPRL_PowerLawFit, h_MCPrimaryMomentum_more)
    for i in range(ReweightingFactor_More.GetNbinsX()):
        ReweightingFactor_More.SetBinContent(i+1, ReweightingFactor_More.GetBinContent(i+1))
        ReweightingFactor_More.SetBinError(i+1, 0)

    #### Save Reweighting in ROOT file
    ReweightingFile = TFile.Open(resultpath + '/antiproton_to_proton_ratio_plot/' + "ReweightingFile_" + MCName + "_Pattern_" + str(pattern) + ".root", "RECREATE")
    ReweightingFactor.Write("ReweightingFactor")
    ReweightingFactor_More.Write("ReweightingFactor_More")

    h_MCPrimaryMomentum.Write("h_MCPrimaryMomentum")
    h_MCPrimaryMomentum_more.Write("h_MCPrimaryMomentum_more")

    h_ProtonPRL.Write("h_ProtonPRL")
    h_ProtonPRL_PowerLawFit.Write("h_ProtonPRL_PowerLawFit")
    h_ProtonPRL_AllRange.Write("h_ProtonPRL_AllRange")
    ReweightingFile.Close()

    #### Plot
    ReweightingFacter_function.Plot_MCPrimaryMomentum(resultpath, MCName, pattern, h_MCPrimaryMomentum   , h_MCPrimaryMomentum_more)
    ReweightingFacter_function.Plot_ProtonReference(resultpath  , MCName, pattern, h_ProtonPRL           , h_ProtonPRL_PowerLawFit, h_ProtonPRL_AllRange)
    ReweightingFacter_function.Plot_ReweightFactor(resultpath   , MCName, pattern, ReweightingFactor     , "AntiprotonBins")
    ReweightingFacter_function.Plot_ReweightFactor(resultpath   , MCName, pattern, ReweightingFactor_More, "MoreBins")
    ReweightingFacter_function.Plot_MCWeight(resultpath         , MCName, pattern, MCWeight              , MCPrimaryMomentum, LowEdge, HighEdge, BinNumber)

if __name__ == '__main__':
    #### Parser Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--pattern', help='which tracker pattern you choose')
    arguments = parser.parse_args()

    if (arguments.pattern):
        pattern = arguments.pattern
    else:
        print("You need to choose a tracker pattern!")
        os._exit(0)

    binnenCenter = binning.Newbinnings_525_center_zhili[26:-1]
    binningWidth = binning.NewbinningWidth_525_zhili[26:-1]

    highpath = os.getenv('HPCHIGHENERGYDATADIR')
    resultpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis'

    #MCName = 'B1042_pr.pl1.1800_7.6_all'
    MCName = 'B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all'
    #MCName = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all'

    if MCName == 'B1042_pr.pl1.1800_7.6_all':
        LowEdge   = 10
        HighEdge  = 1000
        BinNumber = 500
    elif MCName == 'B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all':
        LowEdge   = 50
        HighEdge  = 16000
        BinNumber = 5000
    elif MCName == 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all':
        LowEdge   = 10
        HighEdge  = 16000
        BinNumber = 5000

    main()








