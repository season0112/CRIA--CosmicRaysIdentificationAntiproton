#!/usr/bin/env python
import ROOT
from ROOT import TMVA, TFile, TTree, TCut, TString, TBranch, gROOT, TChain
import matplotlib.pyplot as plt
import numpy as np
import binning
from root_numpy import array2tree, array2root, tree2array, root2array
import os
import TMVA_Application_model


def main():
    #### Parameters
    datapath = os.getenv('HPCHIGHENERGYDATADIR')

    #Correct_MC_name = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive'                    # Large
    Correct_MC_name = 'B1042_antipr.pl1.1800_7.6_all_Tree'                                       # Default
    #Correct_MC_name = 'B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_positive'          # Test
    #Confused_MC_name = 'B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_negative'          # Test
    Confused_MC_name = 'B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative'                     # Default
    Electron_MC_name = 'B1091_el.pl1.0_25_2000_7.6_all_Tree'                                      # Default

    chain_Proton   = TChain('ExampleAnalysisTree','')
    chain_CCproton = TChain('ExampleAnalysisTree','')
    chain_Electron = TChain('ExampleAnalysisTree','')

    LowMergedList = binning.bins[0:23]
    for i in LowMergedList:
        chain_Proton.AddFile( datapath + "/" + Correct_MC_name + "_" + str(i) + ".root")
        chain_CCproton.AddFile( datapath + "/" + Confused_MC_name + "_" + str(i) + ".root")
        chain_Electron.AddFile( datapath + "/" + Electron_MC_name + "_" + str(i) + ".root")

    chain_Proton.Merge(datapath + "/" + Correct_MC_name + "_" + "14.1_80.5" + ".root")
    chain_CCproton.Merge(datapath + "/" + Confused_MC_name + "_" + "14.1_80.5" + ".root")
    chain_Electron.Merge(datapath + "/" + Electron_MC_name + "_" + "14.1_80.5" + ".root")


if __name__ == "__main__":
    main()



