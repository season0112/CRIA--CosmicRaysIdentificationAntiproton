#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TGraphErrors
from root_numpy import fill_hist
from root_numpy import root2array, tree2array
import PythonPlotDefaultParameters
import binning
import os
import CCEstimator_tool


def main():

    #### Load Data
    ## Fit result Parameters
    # For Aachen group meeting around Apr.8.2021:best setting: 211-259GV:CCCut = "0.70",CCN = "9",TRDN = "11"; 259-330 and 330-525:others, need to try;  binnumber = 20
    '''
    CCCut_all = ['0.00', '0.20', '0.25', '0.30', '0.35', '0.40', '0.45']
    CCN_all   = ['20', '20', '20', '9']
    TRDN_all  = ['20', '16', '12', '11']
    '''
    CCCut_all = ['0.20']
    CCN_all   = ['20']
    TRDN_all  = ['20']    

    ## Free Parameters
    #binnumber_all    = [20, 25, 30, 35, 40, 50] 
    binnumber_all    = [30]

    TRDLikelihoodCut = 0.7
    MarkerSize       = 30
    Rigiditybin = binning.bins # 32, 14.1-525

    #CorrectName  = "MCChargeCorrect_B1042_antipr.pl1.1800_7.6_all_Tree_"
    CorrectName  = "MCChargeCorrect_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_"
    #CorrectName  = "MCChargeCorrect_B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_positive_"

    ConfusedName = "MCChargeConfused_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_"
    #ConfusedName = "MCChargeConfused_B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_negative_"

    #### main loop()
    for binnumber in binnumber_all:
        for CCCut in CCCut_all:
            for CCN, TRDN in zip(CCN_all, TRDN_all):

                    a_number = CCEstimator_tool.LoadFitResult(highpath, CCCut, CCN, TRDN)

                    #### loop last three bins ################################
                    for i in range(-3,0): # 211-259,259-330,330-525
                    #for i in range(-15,0):

                        Correct_MC    = np.load(highpath + "/ISS_anylsis/data/" + CorrectName  + str(Rigiditybin[i]) + ".npy")
                        Confused_MC   = np.load(highpath + "/ISS_anylsis/data/" + ConfusedName + str(Rigiditybin[i]) + ".npy")
                        Correct_data  = np.load(highpath + "/ISS_anylsis/data/plot_ISS_positive_rigidity" + str(Rigiditybin[i]) + "_pass7.8.npy")
                        Confused_data = np.load(highpath + "/ISS_anylsis/data/plot_ISS_negative_rigidity" + str(Rigiditybin[i]) + "_pass7.8.npy")

                        ccpercentage = a_number[i]/Confused_data.shape[0] 
                        print("ccpercentage:" + str(ccpercentage))

                        Correct_MC    = Correct_MC   [np.where(Correct_MC[:,1]   >TRDLikelihoodCut)[0]]
                        Confused_MC   = Confused_MC  [np.where(Confused_MC[:,1]  >TRDLikelihoodCut)[0]]
                        Correct_data  = Correct_data [np.where(Correct_data[:,1] >TRDLikelihoodCut)[0]]
                        Confused_data = Confused_data[np.where(Confused_data[:,1]>TRDLikelihoodCut)[0]]

                        factor_CorrectConfusedCountRatio = Correct_MC.shape[0] / Confused_MC.shape[0]

                        #### Plot CC Estimator
                        CCEstimator_tool.Plot_All(Confused_data, Correct_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, MarkerSize, factor_CorrectConfusedCountRatio, CCCut, CCN, TRDN)

                        #### Plot Rejection Power
                        FractionSignalSet, RejectionBackgroundSet = CCEstimator_tool.CalculateRejectionPower(binnumber, Correct_MC[:,0], Confused_MC[:,0].repeat(factor_CorrectConfusedCountRatio))
                        #print("FractionSignalSet:" + str(FractionSignalSet))
                        CCEstimator_tool.Plot_RejectionPower(FractionSignalSet, RejectionBackgroundSet, highpath, Rigiditybin, i, binnumber, CCCut, CCN, TRDN)


if __name__ == '__main__':

    #highpath = os.getenv('HPCHIGHENERGYDATADIR')
    highpath = "/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0"

    main()






