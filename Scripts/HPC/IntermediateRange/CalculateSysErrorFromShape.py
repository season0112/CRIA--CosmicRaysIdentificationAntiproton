#!/usr/bin/env python
import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F, gStyle
import binning
import argparse 
import os
import matplotlib.pyplot as plt
import PythonPlotDefaultParameters
import CalculateSysErrorFromShape_tool
import CalculateSysErrorFromShape_Plot


def main():

    for splitindex in range(0, SplitTotal, mergestep):  ## If time dependent, loop in time index, for time averaged, it's 1.

        #### Load Fit Result
        if Timemode == "TimeAveraged":
            if RigidityRange == 'low':
                rootfile = TFile(LowPath          + "/totalall/Time_Averaged_ratio_Low/binmerge1/plots/Ratio_" + ISSversion + "_uncertainty.root", "READ")
            elif RigidityRange == 'intermediate':
                rootfile = TFile(IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/intermediate_0124_free_" + ISSversion + "binmerge1_uncertainty.root", "READ") 
        elif Timemode == "TimeDependent":
            if RigidityRange == 'low':
                rootfile = TFile(LowPath          + "/totalall/SysError_Shape/Ratio_TimeIndex" +  str(splitindex) + ".root", "READ")

        ## Create list to hold SysErrorShape 
        SysErrorShape_original = []

        ## Loop in Rigidity and Randomlised Templates
        for RIndex in range(UsedBinningCenter.shape[0]):
            CalculateSysErrorFromShape_tool.FillinHist(rootfile, BinsNumber, Ratio_min, Ratio_max, generatenumber, RigidityRange, SysErrorShape_original, RIndex, UsedBinningEdge, Timemode)
        rootfile.Close()
        ## Manually fix for 1.33-1.51GV
        if Timemode == "TimeAveraged":
            if RigidityRange == 'low':
                SysErrorShape_original[2] = (SysErrorShape_original[3] + SysErrorShape_original[1])/2
            if RigidityRange == 'intermediate':
                SysErrorShape_original[4] = (SysErrorShape_original[3] + SysErrorShape_original[6])/2
                SysErrorShape_original[5] = (SysErrorShape_original[3] + SysErrorShape_original[6])/2
                SysErrorShape_original[9] = (SysErrorShape_original[8] + SysErrorShape_original[10])/2
        ## Plot RMS
        #Covnumber=2 # Defalte
        #CalculateSysErrorFromShape_Plot.Plot_RMS(Covnumber, UsedBinningCenter, SysErrorShape_original, RigidityRange)
        ## Decide use smoothed curve to save
        if Timemode == "TimeAveraged":
            if RigidityRange == 'low':
                SysErrorShape = SysErrorShape_original
            if RigidityRange == 'intermediate':
                 Covnumber = 20
                 SysErrorShape = CalculateSysErrorFromShape_tool.smooth(SysErrorShape_original,Covnumber)
        elif Timemode == "TimeDependent":
            SysErrorShape = SysErrorShape_original
            #Covnumber = 2
            #SysErrorShape = CalculateSysErrorFromShape_tool.smooth(SysErrorShape_original,Covnumber)

        ## Make TGraph and Save into Root File
        SysErrorShape = np.array(SysErrorShape)
        print("SysErrorShape: " + str(SysErrorShape))
        g_SysErrorShape  = TGraph(UsedBinningCenter.shape[0], UsedBinningCenter, SysErrorShape)

        if Timemode == "TimeAveraged":
            if RigidityRange == 'low':
                f_SysErr = TFile(LowPath + "/totalall/Time_Averaged_ratio_Low/binmerge1/plots/SysErr_Shape_" + ISSversion + ".root", "RECREATE")
            elif RigidityRange == 'intermediate':
                f_SysErr = TFile(IntermediatePath + "/total/Time_Averaged_ratio/binmerge1/SysErr_Shape_" + ISSversion + ".root", "RECREATE")
            g_SysErrorShape.Write("g_SysErrorShape")
            f_SysErr.Close()
        elif Timemode == "TimeDependent":
            if RigidityRange == 'low':
                f_SysErr = TFile(LowPath + "/totalall/SysError_Shape/SysError_Shape_TimeIndex" +  str(splitindex) + ".root", "RECREATE")
            g_SysErrorShape.Write("g_SysErrorShape")
            f_SysErr.Close()

if __name__ == '__main__':

    #### Parser Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('--Range'         , help='Rigidity Range, low, intermediate or high')
    parser.add_argument('--issversion'    , help='ISS data version, pass7.8 or 2016paper')
    parser.add_argument('--GenerateNumber', type=int, help='how many tempalte are generated to calculate sys err due to template shape')
    parser.add_argument('--Timemode'      , help='Time averaged or time dependent.')   
    parser.add_argument('--TimeSplitMode' , help='which TimeSplitMode')
    arguments = parser.parse_args()

    if (arguments.Range and arguments.issversion and arguments.GenerateNumber or arguments.Timemode):
        ISSversion     = arguments.issversion
        RigidityRange  = arguments.Range
        generatenumber = arguments.GenerateNumber
        Timemode       = arguments.Timemode
        TimeSplitMode = arguments.TimeSplitMode
    else:
        print("You need to provide all parameters!")
        os._exit(0)

    IntermediatePath = os.getenv('HPCINTERMEDIATEDIR')
    LowPath          = os.getenv('HPCLOWENERGYDIR')

    if RigidityRange == 'intermediate':
        BinningAll       = binning.binslow[10:]
        BinningAll[19]   = '16.6_18'
        BinningCenterAll = binning.published2016binnings_center[9:29]
    elif RigidityRange == 'low':
        BinningAll       = binning.binslow[1:17]
        BinningAll[0]    = '1_1.16'
        BinningCenterAll = binning.published2016binnings_center[0:16]


    ## Binnumbers, R_min and R_max only used for plot, doesn't affect RMS
    ## One bin to show for 1.51-1.71GV
    #Ratio_min  = 2
    #Ratio_max  = 2.6
    Ratio_min  = 0
    Ratio_max  = 50
    BinsNumber = 50


    if Timemode == "TimeAveraged":
        UsedBinningEdge = BinningAll
        UsedBinningCenter = BinningCenterAll
        SplitTotal = 1
        mergestep = 1

    elif Timemode == "TimeDependent":
        if RigidityRange == 'intermediate':
            print("in progress")
            os._exit(0)
        elif RigidityRange == 'low':
            UsedBinningEdge   = binning.BinsEdgeName_TimeDependent_Low
            UsedBinningCenter = binning.BinsCenter_TimeDependent_Low

        if TimeSplitMode == "3BartalRotation":
            SplitTotal = 123
            mergestep = 3
        elif TimeSplitMode == "6BartalRotation":
            SplitTotal = 121
            mergestep = 6
        elif TimeSplitMode == "6Months":
            SplitTotal = 20
            mergestep = 1 

    main()


