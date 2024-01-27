#!/usr/bin/env python
import numpy as np
import uproot
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import matplotlib.dates as md
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import argparse
import binning
import TimeDependentTimeRange
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import numpy.polynomial.polynomial as poly
from scipy import optimize
import PythonPlotDefaultParameters
import load
import draw

import LoadExcel

import sys
sys.path.append('/home/bo791269/Software/AntiprotonAnalysis/Scripts/HPC/others/TimeDependentPlot/plotmodule/')
import drawEleOverPosRatio
import drawCompareWithElectron
import drawChi2
import drawmyresult
import drawCompareIHEPTaiwanResult
import drawCompareWithOtherGroup
import drawStatisticalErrorComparision
import drawBreakDownOfTotalError
import drawProtontoAntiproton
import drawRelativeResult
import drawTRDTOFCompare
import drawProtontoAntiprotonRelativeResult
import drawAllbins
import drawAllbinsRelative
import drawPosOverEleRatio
import drawFitNumbers
import drawPbarRatioInAll
import drawRelativeErrorvsRigidity
import drawModel
import MakeNumpyForTimeDependent

def main():


    ##################################################################
    # You need to manually appoint a electron positron result index! 
    # This one only affect originll electron positron result, for rebinned one, no affect. 
    ##################################################################

    ##### Load Electron Positron Result
    ElectronPositronIndex = 13
    ## Load Niko Electron Positron Result
    TFile_lepton_result      , PosiOverElec, PosiOverElec_error, Electron, Electron_error_relative, Positron, Positron_error_relative = load.LoadElectronPositronResult(ElectronPositronIndex)
    ## Load Maura Electron Positron Result
    TFile_lepton_result_Maura, Electron_Maura, ElectronError_Maura, Positron_Maura, PositronError_Maura                               = load.LoadMauraElectronPositronResultMonthly(ElectronPositronIndex)
 
    #### Load Other Group Result
    Ihep_all, Ihep_error_all, Taiwan_all, Taiwan_error_all, Taiwan_all_6, Taiwan_error_all_6 = load.LoadOtherGroupResult()

    #### Load Time Averaged ratio 
    time_averaged_ratio, time_averaged_ratio_with_effective, time_averaged_error, time_averaged_proton, time_averaged_antiproton = load.LoadTimeAveragedResult(rigidityrange, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge) # need to / 100000.

    #### Load Effective Acceptance
    effective_acceptance_correction_low, effective_acceptance_correction_intermediate = load.LoadEffectiveAcceptance(lowworkpath, rigidityrange, intermediateworkpath)

    #### Load Systematic Error due to Acceptance
    SystematicRelativeErrorAccpetance = load.LoadSystematicErrorFromAcceptance(binmerge) 
    if rigidityrange == "low":
        SysErrorAccTimeIndependentAll = time_averaged_ratio_with_effective/100000 * SystematicRelativeErrorAccpetance[0:8]
    elif rigidityrange == "intermediate":
        SysErrorAccTimeIndependentAll = time_averaged_ratio_with_effective/100000 * SystematicRelativeErrorAccpetance[4:14]

    #### Load Model Prediction
    PbarOverProAslamModel_2426_2480_1BR = LoadExcel.LoadModelFromAslam('ModelPeriod_2426_2480_1BR')
    PbarOverProAslamModel_2480_2506_6BR = LoadExcel.LoadModelFromAslam('ModelPeriod_2480_2506_6BR')
    PbarOverProAslamModel_2481_2520_1BR = LoadExcel.LoadModelFromAslam('ModelPeriod_2481_2520_1BR')

    #### Load TimeDependent MeasuringTime 
    # Note: Rebined = 2; 
    # First index([:,0]) = MeasuringTime in 0.8-1.16 GV, 
    # Second([:,1]) in 1.16-1.51 GV .... 
    # Nineth([:,9]) in 6.47-7.76 GV  
    MeasuringTime_6B = np.load("/rwthfs/rz/cluster/hpcwork/jara0052/sichen/Measuringtime/TimeDependent/TimeDependentMeasuringTime_AllRigidity_6B.npy")


    #### Total value vs Rigidity
    PbarNumber_Merged_6M   = []
    ProtonNumber_Merged_6M = []
    PbarNumber_Merged_6B   = []
    ProtonNumber_Merged_6B = []
    PbarNumber_Merged_3B   = []
    ProtonNumber_Merged_3B = []

    
    #### 1. Main plot loop in Rigidity
    for index in plotrange:
        print("\033[1;31;40mRigidityIndex is: \033[0m" + str((index))) 

        #### Load Fit Result
        draw.LoadFitResult(lowworkpath, intermediateworkpath, trackerpattern, richcut, plotrange, binmerge, time_averaged_ratio, time_averaged_error, effective_acceptance_correction_low, effective_acceptance_correction_intermediate, PosiOverElec, Taiwan_all_6, Ihep_all, Taiwan_all, index, rigidityrange, rigidity_start, SplitTotal_6B, mergestep_6B, modename_6B, SplitTotal_3B, mergestep_3B, modename_3B, SplitTotal_6M, mergestep_6M, modename_6M, pointnumber, SysErrorAccTimeIndependentAll, SystematicRelativeErrorAccpetance)


        #### Draw this analysis
        drawmyresult.drawmyresult_StaErrOnly(lowworkpath, intermediateworkpath, trackerpattern, richcut, rigidityrange, binmerge, yaxisfactor, mode, index, RemovedList_3B)
        #drawmyresult.drawmyresult           (lowworkpath, intermediateworkpath, trackerpattern, richcut, rigidityrange, binmerge, yaxisfactor, mode, index, RemovedList_3B)
        #drawmyresult.drawmyresult_alltime(rigidityrange, lowworkpath, intermediateworkpath, index, binmerge, RemovedList_3B)
        #drawmyresult.drawTimeDependentRatio_vs_Rigditiy(rigidityrange, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, mode, plotrange) ## must run Plot_FullRangeRatio.C, to have CheckResult.root file to hold total error.
        #drawmyresult.draw6Bartel6MonthsCompare(lowworkpath, intermediateworkpath, trackerpattern, richcut, rigidityrange, binmerge, yaxisfactor, mode, index)

        drawFitNumbers.drawFitNumbers(mode, rigidityrange, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, index, time_averaged_proton, time_averaged_antiproton, MeasuringTime_6B, yaxisfactor)
        
        drawChi2.drawChi2(lowworkpath, intermediateworkpath, trackerpattern, richcut, plotrange, rigidityrange, binmerge, index, mode, RemovedList_3B)
        #drawChi2.drawChi2_6Bartel6MonthsCompare(lowworkpath, intermediateworkpath, trackerpattern, richcut, plotrange, rigidityrange, binmerge, index, mode)

        drawBreakDownOfTotalError.drawStatisticalError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, 'absolute')
        drawBreakDownOfTotalError.drawStatisticalError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, 'relative')
        drawBreakDownOfTotalError.drawBreakDownOfTotalSysError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, 'absolute')
        drawBreakDownOfTotalError.drawBreakDownOfTotalSysError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, 'relative')
        drawBreakDownOfTotalError.drawBreakDownOfTotalError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, 'absolute')
        drawBreakDownOfTotalError.drawBreakDownOfTotalError(mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start, 'relative')
        #drawBreakDownOfTotalError.drawErrorPersentage      (mode, rigidityrange, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, rigidity_start)

        #### Draw model prediction
        drawModel.drawModelAslam           (PbarOverProAslamModel_2426_2480_1BR, PbarOverProAslamModel_2480_2506_6BR, PbarOverProAslamModel_2481_2520_1BR, index, rigidityrange, lowworkpath, intermediateworkpath, binmerge, yaxisfactor) 
        drawModel.drawCompareWithModelAslam(PbarOverProAslamModel_2426_2480_1BR, PbarOverProAslamModel_2480_2506_6BR, PbarOverProAslamModel_2481_2520_1BR, index, rigidityrange, lowworkpath, intermediateworkpath, binmerge, yaxisfactor)       
 
        #### Compare with Electron
        #drawEleOverPosRatio.drawEleOverPosRatio(TFile_lepton_result, TFile_lepton_result_Maura, lowworkpath, intermediateworkpath, rigidityrange, Electron, Positron, Electron_error_relative, Positron_error_relative, yaxisfactor, index, ElectronPositronIndex, binmerge, Electron_Maura, Positron_Maura, ElectronError_Maura, PositronError_Maura) 
        #drawCompareWithElectron.drawCompareWithElectron(rigidityrange, lowworkpath, intermediateworkpath, TFile_lepton_result, TFile_lepton_result_Maura, trackerpattern, richcut, plotrange, binmerge, Electron, Positron, Electron_error_relative, Positron_error_relative, yaxisfactor, index, mode, ElectronPositronIndex, Electron_Maura, Positron_Maura, ElectronError_Maura, PositronError_Maura, "Aachen") 
        #drawPosOverEleRatio.drawPosOverEleRatio(TFile_lepton_result, TFile_lepton_result_Maura, lowworkpath, intermediateworkpath, rigidityrange, Electron, Positron, Electron_error_relative, Positron_error_relative, yaxisfactor, index, ElectronPositronIndex, binmerge, PosiOverElec, PosiOverElec_error, Electron_Maura, Positron_Maura, ElectronError_Maura, PositronError_Maura) #(For Error check)
        #drawRelativeResult.drawRelativeResult(rigidityrange, time_averaged_ratio, mode, yaxisfactor, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, index, TFile_lepton_result, TFile_lepton_result_Maura)

        #### Comapre with other groups (Error is statistical only?)
        #if mode == "6Bartels" or mode == "3Bartels":
            #drawCompareWithOtherGroup.drawCompareWithOtherGroup(lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge, Taiwan_all, Ihep_all, Taiwan_all_6, Taiwan_error_all, Taiwan_error_all_6, Ihep_error_all, index, yaxisfactor, mode, rigidityrange) #only 3B and 6B 
        #drawCompareIHEPTaiwanResult.drawCompareIHEPTaiwanResult(rigidityrange, Taiwan_all, Taiwan_error_all, Ihep_all, Ihep_error_all, index, yaxisfactor, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge) # compare 3B
        #drawStatisticalErrorComparision.drawStatisticalErrorComparision(mode, rigidityrange, Taiwan_error_all, Ihep_error_all, Taiwan_error_all_6, Taiwan_all_6, Taiwan_all, Ihep_all, index, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge) 

        #### Obsolete functions
        #drawProtontoAntiproton.drawProtontoAntiproton(rigidityrange) 
        #drawProtontoAntiprotonRelativeResult.drawProtontoAntiprotonRelativeResult(mode, intermediateworkpath, trackerpattern, richcut, binmerge, index) #to be fixed
        #drawAllbins.drawAllbins(mode) #to be fixed
        #drawAllbinsRelative.drawAllbinsRelative(B3_relative) #to be fixed

        ##
        PbarNumber_Merged_6M  .append( sum(draw.M6_PbarNumber_all  ) )
        ProtonNumber_Merged_6M.append( sum(draw.M6_ProtonNumber_all) ) 
        PbarNumber_Merged_6B  .append( sum(draw.B6_PbarNumber_all  ) )
        ProtonNumber_Merged_6B.append( sum(draw.B6_ProtonNumber_all) )
        PbarNumber_Merged_3B  .append( sum(draw.B3_PbarNumber_all  ) )
        ProtonNumber_Merged_3B.append( sum(draw.B3_ProtonNumber_all) )
    
    
    #### Make MakeNumpyForTimeDependent  (save the numpy file in LowResult path) (FIXME: acceptance correction)
    MakeNumpyForTimeDependent.SaveTimeDependentResult("6months" , lowworkpath, intermediateworkpath, "2", "0124", "free", SplitTotal_6M, mergestep_6M, modename_6M, effective_acceptance_correction_low, effective_acceptance_correction_intermediate)
    MakeNumpyForTimeDependent.SaveTimeDependentResult("3Bartels", lowworkpath, intermediateworkpath, "2", "0124", "free", SplitTotal_3B, mergestep_3B, modename_3B, effective_acceptance_correction_low, effective_acceptance_correction_intermediate)
    MakeNumpyForTimeDependent.SaveTimeDependentResult("6Bartels", lowworkpath, intermediateworkpath, "2", "0124", "free", SplitTotal_6B, mergestep_6B, modename_6B, effective_acceptance_correction_low, effective_acceptance_correction_intermediate)
    

    #### 2. Plot for all Rigidity

    PbarNumber_Merged_6M   = np.array( PbarNumber_Merged_6M   )
    ProtonNumber_Merged_6M = np.array( ProtonNumber_Merged_6M )
    PbarNumber_Merged_6B   = np.array( PbarNumber_Merged_6B   )
    ProtonNumber_Merged_6B = np.array( ProtonNumber_Merged_6B )
    PbarNumber_Merged_3B   = np.array( PbarNumber_Merged_3B   )
    ProtonNumber_Merged_3B = np.array( ProtonNumber_Merged_3B )     
    
    #drawTRDTOFCompare.drawTRDTOFCompare(intermediateworkpath, lowworkpath, yaxisfactor, binmerge, index)

    drawFitNumbers.drawBinmergedRatio(time_averaged_antiproton, time_averaged_proton, lowworkpath, intermediateworkpath, binmerge, PbarNumber_Merged_6M, ProtonNumber_Merged_6M, PbarNumber_Merged_6B, ProtonNumber_Merged_6B, PbarNumber_Merged_3B, ProtonNumber_Merged_3B, mode, rigidityrange, trackerpattern, richcut) ## must run Plot_FullRangeRatio.C, to have CheckResult.root file to hold total error.

    drawPbarRatioInAll.drawPbarRatioInAllRigidity(lowworkpath, intermediateworkpath, mode, binmerge)           ## Need to load the numpy file in previous generated file

    drawRelativeErrorvsRigidity.drawRelativeErrorvsRigidity(lowworkpath, intermediateworkpath, mode, binmerge) ## Need to load the numpy file in previous generated file

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-C','--cluster'          , help='which cluster you choose')
    parser.add_argument('-M','--timemode'         , help='which timemode you choose.')
    parser.add_argument('-B','--binmerge'         , type=str, help='how many bins you used for each fit.')
    parser.add_argument('-A','--ratiopattern'     , type=str, help='antiproton to proton ratio or proton to antiproton ratio, only valid for relative plots.')
    parser.add_argument('-R','--rigidityrange'    , help='which rigidityrange you choose:intermediate or low')

    parser.add_argument('-P','--trackerpattern'   , type=str, help='which trackerpattern you choose. Only for intermediate range')
    parser.add_argument('-I','--richcut', type=str, help='what richcut you choose. Only for intermediate range')
    arguments = parser.parse_args()


    if arguments.cluster == "JUAMS":
        intermediateworkpath = os.getenv('JUAMSINTERMEDIATEENERGYDIR')
        print("in progress")
    elif arguments.cluster == "HPC":
        intermediateworkpath = os.getenv('HPCINTERMEDIATEDIR')
        lowworkpath          = os.getenv('HPCLOWENERGYDIR')


    trackerpattern = arguments.trackerpattern
    richcut        = arguments.richcut
    binmerge       = arguments.binmerge
    ratiopattern   = arguments.ratiopattern
    mode           = arguments.timemode
    rigidityrange  = arguments.rigidityrange


    yaxisfactor = 1.1

    SplitTotal_3B = TimeDependentTimeRange.SplitTotal_3B
    mergestep_3B  = TimeDependentTimeRange.mergestep_3B
    modename_3B   = "3BartalRotation"

    SplitTotal_6B = TimeDependentTimeRange.SplitTotal_6B
    mergestep_6B  = TimeDependentTimeRange.mergestep_6B
    modename_6B   = "6BartalRotation"

    SplitTotal_6M = TimeDependentTimeRange.SplitTotal_6M
    mergestep_6M  = TimeDependentTimeRange.mergestep_6M
    modename_6M   = "6months"

    if rigidityrange == "low":
        rigidity_start = 1 # correspondent to 1.16
        #rigidity_start = 3 # correspondent to 1.51
        #rigidity_end = 15 # correspondent to 5.37
        rigidity_end = 17
    elif rigidityrange == "intermediate":
        rigidity_start = 9 # correspondent to 2.97
        rigidity_end = 29 # correspondent to 18.0

    pointnumber = int((rigidity_end-rigidity_start)/int(binmerge))

    plotrange = range(rigidity_start, rigidity_end, int(binmerge))
    #plotrange = [11,13,15,17,19,21,23] # 3.64, 4.4, 5.37, 6.47, 7.76, 9.26, 11.0

    RemovedList_3B = [15, 38, 43]

    main()



