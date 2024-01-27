#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
import matplotlib.dates as md
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import argparse
import binning
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import numpy.polynomial.polynomial as poly
from scipy import optimize
import uproot

def LoadElectronPositronResult(GetYindex):
    time_array = binning.Bartals1Unixtime_edge
    lepton_result = TFile("$MY_ANALYSIS/ReferenceFiles/TimeDependentAnalysis/LeptonFluxes_Zimmermann_04_09_2018_2D.root")

    ## 52 points. GetX/Y index:0-51. result is from index 2. first two is under flow.  Attention: index of GetXaxis()->GetBinLowEdge is different. so if you want to know the binning, please look at PRL paper.
    ## binindex:
    #(-6). 1.01-1.22  GetX/Y is 2
    #(-5). 1.22-1.46  GetX/Y is 3
    #(-4). 1.46-1.72  GetX/Y is 4
    #(-3). 1.72-2.00  GetX/Y is 5
    #(-2). 2.00-2.31  GetX/Y is 6
    #(-1). 2.31-2.65  GetX/Y is 7
    #(0). 2.65-3.00  GetX/Y is 8
    #(1). 3.00-3.36 index 7-8.  GetX/Y is 9(since underflows).
    #(2). 3.73-4.12. index 9-10. GetX/Y is 11(since underflows).
    #(3). 4.54-5.00. index 11-12. GetX/Y is 13.
    #(4). 5.49-6.00. index 13-14. GetX/Y is 15.
    #(5). 6.54-7.10. index 15-16. GetX/Y is 17.
    #(6). 7.10-7.69. index 16-17. GetX/Y is 18.
    #(7). 8.30-8.95. index 18-19. GetX/Y is 20.
    #(8). 9.62-10.32. index 20-21. GetX/Y is 22.
    #(9). 11.04-11.8. index 22-23. GetX/Y is 24.
    #(10).11.8-12.59. index 23-24. GetX/Y is 25.
    #(11).13.41-14.25. index 25-26. GetX/Y is 27.
    #(12).16.05-17.0. index 28-29. GetX/Y is 30.

    PosiOverElec = []
    PosiOverElec_error = []
    Electron=[]
    Positron=[]
    Electron_error_relative=[]
    Positron_error_relative=[]
    for i in range(0,88):   ## Paper published for 79 Bartels. Nico's last result is for 88 Bartels.  # 86 Bartels have results. [46,47] are empty.
        PosiOverElec.append(lepton_result.Get("SxRatioTools_PosiOverElec_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))).Get("grPosiOverElec_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))+"_PosiOverElec").GetY()[GetYindex])
        PosiOverElec_error.append(lepton_result.Get("SxRatioTools_PosiOverElec_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))).Get("grPosiOverElec_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))+"_PosiOverElec").GetErrorY(GetYindex))
        Electron.append(lepton_result.Get("SxFluxTools_Electron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))).Get("grElectron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))+"_FinalFlux").GetY()[GetYindex])
        Electron_error_relative.append(lepton_result.Get("SxFluxTools_Electron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))).Get("grElectron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))+"_RelErrorStat").GetY()[GetYindex])
        Positron.append(lepton_result.Get("SxFluxTools_Positron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))).Get("grPositron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))+"_FinalFlux").GetY()[GetYindex])
        Positron_error_relative.append(lepton_result.Get("SxFluxTools_Positron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))).Get("grPositron_"+str(int(time_array[i]))+"_"+str(int(time_array[i+1]))+"_RelErrorStat").GetY()[GetYindex])
    #E_ave = lepton_result.Get("SxRatioTools_PosiOverElec_RWTH").Get("grPosiOverElec_RWTH_PosiOverElec").GetY()[GetYindex]
    #E_ave_error = lepton_result.Get("SxRatioTools_PosiOverElec_RWTH").Get("grPosiOverElec_RWTH_PosiOverElec").GetErrorY(GetYindex)

    PosiOverElec = np.array(PosiOverElec)
    PosiOverElec_error = np.array(PosiOverElec_error)
    Electron = np.array(Electron)
    Positron = np.array(Positron)
    Electron_error_relative = np.array(Electron_error_relative)
    Positron_error_relative = np.array(Positron_error_relative)
    #PosiOverElec_relative_error = np.sqrt( (PosiOverElec_error/E_ave)**2 + (PosiOverElec*E_ave_error/E_ave**2)**2 )
    #PosiOverElec_relative_error = PosiOverElec_error/E_ave

    return lepton_result, PosiOverElec, PosiOverElec_error, Electron, Electron_error_relative, Positron, Positron_error_relative


def LoadMauraElectronPositronResultMonthly(GetYindex):
    lepton_result_Maura = TFile("$MY_ANALYSIS/ReferenceFiles/TimeDependentAnalysis/MauraEleOverPos_Nov2017/ele_pos_vs_time_TGraph.root")

    ## binindex:
    #(-6). 1.01-1.22  GetYindex is 2
    #(-5). 1.22-1.46  GetYindex is 3
    #(-4). 1.46-1.72  GetYindex is 4
    #(-3). 1.72-2.00  GetYindex is 5

    g_positronflux = lepton_result_Maura.Get("pos_maura_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[GetYindex+1])) + "_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[GetYindex+2])) + ";1")
    g_electronflux = lepton_result_Maura.Get("ele_maura_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[GetYindex+1])) + "_" + str("{:.2f}".format(binning.LeptonLowEnergyBin[GetYindex+2])) + ";1")

    #g_positronflux = lepton_result_Maura.Get("pos_maura_" + "1.72" + "_" + "2.00" + ";1")
    #g_electronflux = lepton_result_Maura.Get("ele_maura_" + "1.72" + "_" + "2.00" + ";1")

    Electron_Maura=[]
    Positron_Maura=[]
    ElectronError_Maura=[]
    PositronError_Maura=[]

    for i in range(g_positronflux.GetN()):
        Electron_Maura.append(g_electronflux.GetY()[i])
        Positron_Maura.append(g_positronflux.GetY()[i])
        ElectronError_Maura.append(g_electronflux.GetErrorY(i))
        PositronError_Maura.append(g_positronflux.GetErrorY(i))

    return lepton_result_Maura, Electron_Maura, ElectronError_Maura, Positron_Maura, PositronError_Maura

def LoadOtherGroupResult():   
    ## 1.IHEP result
    ## For old version,
    ## 1.0-15.3 GetBinContent(timestampt,index) index:2-16. 15 points in total.
    # 1.0-1.51: index2 (3 bins)
    # 1.51-1.92 index3 (2 bins)
    #1.92-2.97 index4 (4 bins)
    #2.97-3.64 index5 (2 bins)
    #3.64-4.43 index6 (2 bins)
    #4.43-5.37 index7 (2 bins)
    #5.37-5.9 index8 (1 bin)
    #5.9-6.47 index9 (1 bin)
    #6.47-7.09 index10 (1 bin)
    #7.09-7.76 index11 (1 bin)
    #7.76-8.48 index12 (1 bin)
    #8.48-9.26 index13 (1 bin)
    #9.26-11.0 index14 (2 bins)
    #11.0-13.0 index15 (2 bins)
    #13.0-15.3 index16 (2 bins)

    #### For new version(to oct 2019): 1.51-33.5 GV, index:5-20.  16 points in total.
    # 0.7-0.8 index1
    # 0.8-0.9 index2
    # 0.9-1.0 index3
    # 1.0-1.51 index4
    # 1.51-1.92: index5 (2 bins) From this bin, start to have value.
    # 1.92-2.97: index6 (4 bins)
    # 2.97-3.64: index7 (2 bins)
    # 3.64-4.43: index8 (2 bins)
    # 4.43-5.37: index9 (2 bins)
    # 5.37-6.47: index10 (2 bins)
    # 6.47-7.76: index11 (2 bins)
    # 7.76-9.26: index12 (2 bins)
    # 9.26-11: index13 (2 bins)
    # 11-13: index14 (2 bins)
    # 13-15.3: index15 (2 bins)
    # 15.3-18: index16 (2 bins)
    # 18-21.1: index17 (2 bins)
    # 21.1-24.7: index18 (2 bins)
    # 24.7-28.8: index19 (2 bins)
    # 28.8-33.5: index20 (2 bins)

    #IHEP_result = TFile("$MY_ANALYSIS/ReferenceFiles/TimeDependentAnalysis/AntiprotonFlux_AntiPBin_3Months_IHEP_old.root")
    IHEP_result = TFile("$MY_ANALYSIS/ReferenceFiles/TimeDependentAnalysis/Antiproton_pbin_3Months_IHEP_Oct_2019.root")  #from min:06/11/2011 to max:11/14/2019
    hh_Ratio_StatErr = IHEP_result.Get("hh_Ratio_StatErr")
    Ihep_all = [[],[],[],[],[]]
    Ihep_error_all = [[],[],[],[],[]]
    for column in range(5,21): #5-20: 5:1.5-1.92, 6:1.92-2.97.......... 20:28.8-33.5
        Ihep_all.append([])
        Ihep_error_all.append([])
        for row in range(1,39): ## extension Oct 2019:1-38 points.
            Ihep_all[column].append(hh_Ratio_StatErr.GetBinContent(row,column))
            Ihep_error_all[column].append(hh_Ratio_StatErr.GetBinError(row,column))
    del Ihep_all[0]
    del Ihep_error_all[0]
    del Ihep_all[0]
    del Ihep_error_all[0]
    del Ihep_all[0]
    del Ihep_error_all[0]
    del Ihep_all[0]
    del Ihep_error_all[0]
    del Ihep_all[0]
    del Ihep_error_all[0]

    #### 2.1 Taiwan group result (Fist bartel rotation removed, last bin max is 11.14.2019)
    #### To oct 2019: 1.0-45.1GV
    # 1.0-1.51: index1 (3bins)
    # 1.51-1.92: index2 (2bins)
    # 1.92-2.97: index3 (4bins)
    # 2.97-3.64: index4 (2bins)
    # 3.64-4.43: index5 (2bins)
    # 4.43-5.37: index6 (2bins)
    # 5.37-6.47: index7 (2bins)
    # 6.47-7.76: index8 (2bins)
    # 7.76-9.26: index9 (2bins)
    # 9.26-11.0: index10 (2bins)
    # 11.0-13.0: index11 (2bins)
    # 13.0-15.3: index12 (2bins)
    # 15.3-18.0: index13 (2bins)
    # 18.0-21.0: index14 (2bins)
    # 21.0-24.7: index15 (2bins)
    # 24.7-28.8: index16 (2bins)
    # 28.8-33.5: index17 (2bins)
    # 33.5-38.9: index18 (2bins)
    # 38.9-45.1: index19 (2bins)
    Taiwan_all = [[]]
    Taiwan_error_all = [[]]

    Taiwan_result = TFile("$MY_ANALYSIS/ReferenceFiles/TimeDependentAnalysis/apflux_tme3M.root")
    hTAppStat_3 = Taiwan_result.Get("hTAppStat")
    hTAppTotl_3 = Taiwan_result.Get("hTAppTotl")
    for column in range(1,20):
        Taiwan_all.append([])
        Taiwan_error_all.append([])
        for row in range(1,39): ## extension: 1-38 points.
            Taiwan_all[column].append(hTAppStat_3.GetBinContent(row,column))
            Taiwan_error_all[column].append(hTAppStat_3.GetBinError(row,column))
    del Taiwan_all[0]
    del Taiwan_error_all[0]
    Taiwan_all = np.array(Taiwan_all)
    Taiwan_error_all = np.array(Taiwan_error_all)

    #### 2.2 Taiwan group result (Fist bartel rotation included, up to 10.18.2019)
    #### To oct 2019: 1.0-45.1GV
    # 1.0-1.92: index1 (5bins)
    # 1.92-2.97: index2 (4bins)
    # 2.97-3.64: index3 (2bins)
    # 3.64-4.43: index4 (2bins)
    # 4.43-5.37: index5 (2bins)
    # 5.37-6.47: index6 (2bins)
    # 6.47-7.76: index7 (2bins)
    # 7.76-9.26: index8 (2bins)
    # 9.26-11.0: index9 (2bins)
    # 11.0-13.0: index10 (2bins)
    # 13.0-15.3: index11 (2bins)
    # 15.3-18.0: index12 (2bins)
    # 18.0-21.0: index13 (2bins)
    # 21.0-24.7: index14 (2bins)
    # 24.7-28.8: index15 (2bins)
    # 28.8-33.5: index16 (2bins)
    # 33.5-38.9: index17 (2bins)
    # 38.9-45.1: index18 (2bins)
    Taiwan_result_6 = TFile("$MY_ANALYSIS/ReferenceFiles/TimeDependentAnalysis/apflux_tme6M.root")
    hTAppStat_6 = Taiwan_result_6.Get("hTAppStat")
    hTAppTotl_6 = Taiwan_result_6.Get("hTAppTotl")
    Taiwan_all_6 = [[]]
    Taiwan_error_all_6 = [[]]

    for column in range(1,19):
        Taiwan_all_6.append([])
        Taiwan_error_all_6.append([])
        for row in range(1,20): ## extension: 1-19 points.
            Taiwan_all_6[column].append(hTAppStat_6.GetBinContent(row,column))
            Taiwan_error_all_6[column].append(hTAppStat_6.GetBinError(row,column))
    del Taiwan_all_6[0]
    del Taiwan_error_all_6[0]

    Ihep_all = np.array(Ihep_all)
    Ihep_error_all = np.array(Ihep_error_all)
    Taiwan_all = np.array(Taiwan_all)
    Taiwan_error_all = np.array(Taiwan_error_all)
    Taiwan_all_6 = np.array(Taiwan_all_6)
    Taiwan_error_all_6 = np.array(Taiwan_error_all_6)

    return Ihep_all, Ihep_error_all, Taiwan_all, Taiwan_error_all, Taiwan_all_6, Taiwan_error_all_6

 
def LoadTimeAveragedResult(rigidityrange, lowworkpath, intermediateworkpath, trackerpattern, richcut, binmerge):   
    time_averaged_ratio                = np.array([])
    time_averaged_ratio_with_effective = np.array([])
    time_averaged_error                = np.array([])
    time_averaged_proton               = np.array([])
    time_averaged_antiproton           = np.array([])

    if rigidityrange == "low":
        f_averaged = TFile(lowworkpath + "/totalall/Time_Averaged_ratio_Low/binmerge" + str(binmerge) + "/plots/Ratio_pass7.8.root")
        g_ratio                = f_averaged.Get("ratio_tof_TRDeff_0.94_TOFeff_0.95")
        g_ratio_with_effective = f_averaged.Get("ratio_tof_with_effective_TRDeff_0.94_TOFeff_0.95")
        g_error                = f_averaged.Get("g_error_TRDeff_0.94_TOFeff_0.95")
        g_proton               = f_averaged.Get("g_proton_number_TRDeff_0.94_TOFeff_0.95")
        g_antiproton           = f_averaged.Get("g_antiproton_number_TRDeff_0.94_TOFeff_0.95")
    elif rigidityrange == "intermediate":
        f_averaged = TFile(intermediateworkpath + "/total/Time_Averaged_ratio/binmerge" + str(binmerge) + "/intermediate_"+ str(trackerpattern) + "_" + str(richcut) + "_" + "pass7.8" + "binmerge" + str(binmerge) + ".root")
        g_ratio                = f_averaged.Get("g_ratio_TRDeff_0.99")
        g_ratio_with_effective = f_averaged.Get("g_ratio_with_effective_acceptance_TRDeff_0.99")
        g_error                = f_averaged.Get("g_error_TRDeff_0.99")
        g_proton               = f_averaged.Get("g_proton_TRDeff_0.99")
        g_antiproton           = f_averaged.Get("g_antiproton_TRDeff_0.99")

    for i in range(g_ratio_with_effective.GetN()):
        time_averaged_ratio                = np.append(time_averaged_ratio               , g_ratio.GetY()[i])
        time_averaged_ratio_with_effective = np.append(time_averaged_ratio_with_effective, g_ratio_with_effective.GetY()[i])
        time_averaged_error                = np.append(time_averaged_error               , g_error.GetY()[i])
        time_averaged_proton               = np.append(time_averaged_proton              , g_proton.GetY()[i])
        time_averaged_antiproton           = np.append(time_averaged_antiproton          , g_antiproton.GetY()[i])

    return time_averaged_ratio, time_averaged_ratio_with_effective, time_averaged_error, time_averaged_proton, time_averaged_antiproton


def LoadEffectiveAcceptance(lowworkpath, rigidityrange, intermediateworkpath):
    EffectiveAcceptance_antipro_low          = TFile(lowworkpath + "/totalall" + "/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root")
    EffectiveAcceptance_proton_low           = TFile(lowworkpath + "/totalall" + "/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root")
    EffectiveAcceptance_antipro_intermediate = TFile(lowworkpath + "/totalall" + "/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root")
    EffectiveAcceptance_proton_intermediate  = TFile(lowworkpath + "/totalall" + "/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root")

    Acceptance_antiproton_low          = EffectiveAcceptance_antipro_low.Get("QualityCuts/effectiveAcceptanceAfterAllCuts") ## 62 points in total,0.8-1130, GetX()[0]=0.9, GetX()[16]=5.635, GetX()[61]=976,fitst and last points are 0, since MC generated R range.   # Acceptance_antiproton_low.GetN()-1 = 61
    Acceptance_proton_low              = EffectiveAcceptance_proton_low.Get("QualityCuts/effectiveAcceptanceAfterAllCuts")
    Acceptance_antiproton_intermediate = EffectiveAcceptance_antipro_intermediate.Get("QualityCuts/effectiveAcceptanceAfterAllCuts") 
    Acceptance_proton_intermediate     = EffectiveAcceptance_proton_intermediate.Get("QualityCuts/effectiveAcceptanceAfterAllCuts")

    effective_acceptance_correction_low = np.array([])
    effective_acceptance_correction_low = np.append(effective_acceptance_correction_low, 0.0)   ## first point gives 0. ()
    effective_acceptance_correction_intermediate = np.array([])
    effective_acceptance_correction_intermediate = np.append(effective_acceptance_correction_intermediate, 0.0)   

    #for i in range(1,17):   ## from index 1 and 60, index 1 and 61 are manually given as 0.  
    #for i in range(1, Acceptance_antiproton.GetN()-1):
    for i in range(1,30):   ## GetX()[0] = 0.8-1.0 (0/0problem), GetX()[1] = 1.0-1.16, GetX()[2] = 1.16-1.33
        effective_acceptance_correction_low          = np.append(effective_acceptance_correction_low         , Acceptance_proton_low.GetY()[i]/Acceptance_antiproton_low.GetY()[i])
        effective_acceptance_correction_intermediate = np.append(effective_acceptance_correction_intermediate, Acceptance_proton_intermediate.GetY()[i]/Acceptance_antiproton_intermediate.GetY()[i])
        #print('i: ' + str(i) + ", " + "X value:" + str(Acceptance_proton.GetX()[i]) )
    #effective_acceptance_correction = np.append(effective_acceptance_correction, 0.0) ## last point gives 0.



    ''' index inforomation:
    effective_acceptance_correction 0         （0.8-1.0）
    effective_acceptance_correction 1 : 1.08   (1.0-1.16)
    effective_acceptance_correction 2 : 1.245
    effective_acceptance_correction 3 : 1.42 
    effective_acceptance_correction 4 : 1.6099999999999999
    effective_acceptance_correction 5 : 1.815
    effective_acceptance_correction 6 : 2.035
    '''

    return effective_acceptance_correction_low, effective_acceptance_correction_intermediate 


def LoadSystematicErrorFromAcceptance(binmerge):
    f_eff_B1220_antiproton         = TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_antipr.pl1ph.021000.qgsp_bic_ams_7.8_all.root");
    f_eff_B1220_antiproton_minus10 = TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_antipr.pl1ph.021000.qgsp_bic_ams.minus10_7.8_all.root");
    f_eff_B1220_antiproton_plus10  = TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_antipr.pl1ph.021000.qgsp_bic_ams.plus10_7.8_all.root");
    f_eff_B1220_proton             = TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_pr.pl1ph.021000_7.8_all.root");

    Acceptance_B1220_antiproton         = f_eff_B1220_antiproton        .Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    Acceptance_B1220_antiproton_minus10 = f_eff_B1220_antiproton_minus10.Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    Acceptance_B1220_antiproton_plus10  = f_eff_B1220_antiproton_plus10 .Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    Acceptance_B1220_proton             = f_eff_B1220_proton            .Get("QualityCuts/effectiveAcceptanceAfterAllCuts");

    ## Remove last point to avoid NaN
    Acceptance_B1220_antiproton        .RemovePoint(61);
    Acceptance_B1220_antiproton_minus10.RemovePoint(61);
    Acceptance_B1220_antiproton_plus10 .RemovePoint(61);
    Acceptance_B1220_proton            .RemovePoint(61);

    Acceptance_B1220_antiproton        .RemovePoint(60);
    Acceptance_B1220_antiproton_minus10.RemovePoint(60);
    Acceptance_B1220_antiproton_plus10 .RemovePoint(60);
    Acceptance_B1220_proton            .RemovePoint(60);

    Acceptance_B1220_antiproton        .RemovePoint(0);
    Acceptance_B1220_antiproton_minus10.RemovePoint(0);
    Acceptance_B1220_antiproton_plus10 .RemovePoint(0);
    Acceptance_B1220_proton            .RemovePoint(0);

    ## Fill EffectiveAcceptance ratios in TGraphErrors and calculate the Systematic Uncertainties.
    number = Acceptance_B1220_antiproton.GetN();

    Effective_Acceptance_ratio         = TGraphErrors(number-1);
    Effective_Acceptance_ratio_minus10 = TGraphErrors(number-1);
    Effective_Acceptance_ratio_plus10  = TGraphErrors(number-1);
    g_SysUncertaintyRel_ACCratio       = TGraph(number-1);

    SystematicRelativeErrorAccpetance = []

    for p in range(number-1):
        ##(1). Fill EffectiveAcceptance ratios
        ## value taken start from 1, because 0.8-1.0GV is not used in this analysis.
        Aa = Acceptance_B1220_antiproton                .GetY()[p+1];
        Aa_minus10 = Acceptance_B1220_antiproton_minus10.GetY()[p+1];
        Aa_plus10 = Acceptance_B1220_antiproton_plus10  .GetY()[p+1];
        Ap = Acceptance_B1220_proton                    .GetY()[p+1];
        ea = Acceptance_B1220_antiproton                .GetErrorY(p+1);
        ea_minus10 = Acceptance_B1220_antiproton_minus10.GetErrorY(p+1);
        ea_plus10 = Acceptance_B1220_antiproton_plus10  .GetErrorY(p+1);
        ep = Acceptance_B1220_proton                    .GetErrorY(p+1);
        x  = Acceptance_B1220_antiproton                .GetX()[p+1];

        Effective_Acceptance_ratio        .SetPoint     (p, x, Ap/Aa)
        Effective_Acceptance_ratio        .SetPointError(p, 0, np.sqrt( ep**2/Aa**2 + Ap**2/Aa**4*ea**2))
        Effective_Acceptance_ratio_minus10.SetPoint     (p, x, Ap/Aa_minus10)
        Effective_Acceptance_ratio_minus10.SetPointError(p, 0, np.sqrt(ep**2/Aa_minus10**2 + Ap**2/Aa_minus10**4*ea_minus10**2))
        Effective_Acceptance_ratio_plus10 .SetPoint     (p, x, Ap/Aa_plus10)
        Effective_Acceptance_ratio_plus10 .SetPointError(p, 0, np.sqrt( ep**2/Aa_plus10**2 + Ap**2/Aa_plus10**4*ea_plus10**2 ))

        ##(2). calculate the Systematic Uncertainties due to ACC.
        SystematicRelativeErrorAccpetance.append( Ap/Aa**2*((Aa-Aa_plus10) + (Aa_minus10-Aa))/2 )

    SystematicRelativeErrorAccpetance_merged = []
    for q in range(0, len(SystematicRelativeErrorAccpetance), int(binmerge)):
        SystematicRelativeErrorAccpetance_merged.append( (SystematicRelativeErrorAccpetance[q] + SystematicRelativeErrorAccpetance[q+1])/2 )

    return SystematicRelativeErrorAccpetance_merged


