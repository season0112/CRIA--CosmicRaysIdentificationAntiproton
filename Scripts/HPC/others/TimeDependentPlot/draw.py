#!/usr/bin/env python

#### In this module, we have:
# PRL_func, Fermi_func for fit functions.
# LoadFitResult

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
from uncertainties import ufloat
import model


def PRL_func(x, x0, c, miu, kt):
    x=x/1305417600
    miu=1.0515005941393774
    kt=0.012513472713317396
    return x0 * (1+ c/(np.exp(-(x-miu)/kt)+1))


def Fermi_func(x, x0, c, miu, kt):
    x=x/1305417600
    return x0 * (1+ c/(np.exp(-(x-miu)/kt)+1))


def LoadFitResult(lowworkpath, intermediateworkpath, trackerpattern, richcut, plotrange, binmerge, time_averaged_ratio, time_averaged_error, effective_acceptance_correction_low, effective_acceptance_correction_intermediate, PosiOverElec, Taiwan_all_6, Ihep_all, Taiwan_all, index, rigidityrange, rigidity_start, SplitTotal_6B, mergestep_6B, modename_6B, SplitTotal_3B, mergestep_3B, modename_3B, SplitTotal_6M, mergestep_6M, modename_6M, pointnumber, ErrorAll_Acc, SystematicRelativeErrorAccpetance):

    global dates_1Bartel_unixsecond, dates_1Bartel, dates_1Bartel_unixsecond_Aslam, dates_1Bartel_Aslam   

    global datesM6, M6result, ErrorM6, Error_Sys6M, Error_SysAcc6M, Error_TotalSys6M, Error_TotalTimeDependent6M, Error_Total6M, M6Chi2, datesM6_unixsecond
    global dates6B, B6result, Error6b, Error_Sys6B, Error_SysAcc6B, Error_TotalSys6B, Error_TotalTimeDependent6B, Error_Total6B, B6Chi2, datesB6_unixsecond, dates6B_TaiwanRange
    global dates,   B3result, Error3b, Error_Sys3B, Error_SysAcc3B, Error_TotalSys3B, Error_TotalTimeDependent3B, Error_Total3B, B3Chi2, dates_unixsecond  , dates_IhepRange    , dates_TaiwanRange
    global Error_SysAccAveraged


    global RelativeStatisticalError_6B    , RelativeStatisticalError_3B    , RelativeStatisticalError_6M
    global RelativeSystematicError_6B     , RelativeSystematicError_3B     , RelativeSystematicError_6M
    global RelativeStatisticalError_6B_all, RelativeStatisticalError_3B_all, RelativeStatisticalError_6M_all 
    global RelativeSystematicError_6B_all , RelativeSystematicError_3B_all , RelativeSystematicError_6M_all

    global B3_relative,  B3relative_error, B3_relative_all, B3relative_error_all
    global B6_relative,  B6relative_error, B6_relative_all, B6relative_error_all
    global M6_relative,  M6relative_error, M6_relative_all, M6relative_error_all

    global B3result_all, Error3b_all
    global B6result_all, Error6b_all
    global M6result_all, Error6m_all

    global B3_PbarNumber,   B6_PbarNumber,   M6_PbarNumber  , B3_PbarNumber_all  , B6_PbarNumber_all  , M6_PbarNumber_all
    global B3_ProtonNumber, B6_ProtonNumber, M6_ProtonNumber, B3_ProtonNumber_all, B6_ProtonNumber_all, M6_ProtonNumber_all

    global B3_P_Anti_Error, B3_relative_proton_to_antiproton
    global lepton_result

    if index ==  plotrange[0]:
        B3_relative_all      = []
        B3relative_error_all = []    
        B6_relative_all      = []
        B6relative_error_all = []
        M6_relative_all      = []
        M6relative_error_all = []

        B3result_all = []
        Error3b_all  = []
        B6result_all = []
        Error6b_all  = []
        M6result_all = []
        Error6m_all  = []

        RelativeStatisticalError_6B_all = []
        RelativeStatisticalError_3B_all = []
        RelativeStatisticalError_6M_all = []
        RelativeSystematicError_6B_all = []
        RelativeSystematicError_3B_all = []
        RelativeSystematicError_6M_all = []


    B3_PbarNumber_all   = [] 
    B6_PbarNumber_all   = []  
    M6_PbarNumber_all   = []
    B3_ProtonNumber_all = []
    B6_ProtonNumber_all = []
    M6_ProtonNumber_all = []



    #### Load Fit Result in 6 months

    ## fit result
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as f6m:
            lines6m=f6m.readlines()
        M6result = np.array(list(map(lambda s: float(s.strip()), lines6m)))
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + "_6months.txt") as f6m:
            lines6m=f6m.readlines()
        M6result = np.array(list(map(lambda s: float(s.strip()), lines6m)))

    ## fit result error
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as e6m:
            error6m = e6m.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + "_6months.txt") as e6m:
            error6m = e6m.readlines()
    ErrorM6 = np.array(list(map(lambda s: float(s.strip()), error6m)))

    ## fit chi2
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/fit_chi2dof_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as chi2_6m:
            chi2dof6m = chi2_6m.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_chi2dof_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as chi2_6m:
            chi2dof6m = chi2_6m.readlines()
    M6Chi2 = np.array(list(map(lambda s: float(s.strip()), chi2dof6m)))

    ## relative error
    if rigidityrange == "low":
        #M6_relative     = M6result / np.mean(M6result)
        M6_relative      = M6result / time_averaged_ratio[int((index-1)/int(binmerge))]
        M6relative_error = np.sqrt( ( ErrorM6/time_averaged_ratio[int((index-1)/int(binmerge))] )**2 + ( M6result * time_averaged_error[int((index-1)/int(binmerge))]/time_averaged_ratio[int((index-1)/int(binmerge))]**2 )**2 )
    elif rigidityrange == "intermediate":
        #M6_relative     = M6result / np.mean(M6result)        
        M6_relative      = M6result / time_averaged_ratio[int((index-9)/int(binmerge))]
        M6relative_error = np.sqrt( ( ErrorM6/time_averaged_ratio[int((index-9)/int(binmerge))] )**2 + ( M6result * time_averaged_error[int((index-9)/int(binmerge))]/time_averaged_ratio[int((index-9)/int(binmerge))]**2 )**2 )

    ## Antiproton number
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/antiprotonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as pbarNum_6m:
            PbarNum_6m = pbarNum_6m.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/antiprotonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as pbarNum_6m:
            PbarNum_6m = pbarNum_6m.readlines()
    M6_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_6m)))

    ## Proton number
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/protonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as protonNum_6m:
            ProtonNum_6m = protonNum_6m.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/protonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6months.txt") as protonNum_6m:
            ProtonNum_6m = protonNum_6m.readlines()
    M6_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_6m)))


    #### Load Fit Result in 6 Bartel Rotations

    ## fit result
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as f6b:
            lines6b=f6b.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total/" + "results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + "_6BartalRotation.txt") as f6b:
            lines6b=f6b.readlines()
    B6result = np.array(list(map(lambda s: float(s.strip()), lines6b)))

    ## fit result error
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as e6b:
            error6b = e6b.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index++int(binmerge)]) + "_6BartalRotation.txt") as e6b:
            error6b = e6b.readlines()
    Error6b = np.array(list(map(lambda s: float(s.strip()), error6b)))

    ## fit chi2
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/fit_chi2dof_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as chi2_6b:
            chi2dof6b = chi2_6b.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_chi2dof_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as chi2_6b:
            chi2dof6b = chi2_6b.readlines()
    B6Chi2 = np.array(list(map(lambda s: float(s.strip()), chi2dof6b)))

    ## relative error
    if rigidityrange == "low":
        #B6_relative     = B6result / np.mean(B6result)
        B6_relative      = B6result / (time_averaged_ratio[int((index-1)/int(binmerge))]/10**5)
        B6relative_error = np.sqrt( ( Error6b / (time_averaged_ratio[int((index-1)/int(binmerge))]/10**5) )**2 + ( B6result * (time_averaged_error[int((index-1)/int(binmerge))]/10**5) / (time_averaged_ratio[int((index-1)/int(binmerge))]/10**5)**2 )**2 )
    elif rigidityrange == "intermediate":
        #B6_relative     = B6result / np.mean(B6result)
        B6_relative      = B6result / time_averaged_ratio[int((index-9)/int(binmerge))]
        B6relative_error = np.sqrt( ( Error6b/time_averaged_ratio[int((index-9)/int(binmerge))] )**2 + ( B6result * time_averaged_error[int((index-9)/int(binmerge))]/time_averaged_ratio[int((index-9)/int(binmerge))]**2 )**2 )

    ## Antiproton number
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/antiprotonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as pbarNum_6B:
            PbarNum_6B = pbarNum_6B.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/antiprotonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as pbarNum_6B:
            PbarNum_6B = pbarNum_6B.readlines()
    B6_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_6B)))

    ## Proton number
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/protonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as protonNum_6B:
            ProtonNum_6B = protonNum_6B.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/protonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_6BartalRotation.txt") as protonNum_6B:
            ProtonNum_6B = protonNum_6B.readlines()
    B6_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_6B)))


    #### Load Fit Result in 3 Bartel Rotations

    ## fit result
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall/" + "results/fit_results_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as f3b:
            lines3b=f3b.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total/" + "results/fit_results_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as f3b:
            lines3b=f3b.readlines()
    B3result = np.array(list(map(lambda s: float(s.strip()), lines3b)))

    ## fit result error
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/fit_results_error_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as e3b:
            error3b = e3b.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_results_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as e3b:
            error3b = e3b.readlines()
    Error3b = np.array(list(map(lambda s: float(s.strip()), error3b)))

    ## fit chi2
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/fit_chi2dof_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as chi2_3b:
            chi2dof3b = chi2_3b.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_chi2dof_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as chi2_3b:
            chi2dof3b = chi2_3b.readlines()
    B3Chi2 = np.array(list(map(lambda s: float(s.strip()), chi2dof3b)))

    ## relative error
    if rigidityrange == "low":
        #B3_relative     = B3result / np.mean(B3result)
        B3_relative      = B3result / time_averaged_ratio[int((index-1)/int(binmerge))]
        B3relative_error = np.sqrt( ( Error3b/time_averaged_ratio[int((index-1)/int(binmerge))] )**2 + ( B3result * time_averaged_error[int((index-1)/int(binmerge))]/time_averaged_ratio[int((index-1)/int(binmerge))]**2 )**2 )
    elif rigidityrange == "intermediate":
        #B3_relative     = B3result / np.mean(B6result)
        B3_relative      = B3result / time_averaged_ratio[int((index-9)/int(binmerge))]
        B3relative_error = np.sqrt( ( Error3b/time_averaged_ratio[int((index-9)/int(binmerge))] )**2 + ( B3result * time_averaged_error[int((index-9)/int(binmerge))]/time_averaged_ratio[int((index-9)/int(binmerge))]**2 )**2 )

    # Proton to Antiproton error
    if rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/fit_results_proton_antiproton_error_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as b3_proton_antiproton_error:
            b3_p_anti_error = b3_proton_antiproton_error.readlines()
        B3_P_Anti_Error = np.array(list(map(lambda s: float(s.strip()), b3_p_anti_error)))
        B3_relative_proton_to_antiproton = B3_P_Anti_Error/(1./time_averaged_ratio[int((index-9)/int(binmerge))])

    ## Antiproton number
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/antiprotonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as pbarNum_3B:
            PbarNum_3B = pbarNum_3B.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/antiprotonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as pbarNum_3B:
            PbarNum_3B = pbarNum_3B.readlines()
    B3_PbarNumber = np.array(list(map(lambda s: float(s.strip()), PbarNum_3B)))

    ## Proton number
    if rigidityrange == "low":
        with open (lowworkpath + "/totalall" + "/results/protonnumber_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as protonNum_3B:
            ProtonNum_3B = protonNum_3B.readlines()
    elif rigidityrange == "intermediate":
        with open (intermediateworkpath + "/total" + "/results/protonnumber_" + trackerpattern + "_" + richcut + "_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+int(binmerge)]) + "_3BartalRotation.txt") as protonNum_3B:
            ProtonNum_3B = protonNum_3B.readlines()
    B3_ProtonNumber = np.array(list(map(lambda s: float(s.strip()), ProtonNum_3B)))


    #### Acceptance correction for time dependent result 
    '''
    Low:
    RigidityIndex is 1 1.16_1.51   accindex=(2+3)/2
    RigidityIndex is 3 1.51_1.92   accindex=(4+5)/2
    RigidityIndex is 5 1.92_2.4    accindex=(6+7)/2
    RigidityIndex is n             accindex=(n+1 + n+2)/2
    Intermediate:
    RigidityIndex is 9  2.97_3.64
    RigidityIndex is 11 3.64_4.43
    RigidityIndex is 13 
    RigidityIndex is 15
    RigidityIndex is n             accindex=(n+1 + n+2)/2
    '''
    if rigidityrange == "low":
        effective_acceptance_correction = effective_acceptance_correction_low 
    elif rigidityrange == "intermediate":
        effective_acceptance_correction = effective_acceptance_correction_intermediate

    M6result = M6result * (effective_acceptance_correction[index+1] + effective_acceptance_correction[index + 1 + (int(binmerge)-1)])/2
    B6result = B6result * (effective_acceptance_correction[index+1] + effective_acceptance_correction[index + 1 + (int(binmerge)-1)])/2 
    B3result = B3result * (effective_acceptance_correction[index+1] + effective_acceptance_correction[index + 1 + (int(binmerge)-1)])/2 


    #### Save date time and numpy for All rigidity
    if len(B3result_all) == 0:
        B3result_all         = B3result
        Error3b_all          = Error3b
        B3_relative_all      = B3_relative
        B3relative_error_all = B3relative_error
        #B3_PbarNumber_all    = B3_PbarNumber
        #B3_ProtonNumber_all  = B3_ProtonNumber

        B6result_all         = B6result
        Error6b_all          = Error6b
        B6_relative_all      = B6_relative
        B6relative_error_all = B6relative_error
        #B6_PbarNumber_all    = B6_PbarNumber
        #B6_ProtonNumber_all  = B6_ProtonNumber

        M6result_all         = M6result
        Error6m_all          = ErrorM6
        M6_relative_all      = M6_relative
        M6relative_error_all = M6relative_error
        #M6_PbarNumber_all    = M6_PbarNumber
        #M6_ProtonNumber_all  = M6_ProtonNumber

    else:
        B3result_all         = np.row_stack((B3result_all        , B3result))
        Error3b_all          = np.row_stack((Error3b_all         , Error3b))
        B3_relative_all      = np.row_stack((B3_relative_all     , B3_relative))
        B3relative_error_all = np.row_stack((B3relative_error_all, B3relative_error))
        #B3_PbarNumber_all    = np.row_stack((B3_PbarNumber_all   , B3_PbarNumber))
        #B3_ProtonNumber_all  = np.row_stack((B3_ProtonNumber_all , B3_ProtonNumber))

        B6result_all         = np.row_stack((B6result_all        , B6result))
        Error6b_all          = np.row_stack((Error6b_all         , Error6b))
        B6_relative_all      = np.row_stack((B6_relative_all     , B6_relative))
        B6relative_error_all = np.row_stack((B6relative_error_all, B6relative_error))
        #B6_PbarNumber_all    = np.row_stack((B6_PbarNumber_all   , B6_PbarNumber))
        #B6_ProtonNumber_all  = np.row_stack((B6_ProtonNumber_all , B6_ProtonNumber))

        M6result_all         = np.row_stack((M6result_all        , M6result))
        Error6m_all          = np.row_stack((Error6m_all         , ErrorM6))
        M6_relative_all      = np.row_stack((M6_relative_all     , M6_relative))
        M6relative_error_all = np.row_stack((M6relative_error_all,M6relative_error))
        #M6_PbarNumber_all    = np.row_stack((M6_PbarNumber_all   , M6_PbarNumber))
        #M6_ProtonNumber_all  = np.row_stack((M6_ProtonNumber_all , M6_ProtonNumber))

    B3_PbarNumber_all    = B3_PbarNumber
    B3_ProtonNumber_all  = B3_ProtonNumber
    B6_PbarNumber_all    = B6_PbarNumber
    B6_ProtonNumber_all  = B6_ProtonNumber
    M6_PbarNumber_all    = M6_PbarNumber
    M6_ProtonNumber_all  = M6_ProtonNumber




    #### Load Time dates
    dates_unixsecond = binning.Bartals3Unixtime[0:B3result.shape[0]]                                           # first bartel rotation included.
    dates            = [datetime.fromtimestamp(i) for i in dates_unixsecond]                                   # first bartel rotation included.
    dates_IhepRange       = [datetime.fromtimestamp(i) for i in binning.Bartals3Unixtime[0:Ihep_all[0].shape[0]]]   # first bartel rotation included.
    dates_TaiwanRange     = [datetime.fromtimestamp(i) for i in binning.Bartals3Unixtime[0:Taiwan_all[0].shape[0]]] # first bartel rotation included.

    dates_unixsecond_FirstRemoved = binning.Bartals3Unixtime_WithoutFirstBartel[0:B3result.shape[0]]   # first bartel rotation removed.
    dates_FirstRemoved            = [datetime.fromtimestamp(i) for i in dates_unixsecond_FirstRemoved] # first bartel rotation removed.
    dates_IhepRange_FirstRemoved       = [datetime.fromtimestamp(i) for i in binning.Bartals3Unixtime_WithoutFirstBartel[0:Ihep_all[0].shape[0]]]   # first bartel rotation removed.
    dates_TaiwanRange_FirstRemoved     = [datetime.fromtimestamp(i) for i in binning.Bartals3Unixtime_WithoutFirstBartel[0:Taiwan_all[0].shape[0]]] # first bartel rotation removed.

    dates_1Bartel_unixsecond       = binning.Bartals1Unixtime[0:PosiOverElec.shape[0]]
    dates_1Bartel                  = [datetime.fromtimestamp(i) for i in dates_1Bartel_unixsecond]
    dates_1Bartel_unixsecond_Aslam = binning.Bartals1Unixtime
    dates_1Bartel_Aslam            = [datetime.fromtimestamp(i) for i in dates_1Bartel_unixsecond_Aslam]

    datesM6_unixsecond = binning.Monthes6Unixtime[0:M6result.shape[0]]
    datesM6            = [datetime.fromtimestamp(i) for i in datesM6_unixsecond]

    dates6B             = [datetime.fromtimestamp(i) for i in binning.Bartals6Unixtime[0:B6result.shape[0]]]
    dates6B_TaiwanRange = [datetime.fromtimestamp(i) for i in binning.Bartals6Unixtime[0:Taiwan_all_6[0].shape[0]]]


    #### Load Systematic Error
    SystematicErrorAll_3B = []
    SystematicErrorAll_6B = []
    SystematicErrorAll_6M = []

    ## 6 months
    for timeindex in range(0, SplitTotal_6M, mergestep_6M):
        SystematicError_6M = []
        for rigidityindex in range(pointnumber):
            if rigidityrange == "low":
                rootfile = uproot.open(lowworkpath + "/totalall/SysError_Eff/SysErr_SigEff_" + modename_6M + '_' + str(timeindex) + ".root")
                SystematicError_6M.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])
            elif rigidityrange == "intermediate":
                rootfile = uproot.open(intermediateworkpath + "/total/SysError_Eff/SysErr_SigEff_" + modename_6M + '_' + str(timeindex) + ".root")
                SystematicError_6M.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])
        SystematicErrorAll_6M.append(SystematicError_6M)
    SystematicErrorAll_6M = np.array(SystematicErrorAll_6M)
    Error_Sys6M = SystematicErrorAll_6M[:,int((index-rigidity_start)/2)]/100000

    ## 3 Bartels
    for timeindex in range(0, SplitTotal_3B, mergestep_3B):
        SystematicError_3B = []   
        ## Remove 123-126 point
        if timeindex == 123:
            for rigidityindex in range(pointnumber):
                SystematicError_3B.append(0)
        else:
            for rigidityindex in range(pointnumber):
                if rigidityrange == "low":
                    rootfile = uproot.open(lowworkpath + "/totalall/SysError_Eff/SysErr_SigEff_" + modename_3B + '_' + str(timeindex) + ".root")
                    SystematicError_3B.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])
                elif rigidityrange == "intermediate":
                    rootfile = uproot.open(intermediateworkpath + "/total/SysError_Eff/SysErr_SigEff_" + modename_3B + '_' + str(timeindex) + ".root")
                    SystematicError_3B.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])
        SystematicErrorAll_3B.append(SystematicError_3B)
    SystematicErrorAll_3B = np.array(SystematicErrorAll_3B)
    Error_Sys3B = SystematicErrorAll_3B[:,int((index-rigidity_start)/2)]/100000

    ## 6 Bartels
    for timeindex in range(0, SplitTotal_6B, mergestep_6B):
        SystematicError_6B = []
        for rigidityindex in range(pointnumber):
            if rigidityrange == "low":
                rootfile = uproot.open(lowworkpath + "/totalall/SysError_Eff/SysErr_SigEff_" + modename_6B + '_' + str(timeindex) + ".root")
                SystematicError_6B.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])
            elif rigidityrange == "intermediate":
                rootfile = uproot.open(intermediateworkpath + "/total/SysError_Eff/SysErr_SigEff_" + modename_6B + '_' + str(timeindex) + ".root")
                SystematicError_6B.append(rootfile['g_SysError_SigEff'].values()[1][rigidityindex])
        SystematicErrorAll_6B.append(SystematicError_6B)
    SystematicErrorAll_6B = np.array(SystematicErrorAll_6B)
    Error_Sys6B = SystematicErrorAll_6B[:,int((index-rigidity_start)/2)]/100000


    #### Load Systematic Error due to Acceptance
    Error_SysAccAveraged      = ErrorAll_Acc[ int( (index - rigidity_start) / int(binmerge)) ]
    if rigidityrange == "low":
        Error_SysAcc6M = SystematicRelativeErrorAccpetance[int( (index - rigidity_start)/int(binmerge) )] * M6result
        Error_SysAcc6B = SystematicRelativeErrorAccpetance[int( (index - rigidity_start)/int(binmerge) )] * B6result
        Error_SysAcc3B = SystematicRelativeErrorAccpetance[int( (index - rigidity_start)/int(binmerge) )] * B3result
    elif rigidityrange == "intermediate":
        Error_SysAcc6M = SystematicRelativeErrorAccpetance[int( (index - rigidity_start)/int(binmerge) ) + 4] * M6result
        Error_SysAcc6B = SystematicRelativeErrorAccpetance[int( (index - rigidity_start)/int(binmerge) ) + 4] * B6result
        Error_SysAcc3B = SystematicRelativeErrorAccpetance[int( (index - rigidity_start)/int(binmerge) ) + 4] * B3result


    #### Calculate Total Error
    ## Time-Dependent Total Error 
    Error_TotalTimeDependent6M = np.sqrt( Error_Sys6M**2 + (ErrorM6)**2 )
    Error_TotalTimeDependent6B = np.sqrt( Error_Sys6B**2 + (Error6b)**2 )
    Error_TotalTimeDependent3B = np.sqrt( Error_Sys3B**2 + (Error3b)**2 )
    ## All Systematic Error and All Error
    Error_TotalSys6M = np.sqrt( Error_Sys6M**2 + Error_SysAcc6M**2 )
    Error_TotalSys6B = np.sqrt( Error_Sys6B**2 + Error_SysAcc6B**2 )
    Error_TotalSys3B = np.sqrt( Error_Sys3B**2 + Error_SysAcc3B**2 )
    Error_Total6M = np.sqrt( Error_Sys6M**2 + (ErrorM6)**2 + (Error_SysAcc6M)**2 )
    Error_Total6B = np.sqrt( Error_Sys6B**2 + (Error6b)**2 + (Error_SysAcc6B)**2 )
    Error_Total3B = np.sqrt( Error_Sys3B**2 + (Error3b)**2 + (Error_SysAcc3B)**2 )
    

    RelativeStatisticalError_6B = Error6b          / B6result
    RelativeSystematicError_6B  = Error_TotalSys6B / B6result
    RelativeStatisticalError_3B = Error3b          / B3result
    RelativeSystematicError_3B  = Error_TotalSys3B / B3result
    RelativeStatisticalError_6M = ErrorM6          / M6result
    RelativeSystematicError_6M  = Error_TotalSys6M / M6result



    if len(RelativeStatisticalError_3B_all) == 0:
        RelativeStatisticalError_3B_all = RelativeStatisticalError_3B
        RelativeSystematicError_3B_all  = RelativeSystematicError_3B
        RelativeStatisticalError_6B_all = RelativeStatisticalError_6B
        RelativeSystematicError_6B_all  = RelativeSystematicError_6B
        RelativeStatisticalError_6M_all = RelativeStatisticalError_6M
        RelativeSystematicError_6M_all  = RelativeSystematicError_6M
    else:
        RelativeStatisticalError_3B_all = np.row_stack((RelativeStatisticalError_3B_all, RelativeStatisticalError_3B))
        RelativeSystematicError_3B_all  = np.row_stack((RelativeSystematicError_3B_all , RelativeSystematicError_3B))
        RelativeStatisticalError_6B_all = np.row_stack((RelativeStatisticalError_6B_all, RelativeStatisticalError_6B))
        RelativeSystematicError_6B_all  = np.row_stack((RelativeSystematicError_6B_all , RelativeSystematicError_6B))
        RelativeStatisticalError_6M_all = np.row_stack((RelativeStatisticalError_6M_all, RelativeStatisticalError_6M))
        RelativeSystematicError_6M_all  = np.row_stack((RelativeSystematicError_6M_all , RelativeSystematicError_6M))












