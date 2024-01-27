import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F
import uproot
import scipy as sp
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import binning
import argparse
import math
from iminuit import Minuit
from iminuit.util import describe, make_func_code
from iminuit.cost import LeastSquares
import uncertainties
import ParameterizationTemplates_function
import inspect, re


#### TRD

def TRD_InitialValues_OneGaussOneNovo(TemplateName, RigidityIndex, RigidityRange):

    if RigidityRange == 'intermediate':
        if TemplateName == "Antiproton":
            if RigidityIndex == 7:
                a_0 = 20000
                b_0 = 800000
                mean_gauss_0   = -1.0
                sigma_gauss_0  = 0.5
                miu_Novo_0     = -1.0    # (peak position)
                sigma_Novo_0   = 0.2     # (width of the peak)
                tau_Novo_0     = 0.01    # (parameter describing the tail of the distribution)
            else:
                a_0 = 500000
                b_0 = 500000
                mean_gauss_0   = -1.0
                sigma_gauss_0  = 0.5
                miu_Novo_0     = -1.0    # (peak position)
                sigma_Novo_0   = 0.2     # (width of the peak)
                tau_Novo_0     = 0.01    # (parameter describing the tail of the distribution)

        elif TemplateName == "Electron":

            '''
            a_0 = 10000
            b_0 = 10000
            mean_gauss_0   = -0.35
            sigma_gauss_0  = 0.1
            miu_Novo_0     = -0.35    # (peak position)
            sigma_Novo_0   = 0.2     # (width of the peak)
            tau_Novo_0     = 0.1    # (parameter describing the tail of the distribution)
            '''
            
            if RigidityIndex < 17:
                a_0 = 300
                b_0 = 7000
                mean_gauss_0   = -0.3
                sigma_gauss_0  = 0.1
                miu_Novo_0     = -0.35    # (peak position)
                sigma_Novo_0   = 0.2     # (width of the peak)
                tau_Novo_0     = 0.1    # (parameter describing the tail of the distribution)
            else:
                a_0 = 284
                b_0 = 13911
                mean_gauss_0   = -0.26
                sigma_gauss_0  = 0.05
                miu_Novo_0     = -0.36    # (peak position)
                sigma_Novo_0   = 0.09     # (width of the peak)
                tau_Novo_0     = -0.22    # (parameter describing the tail of the distribution)
            

            '''
            if RigidityIndex == 0:
                a_0 = 2000
                b_0 = 22000
                mean_gauss_0   = -0.4
                sigma_gauss_0  = 0.1
                miu_Novo_0     = -0.3    # (peak position)
                sigma_Novo_0   = 0.2     # (width of the peak)
                tau_Novo_0     = 0.1    # (parameter describing the tail of the distribution)
            elif RigidityIndex == 2:
                a_0 = 315
                b_0 = 26026
                mean_gauss_0   = -0.4
                sigma_gauss_0  = 0.09
                miu_Novo_0     = -0.3    # (peak position)
                sigma_Novo_0   = -0.09     # (width of the peak)
                tau_Novo_0     = 0.3    # (parameter describing the tail of the distribution)
            elif RigidityIndex == 1 or 2<RigidityIndex<=12:
                a_0 = 2000
                b_0 = 22000
                mean_gauss_0   = -0.4
                sigma_gauss_0  = 0.1
                miu_Novo_0     = -0.3    # (peak position)
                sigma_Novo_0   = -0.09     # (width of the peak)
                tau_Novo_0     = 0.3    # (parameter describing the tail of the distribution)
            elif 13<=RigidityIndex<=18:
                a_0 = 1000
                b_0 = 9000 #(14000)
                mean_gauss_0   = -0.4
                sigma_gauss_0  = 0.1
                miu_Novo_0     = -0.3    # (peak position)
                sigma_Novo_0   = -0.09     # (width of the peak)
                tau_Novo_0     = 0.3    # (parameter describing the tail of the distribution)
            elif RigidityIndex == 19:
                a_0 = 300
                b_0 = 7000
                mean_gauss_0   = -0.3
                sigma_gauss_0  = 0.1
                miu_Novo_0     = -0.35    # (peak position)
                sigma_Novo_0   = 0.2     # (width of the peak)
                tau_Novo_0     = 0.1    # (parameter describing the tail of the distribution)
            '''
        elif TemplateName == "Pion":
            a_0 = 500000
            b_0 = 500000
            mean_gauss_0   = -0.4
            sigma_gauss_0  = 0.5
            miu_Novo_0     = -1.0    # (peak position)
            sigma_Novo_0   = 0.2     # (width of the peak)
            tau_Novo_0     = 0.01    # (parameter describing the tail of the distribution)

    elif RigidityRange == 'low':
        if TemplateName == "Antiproton":
            if RigidityIndex == 2:
                a_0 = 200000
                b_0 = 5000
                mean_gauss_0   = -0.9
                sigma_gauss_0  = 0.1
                miu_Novo_0     = -0.9    # (peak position)
                sigma_Novo_0   = 0.2     # (width of the peak)
                tau_Novo_0     = 0.01    # (parameter describing the tail of the distribution)
            else:
                a_0 = 2000
                b_0 = 5437499
                mean_gauss_0   = -1.0
                sigma_gauss_0  = 0.5
                miu_Novo_0     = -1.0    # (peak position)
                sigma_Novo_0   = 0.2     # (width of the peak)
                tau_Novo_0     = 0.01    # (parameter describing the tail of the distribution)
        elif TemplateName == "Electron":
            a_0 = 7000000
            b_0 = 5437499
            mean_gauss_0   = -0.44
            sigma_gauss_0  = 0.3
            miu_Novo_0     = -1.0    # (peak position)
            sigma_Novo_0   = 0.2     # (width of the peak)
            tau_Novo_0     = 0.01    # (parameter describing the tail of the distribution)
        elif TemplateName == "Pion":
            if RigidityIndex < 2:
                a_0 = 5000
                b_0 = 5000
                mean_gauss_0   = 0.85
                sigma_gauss_0  = 0.1
                miu_Novo_0     = 0.85    # (peak position)
                sigma_Novo_0   = 0.1     # (width of the peak)
                tau_Novo_0     = 0.01    # (parameter describing the tail of the distribution)
            elif RigidityIndex == 2:
                a_0 = 2081
                b_0 = 2000
                mean_gauss_0   = 0.94
                sigma_gauss_0  = 0.06
                miu_Novo_0     = 0.94    # (peak position)
                sigma_Novo_0   = 0.05     # (width of the peak)
                tau_Novo_0     = 6.90004262e-05    # (parameter describing the tail of the distribution)
            elif RigidityIndex > 3 and RigidityIndex < 11:
                a_0 = 2081
                b_0 = 15463
                mean_gauss_0   = 0.94
                sigma_gauss_0  = 0.12
                miu_Novo_0     = 0.036    # (peak position)
                sigma_Novo_0   = 0.4     # (width of the peak)
                tau_Novo_0     = 6.90004262e-05    # (parameter describing the tail of the distribution)
            elif RigidityIndex == 10:
                #Para_TRD:[1.50571329e+03 6.68920365e+03 9.40554689e-01 1.58781537e-01
                # 4.33409145e-02 4.87825005e-01 1.75288267e-04]
                a_0 = 1.50571329e+03
                b_0 = 6.68920365e+03
                mean_gauss_0   = 9.40554689e-01
                sigma_gauss_0  = 1.58781537e-01
                miu_Novo_0     = 4.33409145e-02    # (peak position)
                sigma_Novo_0   = 4.87825005e-01     # (width of the peak)
                tau_Novo_0     = 1.75288267e-04    # (parameter describing the tail of the distribution)
            else:
                a_0 = 998
                b_0 = 2846
                mean_gauss_0   = 0.89
                sigma_gauss_0  = 0.17
                miu_Novo_0     = 0.02    # (peak position)
                sigma_Novo_0   = 0.49     # (width of the peak)
                tau_Novo_0     = 6.76153257e-04    # (parameter describing the tail of the distribution)

    return a_0, b_0, mean_gauss_0, sigma_gauss_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0


def TRD_InitialValues_OneGauss(TemplateName, RigidityIndex, RigidityRange):

    if RigidityRange == 'low':
        if TemplateName == "Antiproton":
            if RigidityIndex == 2:
                a_0 = 200000
                mean_gauss_0   = 0.9
                sigma_gauss_0  = 0.1
            else:
                a_0 = 1.03455290e+05
                mean_gauss_0   = 9.33482797e-01
                sigma_gauss_0  = 1.34921533e-01
        elif TemplateName == "Electron":
            a_0 = 7000000
            mean_gauss_0   = -0.44
            sigma_gauss_0  = 0.3
        elif TemplateName == "Pion":
            a_0 = 2000
            mean_gauss_0   = -1.0
            sigma_gauss_0  = 0.5

    return a_0, mean_gauss_0, sigma_gauss_0


def TRD_InitialValues_OneNovo(TemplateName, RigidityIndex, RigidityRange):

    if RigidityRange == 'low':
        if TemplateName == "Pion":
            b_0          = 12000
            miu_Novo_0   = 0.85
            sigma_Novo_0 = 0.1
            tau_Novo_0   = 0.01

    return b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0


def TRD_InitialValues_ExponentialFun(TemplateName, RigidityIndex, RigidityRange):
    #Para_TRD:[  3.94041455 -11.06070953  14.16457486   2.36844181]
    #Para_TRD:[  9.09120702 -11.37522144  13.96094275   2.03454667]
    if RigidityRange == 'low':
        if TemplateName == "Electron":
            k1 = 4
            k2 = -11
            b  = 14
            c  = 2
    return k1, k2, b, c


#### TOF
def TOF_InitialValues_OneGauss(TemplateName, RigidityIndex, RigidityRange):

    if RigidityRange == 'low':
        if TemplateName == "Antiproton":
            a_0 = 80000
            mean_gauss_0   = 0.0
            sigma_gauss_0  = 0.1
        elif TemplateName == "Electron":
            a_0 = 2000
            mean_gauss_0   = -0.35
            sigma_gauss_0  = 0.2

        '''
            elif TemplateName == "Pion":
            a_0 = 2000
            mean_gauss_0   = -1.0
            sigma_gauss_0  = 0.5
        '''

    return a_0, mean_gauss_0, sigma_gauss_0


def TOF_InitialValues_OneNovo(TemplateName, RigidityIndex, RigidityRange):

    if RigidityRange == 'low':
        if TemplateName == "Electron":
            if RigidityIndex < 2:
                b_0          = 2000
                miu_Novo_0   = -0.32
                sigma_Novo_0 = 0.1
                tau_Novo_0   = 0.01
            else:
                b_0          = 300
                miu_Novo_0   = -0.1
                sigma_Novo_0 = 0.1
                tau_Novo_0   = 0.01

        elif TemplateName == "Pion":
            if RigidityIndex < 2:
                b_0          = 10000
                miu_Novo_0   = -0.3
                sigma_Novo_0 = 0.1
                tau_Novo_0   = 0.01
            else:
                b_0          = 7000
                miu_Novo_0   = -0.1
                sigma_Novo_0 = 0.1
                tau_Novo_0   = 0.01

    return b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0


#### additional tool

def FixForNaf(dic):

    for key, value in dic.items():
        if key == "variable":
            for EachValue in value:
                for i in range(len(EachValue)):
                    if EachValue[i] == 0 or math.isnan(EachValue[i]) == True:
                        EachValue[i] = 1
        elif key == "error":
            for EachValue in value:
                for i in range(len(EachValue)):
                    if EachValue[i] == 0 or math.isnan(EachValue[i]) == True:
                        EachValue[i] = 1

