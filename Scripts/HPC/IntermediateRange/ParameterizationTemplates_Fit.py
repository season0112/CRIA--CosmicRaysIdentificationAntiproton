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

#### OneGaussOneNovo

def CurveFit_OneGaussOneNovo(Dimension, RigidityRange, value, x_center, y_value, y_error, a_0, b_0, mean_gauss_0, sigma_gauss_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0, index):
    UsedFunction = ParameterizationTemplates_function.FitFunction
    InivialValue = [a_0, b_0, mean_gauss_0, sigma_gauss_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0]

    if RigidityRange == 'intermediate':
        if Dimension == "TRD":
            if value == "Antiproton":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Electron":
                #popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000, bounds=([0, 0, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]) )
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Pion":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)  ## bounds=([0, 0, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])

    elif RigidityRange == 'low':
        if Dimension == "TRD":
            if value == "Antiproton":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Electron":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000, sigma = y_error)
                #popt, pcov = curve_fit(UsedFunction, x_center[0:5], y_value[0:5], p0 = InivialValue, maxfev = 800000, sigma = y_error[0:5])
            elif value == "Pion":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)

    ParameterError = np.sqrt(np.diag( abs(pcov) ))  ## Check:  abs(pcov)

    return popt, pcov, ParameterError


def MinuitFit_OneGaussOneNovo(Dimension, RigidityRange, value, FitFunction, x_center, y_value, y_error, a_0, b_0, mean_gauss_0, sigma_gauss_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0, index):

    if RigidityRange == 'intermediate':
        if Dimension == "TRD":

            lsq_rough = LeastSquares(x_center, y_value, y_error, ParameterizationTemplates_function.FitFunction)
            #lsq_rough = LeastSquares(x_center[0:7], y_value[0:7], y_error[0:7], ParameterizationTemplates_function.FitFunction)
            m_rough = Minuit(lsq_rough, a=a_0, b=b_0, mean_gauss=mean_gauss_0, sigma_gauss=sigma_gauss_0, miu_Novo=miu_Novo_0, sigma_Novo=sigma_Novo_0, tau_Novo=tau_Novo_0)
            '''
            m_rough.limits['a'] = (0, np.inf)
            m_rough.limits['b'] = (0, np.inf)
            m_rough.limits['mean_gauss'] = (0, np.inf)
            m_rough.limits['sigma_gauss'] = (0, np.inf)
            m_rough.limits['miu_Novo'] = (0, np.inf)
            m_rough.limits['sigma_Novo'] = (0, np.inf)
            m_rough.limits['tau_Novo'] = (0, np.inf)
            '''
            m_rough.errordef = 1.0   #1.0. errordef should be 1.0 for a least-squares cost function and 0.5 for a negative log-likelihood function.
            #m_rough.errordef = 0.5
            m_rough.migrad()
            m_rough.hesse()
            #for p in m_rough.parameters:
            #    print('{} = {}'.format(p, uncertainties.ufloat(m_rough.values[p], m_rough.errors[p])))
            #print("Parameter: " + str(np.array(m_rough.values)))
            #print("Parameter Error: " + str(np.array(m_rough.errors)))
            #print(m_rough.covariance)
    elif RigidityRange == 'low':
        if Dimension == "TRD":
            if value == "Pion":
                lsq_rough = LeastSquares(x_center, y_value, y_error, ParameterizationTemplates_function.FitFunction)
                m_rough = Minuit(lsq_rough, a=a_0, b=b_0, mean_gauss=mean_gauss_0, sigma_gauss=sigma_gauss_0, miu_Novo=miu_Novo_0, sigma_Novo=sigma_Novo_0, tau_Novo=tau_Novo_0)
                m_rough.limits['a'] = (0, np.inf)
                m_rough.limits['b'] = (0, np.inf)
                m_rough.limits['mean_gauss'] = (0, np.inf)
                m_rough.limits['sigma_gauss'] = (0, np.inf)
                m_rough.limits['miu_Novo'] = (0, np.inf)
                m_rough.limits['sigma_Novo'] = (0, np.inf)
                m_rough.limits['tau_Novo'] = (0, np.inf)
                m_rough.errordef = 1.0
                m_rough.migrad()
                m_rough.hesse()


    return np.array(m_rough.values), np.array(m_rough.covariance), np.array(m_rough.errors)


#### OneGauss

def CurveFit_OneGauss(Dimension, RigidityRange, value, x_center, y_value, y_error, a_0, mean_gauss_0, sigma_gauss_0, index):

    UsedFunction = ParameterizationTemplates_function.Gauss_function
    InivialValue = [a_0, mean_gauss_0, sigma_gauss_0]

    if RigidityRange == 'intermediate':
        if Dimension == "TRD":
            if value == "Antiproton":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Electron":
                #popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000, bounds=([0, 0, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]) )
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Pion":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)  ## bounds=([0, 0, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])

    elif RigidityRange == 'low':
        if Dimension == "TRD":
            if value == "Antiproton":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Electron":
                #popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000, sigma = y_error)
                popt, pcov = curve_fit(UsedFunction, x_center[0:11], y_value[0:11], p0 = InivialValue, maxfev = 800000, sigma = y_error[0:11])
            elif value == "Pion":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
        elif Dimension == "TOF":
            if value == "Antiproton":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Electron":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Pion":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)

    ParameterError = np.sqrt(np.diag( abs(pcov) ))  ## Check:  abs(pcov)

    return popt, pcov, ParameterError


def MinuitFit_OneGauss(Dimension, RigidityRange, value, FitFunction, x_center, y_value, y_error, a_0, mean_gauss_0, sigma_gauss_0):

    if RigidityRange == 'intermediate':
        print("in progress")   
    elif RigidityRange == 'low':
        if Dimension == "TRD":
            if value == "Antiproton":
                lsq_rough = LeastSquares(x_center, y_value, y_error, ParameterizationTemplates_function.Gauss_function)
        elif Dimension == "TOF":
            if value == "Antiproton":
                lsq_rough = LeastSquares(x_center, y_value, y_error, ParameterizationTemplates_function.Gauss_function)
            elif value == "Electron":
                lsq_rough = LeastSquares(x_center, y_value, y_error, ParameterizationTemplates_function.Gauss_function)
            elif value == "Pion":
                lsq_rough = LeastSquares(x_center, y_value, y_error, ParameterizationTemplates_function.Gauss_function)

    m_rough = Minuit(lsq_rough, a=a_0, mean_gauss=mean_gauss_0, sigma_gauss=sigma_gauss_0)
    m_rough.errordef = 1.0
    m_rough.migrad()
    m_rough.hesse()
    #print("Parameter: " + str(np.array(m_rough.values)))
    #print("Parameter Error: " + str(np.array(m_rough.errors)))
    #print(m_rough.covariance)

    return np.array(m_rough.values), np.array(m_rough.covariance), np.array(m_rough.errors)


#### ExponentialFun

def CurveFit_ExponentialFun(Dimension, RigidityRange, value, FitFunction, x_center, y_value, y_error, k1, k2, b, c, index):

    UsedFunction = ParameterizationTemplates_function.ExponentialFun
    InivialValue = [k1, k2, b, c]

    if RigidityRange == 'intermediate':
        if Dimension == "TRD":
            print("in progress")

    elif RigidityRange == 'low':
        if Dimension == "TRD":
            if value == "Electron":
                #popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000, sigma = y_errori) 
                #popt, pcov = curve_fit(UsedFunction, x_center[0:4], y_value[0:4], p0 = InivialValue, maxfev = 1000, sigma = y_error[0:4]) ## bounds=([0, 0, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])
                if index == 2:
                    popt, pcov = curve_fit(UsedFunction, x_center[2:], y_value[2:], p0 = InivialValue, maxfev = 800000)
                else:
                    popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)

    ParameterError = np.sqrt(np.diag( abs(pcov) ))  ## Check:  abs(pcov)

    return popt, pcov, ParameterError


def MinuitFit_ExponentialFun(Dimension, RigidityRange, value, FitFunction, x_center, y_value, y_error, k1_0, k2_0, b_0, c_0, index):

    if RigidityRange == 'intermediate':
        if Dimension == "TRD":
            print("in progress")

    elif RigidityRange == 'low':
        if Dimension == "TRD":
            if value == "Electron":
                if index == 0:
                    startindex = 0
                    endindex   = 10
                elif index == 1:
                    startindex = 0
                    endindex   = 8
                elif index == 2:
                    startindex = 4
                    endindex   = 8
                elif index == 7:
                    startindex = 0
                    endindex   = 10
                elif index == 8:
                    startindex = 0
                    endindex   = 10
                elif index == 10:
                    startindex = 0
                    endindex   = 10
                elif index == 11:
                    startindex = 0
                    endindex   = 10
                elif index == 12:
                    startindex = 0
                    endindex   = 10
                elif index == 13:
                    startindex = 1
                    endindex   = 18
                elif index == 14:
                    startindex = 0
                    endindex   = 10
                elif index == 15:
                    startindex = 0
                    endindex   = 10

                else:
                    startindex = 0
                    endindex   = x_center.shape[0]

                lsq_rough = LeastSquares(x_center[startindex:endindex], y_value[startindex:endindex], y_error[startindex:endindex], ParameterizationTemplates_function.ExponentialFun)

                m_rough = Minuit(lsq_rough, k1=k1_0, k2=k2_0, b=b_0, c=c_0)
                m_rough.errordef = 1.0
                m_rough.migrad()
                m_rough.hesse()

    return np.array(m_rough.values), np.array(m_rough.covariance), np.array(m_rough.errors)


#### OneNovo

def CurveFit_OneNovo(Dimension, RigidityRange, value, x_center, y_value, y_error, b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0, index):

    UsedFunction = ParameterizationTemplates_function.Novosibirsk_function
    InivialValue = [b_0, miu_Novo_0, sigma_Novo_0, tau_Novo_0]

    if RigidityRange == 'intermediate':
        if Dimension == "TRD":
            if value == "Antiproton":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Electron":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)

    elif RigidityRange == 'low':
        if Dimension == "TRD":
            if value == "Pion":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
        elif Dimension == "TOF":
            if value == "Electron":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Pion":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)
            elif value == "Antiproton":
                popt, pcov = curve_fit(UsedFunction, x_center, y_value, p0 = InivialValue, maxfev = 800000)


    ParameterError = np.sqrt(np.diag( abs(pcov) ))  ## Check:  abs(pcov)

    return popt, pcov, ParameterError















