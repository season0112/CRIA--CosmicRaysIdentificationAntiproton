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
import ROOT
import uncertainties
from uncertainties import unumpy
import scipy as sp
import scipy.stats as stats
from pylab import *
from kapteyn import kmpfit
import PythonPlotDefaultParameters


#### Load
def LoadResult(workpath, pattern, NNsuffix, ISSversion, Binningversion, SignalEfficiency, CCcut_TF, UsedRatioName):
    #### Two options: (1). Use CCProtonOverProtonRatio_ISS/MC fit Uncertainty band as error; (2). ISSCCLevel/RigidityResolution fit Uncertainty Band as error

    if UsedRatioName == "CCLevelISSToMc":
        f_CCLevel = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")

        CCLevel_ISS      = np.asarray(f_CCLevel.Get("CCLevelISS_Signal_Efficiency_"       + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        CCLevelError_ISS = np.asarray(f_CCLevel.Get("CCLevelISS_Error_Signal_Efficiency_" + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        CCLevel_MC       = np.asarray(f_CCLevel.Get("CCLevelMC_Signal_Efficiency_"        + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        CCLevelError_MC  = np.asarray(f_CCLevel.Get("CCLevellMC_Error_Signal_Efficiency_" + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        CCLevel_ISSOverMC       = np.asarray(f_CCLevel.Get("CLevelRatio_Signal_Efficiency_"      + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        CCLevel_ISSOverMC_Error = np.asarray(f_CCLevel.Get("CLevelRatioError_Signal_Efficiency_" + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        RigidityBinningCenter   = np.asarray(f_CCLevel.Get("RigidityBinningCenter"))
        f_CCLevel.Close()

        ValueA          = CCLevel_ISS
        ValueAError     = CCLevelError_ISS
        ValueB          = CCLevel_MC
        ValueBError     = CCLevelError_MC
        RatioUserd      = CCLevel_ISSOverMC
        RatioErrorUserd = CCLevel_ISSOverMC_Error

    elif UsedRatioName == "CCProtonOverProtonRatio":
        f_CCProtonOverProton_ISSOverMC = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")

        CCProtonOverProtonRatio_ISS             = np.asarray(f_CCProtonOverProton_ISSOverMC.Get("CCProtonOverProtonISS_Signal_Efficiency_"       + str(SignalEfficiency) + "_CCcut_TF_"+ str(CCcut_TF)))
        CCProtonOverProtonRatioError_ISS        = np.asarray(f_CCProtonOverProton_ISSOverMC.Get("CCProtonOverProtonISS_Error_Signal_Efficiency_" + str(SignalEfficiency) + "_CCcut_TF_"+ str(CCcut_TF)))
        CCProtonOverProtonRatio_MC              = np.asarray(f_CCProtonOverProton_ISSOverMC.Get("CCProtonOverProtonMC_Signal_Efficiency_"        + str(SignalEfficiency) + "_CCcut_TF_"+ str(CCcut_TF)))
        CCProtonOverProtonRatioError_MC         = np.asarray(f_CCProtonOverProton_ISSOverMC.Get("CCProtonOverProtonMC_Error_Signal_Efficiency_"  + str(SignalEfficiency) + "_CCcut_TF_"+ str(CCcut_TF)))
        CCProtonOverProtonRatio_ISSOverMC       = np.asarray(f_CCProtonOverProton_ISSOverMC.Get("CCProtonOverProtonRatio_ISSOverMC_Signal_Efficiency_"       + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        CCProtonOverProtonRatio_ISSOverMC_Error = np.asarray(f_CCProtonOverProton_ISSOverMC.Get("CCProtonOverProtonRatio_ISSOverMC_Error_Signal_Efficiency_" + str(SignalEfficiency) + "_CCcut_TF_" + str(CCcut_TF)))
        RigidityBinningCenter                   = np.asarray(f_CCProtonOverProton_ISSOverMC.Get("RigidityBinningCenter"))
        f_CCProtonOverProton_ISSOverMC.Close()

        ValueA          = CCProtonOverProtonRatio_ISS
        ValueAError     = CCProtonOverProtonRatioError_ISS
        ValueB          = CCProtonOverProtonRatio_MC
        ValueBError     = CCProtonOverProtonRatioError_MC
        RatioUserd      = CCProtonOverProtonRatio_ISSOverMC
        RatioErrorUserd = CCProtonOverProtonRatio_ISSOverMC_Error 

    elif UsedRatioName == "CCLevelOverRigidityResolution":
        f_CCLevelISSOverResolution = TFile.Open(workpath + "/templatefit/negative/FixedCC/ISSCCLevelOverRigidityResolution_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")
        ISSCCLevelOverRigidityResolution      = np.asarray(f_CCLevelISSOverResolution.Get("CCLevelISSOverResolution"))
        ISSCCLevelOverRigidityResolutionError = np.asarray(f_CCLevelISSOverResolution.Get("CCLevelISSOverResolution_Error"))
        RigidityBinningCenter                 = np.asarray(f_CCLevelISSOverResolution.Get("RigidityBinningCenter"))
        ## Rescale to mean. FIXME
        ISSCCLevelOverRigidityResolution = ISSCCLevelOverRigidityResolution / mean(ISSCCLevelOverRigidityResolution)
        f_CCLevelISSOverResolution.Close()
        RatioUserd      = ISSCCLevelOverRigidityResolution
        RatioErrorUserd = ISSCCLevelOverRigidityResolutionError

    return ValueA, ValueAError, ValueB, ValueBError, RatioUserd, RatioErrorUserd, RigidityBinningCenter


#### Manually Calculate fit and uncertainty band

def Calculatet_CI(t, s_err, n, x, x2):   #Calculatet_CI(t, s_err, n, NNpoint, x_more)
    # ////
    # The Procesure is described in: https://www.coder.work/article/6527836
    # References
    #    [1] M. Duarte.  "Curve fitting," Jupyter Notebook.
    #    http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb
    ci = t * s_err * np.sqrt(1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
    return ci

def GeneralPolynomial(a, x):
    """Return a 1D polynomial."""
    return np.polyval(a, x)

def LinearFit(NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, RigidityBinningCenter):

    ## Polyfit
    #p, cov = np.polyfit(RigidityBinningCenter, RatioUserd, 1, cov=True)
    p, cov = np.polyfit(RigidityBinningCenter[PointsForFit_Begin:], RatioUserd[PointsForFit_Begin:], 0, w=1/RatioErrorUserd, cov='unscaled')  ## Const Fit

    y_model = GeneralPolynomial(p, RigidityBinningCenter[PointsForFit_Begin:])
    print( "Fit parameters: " + str(p))
    print( "Fir parameters error: " + str(np.sqrt(np.diag(cov))))
    #print("y_model:" + str(y_model))

    ## Statistics
    # Knowledge: pdf: Probability Density Function, cdf:Cumulative Distribution Function, ppf:Percent Point Function (Inverse of cdf)
    # Knowledge: ppf: Percent Point Function (or Inverse Cumulative Distribution Function), ppf returns the value x of the variable that has a given cumulative distribution probability (cdf). Thus, given the cdf(x) of a x value, ppf returns the value x itself, therefore, operating as the inverse of cdf.
    #n = NNpoint[PointsForFit_Begin:PointsForFit_End].size
    n = RigidityBinningCenter[PointsForFit_Begin:].size
    m = p.size
    dof = n - m

    # Knowledge: stats.t: A Student’s t continuous random variable.  ppf(q, df, loc=0, scale=1): q:Percent point, df=Degree Of Freedom
    t_95 = stats.t.ppf(0.975, dof)  # Because 95% confidence level is 0.975 Percent point
    #t_67 = stats.t.ppf(0.975, dof)    # Because 67% confidence level is 0.83 Percent point
    t_67 = stats.t.ppf(0.83, dof)    # Because 67% confidence level is 0.83 Percent point

    # Estimates of Error in Data/Model
    resid = RatioUserd[PointsForFit_Begin:] - y_model
    #print('Residual**2:' + str(resid**2))
    #print('Sigma**2:' + str(RatioErrorUserd[PointsForFit_Begin:]**2))
    chi2 = np.sum(resid**2 / RatioErrorUserd[PointsForFit_Begin:]**2) 
    print('chi2:' + str(chi2))
    print('dof:' + str(dof))
    chi2_red = chi2 / dof                                      # reduced chi-squared; measures goodness of fit
    print("Chi2/dof: " + str(chi2_red))
    s_err = np.sqrt(np.sum(resid**2) / dof)                    # standard deviation of the error

    # Here use fit function to extrapolate into all Rigidity range.
    y_more = GeneralPolynomial(p, x_more)
    #y_NNpoint = GeneralPolynomial(p, NNpoint)
    y_NNpoint = GeneralPolynomial(p, RigidityBinningCenter)

    # Calcualte Confidence Interval
    #ConfidenceInterval_show_67 = Calculatet_CI(t_67, s_err, n, NNpoint, x_more)
    #ConfidenceInterval_show_95 = Calculatet_CI(t_95, s_err, n, NNpoint, x_more)
    #ConfidenceInterval         = Calculatet_CI(t_67, s_err, n, NNpoint, NNpoint)
    ConfidenceInterval_show_67 = Calculatet_CI(t_67, s_err, n, RigidityBinningCenter, x_more)
    ConfidenceInterval_show_95 = Calculatet_CI(t_95, s_err, n, RigidityBinningCenter, x_more)
    ConfidenceInterval         = Calculatet_CI(t_67, s_err, n, RigidityBinningCenter, NNpoint)
    #print("LinearFit Confidence Band:")
    #print(str(ConfidenceInterval))

    # Calculate Prediction Interval
    PredictionInterval_show_67 = t_67 * s_err * np.sqrt(1 + 1/n + (x_more - np.mean(RigidityBinningCenter))**2 / np.sum((RigidityBinningCenter - np.mean(RigidityBinningCenter))**2))  
    PredictionInterval_show_95 = t_95 * s_err * np.sqrt(1 + 1/n + (x_more - np.mean(RigidityBinningCenter))**2 / np.sum((RigidityBinningCenter - np.mean(RigidityBinningCenter))**2))

    return y_model, y_more, ConfidenceInterval_show_67, ConfidenceInterval, PredictionInterval_show_67, PredictionInterval_show_95


#### KMP Fit 
def PolyFit0Model_ForCurveFit(x, p):
    return p

def PolyFit2Model_ForCurveFit(x, p1, p2, p3):
    return p1*x**2 + p2*x**1 + p3

def PolyFit0Model(p, x):
    a1 = p
    return a1

def LinearFitModel(p, x):
    a1, a2 = p
    return a1*x**1 + a2

def PolyFit3Model(p, x):
    a1, a2, a3, a4= p
    return a1*x**3 + a2*x**2 + a3*x**1 + a4

def PolyFit4Model(p, x):
    a1, a2, a3, a4, a5 = p
    return a1*x**4 + a2*x**3 + a3*x**2 + a4*x**1 + a5

def dfdpPolyFit0(x):
    dfdp = [1]
    return dfdp

def dfdpLinearFit(x):
    dfdp = [x, 1]
    return dfdp

def dfdpPolyFit3(x):
    dfdp = [x**3, x**2, x, 1]
    return dfdp

def dfdpPolyFit4(x):
    dfdp = [x**4, x**3, x**2, x, 1]
    return dfdp

def KMPFit(FitModel, InitialValue, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, dfdpfunction, RigidityBinningCenter):
    #KMPfit = kmpfit.simplefit(FitModel, InitialValue, NNpoint[PointsForFit_Begin:PointsForFit_End], RatioUserd[PointsForFit_Begin:PointsForFit_End])
    KMPfit = kmpfit.simplefit(FitModel, InitialValue, RigidityBinningCenter, RatioUserd, err=RatioErrorUserd)
    print("KMPFit Result Parameters (" + str(FitModel.__name__) + "): " + str(KMPfit.params))
    print("KMPFit Result Parameter Error:" + str(KMPfit.stderr))

    kmp_x = x_more
    dfdp = dfdpfunction(kmp_x)

    #kmp_x_save = NNpoint[0:32] ## Why 0:32 ???
    kmp_x_save = RigidityBinningCenter
    dfdp_save = dfdpfunction(kmp_x_save)

    # Fit points for illustration
    yhat_67     , upper_67     , lower_67      = KMPfit.confidence_band( x_more, dfdp, 0.67, FitModel, True)
    yhat_95     , upper_95     , lower_95      = KMPfit.confidence_band( x_more, dfdp, 0.95, FitModel, True)
    kmpCL_67      = upper_67 - yhat_67
    # Fit points for saving
    #yhat_67_save, upper_67_save, lower_67_save = KMPfit.confidence_band( NNpoint[0:32], dfdp_save, 0.67, FitModel, True)
    #yhat_95_save, upper_95_save, lower_95_save = KMPfit.confidence_band( NNpoint[0:32], dfdp_save, 0.95, FitModel, True)
    yhat_67_save, upper_67_save, lower_67_save = KMPfit.confidence_band( RigidityBinningCenter, dfdp_save, 0.67, FitModel, True)
    yhat_95_save, upper_95_save, lower_95_save = KMPfit.confidence_band( RigidityBinningCenter, dfdp_save, 0.95, FitModel, True)
    kmpCL_67_save = upper_67_save - yhat_67_save

    # Print KMP fit result
    print('KMP Fit Chi2/dof:' + str(KMPfit.rchi2_min))
    #print("KMP Fit Confidence Band:")
    #print(kmpCL_67)
    print("KMP Fit Confidence Band (Save):")
    print(kmpCL_67_save)

    return kmpCL_67, kmpCL_67_save, yhat_67, upper_67, lower_67


#### Plot

def Plot_ConfidenceInterval_Option1_LinearFit(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, y_model, x_more, y_more, ConfidenceInterval_show_67, PredictionInterval_show_67, PredictionInterval_show_95, RigidityBinningCenter, Name, Scaler, CCcut_TF):

    fig, ax = plt.subplots()

    plt.xlim(10,550)
    plt.ylim(0, 3)

    ## plot raw data
    #plt.errorbar( NNpoint[PointsForFit_Begin:PointsForFit_End], RatioUserd[PointsForFit_Begin:PointsForFit_End], yerr=RatioErrorUserd[PointsForFit_Begin:PointsForFit_End], xerr=0, fmt='o' )
    plt.errorbar( RigidityBinningCenter, RatioUserd, yerr=RatioErrorUserd/Scaler, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )

    ## plot Linear Fit (Option 1)
    plt.hlines(y_model[0], x_more[0], x_more[-1], linewidth=5, alpha=0.5, linestyle='dashed', color="red") # label="Linear Fit option1")
    ax.fill_between(x_more, y_more + ConfidenceInterval_show_67, y_more - ConfidenceInterval_show_67, color="yellow", edgecolor="", alpha=0.9, label='68% Confidence Interval') # label='67% Confidence Interval', color="#b9cfe7"
    #ax.fill_between(x_more, y_more + ConfidenceInterval_show_95, y_more - ConfidenceInterval_show_95, color="#e7b9cf", edgecolor="", alpha=0.5, label='LinearFit: 95% CL')
    #ax.plot(x_more, y_more - PredictionInterval_show_67, "--", color="0.5", label="68% Prediction Limits") #label="67% Prediction Limits"
    #ax.plot(x_more, y_more + PredictionInterval_show_67, "--", color="0.5")

    plt.xlabel('|R| / (GV)', horizontalalignment='right', x=1.0)
    if Name == 'CCLevel':
        plt.ylabel('CCLevel ratio (ISS/MC)'     , horizontalalignment='right', y=1.0)
    elif Name == 'CCProtonOverProton':
        plt.ylabel('CCProtonOverProton (ISS/MC)', horizontalalignment='right', y=1.0)

    plt.legend(fontsize=50)

    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_LinearFit_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LinearX.pdf')
    #plt.xscale('log')
    #plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_LinearFit_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LogX.pdf')

    plt.close()


def Plot_ValueAandB(ValueA, ValueAError, ValueB, ValueBError, RigidityBinningCenter, pattern, NNsuffix, Binningversion, suffix, workpath, Name, CCcut_TF):

    fig, ax = plt.subplots()

    plt.xlim(10,550)
    #plt.ylim(0.00000005, 0.001)

    plt.xlabel('|R| / (GV)' , horizontalalignment='right', x=1.0)
    plt.ylabel(Name         , horizontalalignment='right', y=1.0)

    ## plot raw data
    plt.errorbar( RigidityBinningCenter, ValueA, yerr=ValueAError, xerr=0, fmt='o', markerfacecolor="blue", ecolor="blue", linewidth=5, label='ISS')
    plt.errorbar( RigidityBinningCenter, ValueB, yerr=ValueBError, xerr=0, fmt='o', markerfacecolor="red", ecolor="red"  , linewidth=5, label='MC' )

    plt.legend(loc='best', fontsize=50)
    # get handles
    handles, labels = ax.get_legend_handles_labels()
    # remove the errorbars
    handles = [h[0] for h in handles]
    # use them in the legend
    ax.legend(handles, labels, loc='lower right', fontsize=70)

    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/' + Name + '_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + '_' + Binningversion + suffix + '.pdf')
    plt.yscale('log')
    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/' + Name + '_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + '_' + Binningversion + suffix + '_LogY.pdf')

    plt.close()


def Plot_Ratio(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, RigidityBinningCenter, Name, CCcut_TF):

    fig, ax = plt.subplots(figsize=(30, 18))
    plt.xlim(10,550)
    plt.ylim(0, 3)

    plt.xlabel('|R| / (GV)', fontsize=60, horizontalalignment='right', x=1.0)
    if Name == 'CCLevel':
        plt.ylabel('CCLevel ratio (ISS/MC)', fontsize=60, horizontalalignment='right', y=1.0)
    elif Name == 'CCProtonOverProton':
        plt.ylabel('CCProtonOverProton (ISS/MC)', fontsize=60, horizontalalignment='right', y=1.0)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    #plt.grid(True, which='both', axis='y')
    plt.legend(loc='best', fontsize=50)

    ## plot raw data
    #plt.errorbar( NNpoint[PointsForFit_Begin:PointsForFit_End], RatioUserd[PointsForFit_Begin:PointsForFit_End], yerr=RatioErrorUserd[PointsForFit_Begin:PointsForFit_End], xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )
    #plt.errorbar( NNpoint                                    , RatioUserd                                     , yerr=RatioErrorUserd, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )
    plt.errorbar( RigidityBinningCenter                       , RatioUserd                                     , yerr=RatioErrorUserd, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )

    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/Ratio_' + Name + '_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LinearX.pdf')
    #ax.set_xscale('log')
    #ax.set_xticks([20, 30, 50, 70, 100, 200, 400])
    #ax.get_xaxis().set_major_formatter(ScalarFormatter())
    #plt.ylim(0, 2.2)
    #plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/RatioUserd_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LogX.pdf')

    plt.close()


def Plot_ConfidenceInterval_Option2_PolyFit(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, yhat_67_PX, upper_67_PX, lower_67_PX, PolynomialOrder, RigidityBinningCenter, Name, CCcut_TF):
    fig, ax = plt.subplots(figsize=(18, 9))

    plt.xlim(10,550)
    plt.ylim(0, 3)

    plt.xlabel('|R| / (GV)',fontsize=30)
    #plt.ylabel('CCLevel/RigidityResolution',fontsize=30)
    plt.ylabel('CCLevel (ISS/MC)',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.grid(True, which='both', axis='y')
    plt.legend( loc='best',fontsize=20)

    ## plot raw data
    #plt.errorbar( NNpoint[PointsForFit_Begin:PointsForFit_End], RatioUserd[PointsForFit_Begin:PointsForFit_End], yerr=RatioErrorUserd[PointsForFit_Begin:PointsForFit_End], xerr=0, fmt='o' )
    #plt.errorbar( NNpoint                                    , RatioUserd                                     , yerr=RatioErrorUserd, xerr=0, fmt='o' )
    plt.errorbar( RigidityBinningCenter                       , RatioUserd                                     , yerr=RatioErrorUserd, xerr=0, fmt='o' ) 

    ## plot Fit and Uncertainty Band (x_more[11:] is from 60 GV)
    plot(x_more, yhat_67_PX,  c='g', lw=2,  label="P0 Fit option2" )
    ax.fill_between(x_more, upper_67_PX, lower_67_PX, color="#e7b9cf", edgecolor="", alpha=0.5, label='P0 fit: 67% CL')
    #ax.fill_between(x_more, upper_95_P0, lower_95_P0, color="blue"   , edgecolor="", alpha=0.5, label='P0 fit: 95% CL')
    #plot(RigidityBinningCenter, yhat_67_PX,  c='g', lw=2,  label="P0 Fit option2" )

    #plot(x_more[11:], yhat_67_PX[11:],  c='g', lw=2,  label="P3 Fit" )
    #ax.fill_between(x_more[11:], upper_67_PX[11:], lower_67_PX[11:], color="#e7b9cf", edgecolor="", alpha=0.5, label='P3 Fit: 67% CL')

    plt.legend( loc='best',fontsize=20)
    plt.xscale('linear')
    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_Poly4Fit_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + "_" + PolynomialOrder + '_LinearX.pdf')
    #plt.xscale('log')
    #plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_Poly4Fit_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + "_" + PolynomialOrder + '_LogX.pdf')

    plt.close()


def Plot_ConfidenceInterval_Option3_curvefit (NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, y_model_option3, x_more, y_more_option3, ParameterError, RigidityBinningCenter, Name, chi2dof, Scaler, ScalerBoolForPlot, sigma, CCcut_TF):

    fig, ax = plt.subplots()

    plt.xlim(10,550)
    plt.ylim(0, 3)

    ## plot raw data
    plt.errorbar( RigidityBinningCenter, RatioUserd, yerr=RatioErrorUserd/Scaler, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )

    ## plot Linear Fit (Option 3)
    plt.hlines(y_model_option3, x_more[0], x_more[-1], linewidth=5, alpha=0.5, linestyle='dashed', color="red", label='chi2dof=' + str(chi2dof))
    if ScalerBoolForPlot == True:
        ax.fill_between(x_more, y_more_option3 + ParameterError * Scaler, y_more_option3 - ParameterError * Scaler, color="yellow", edgecolor="", alpha=1.0, label='68% Confidence Interval') # label='67% Confidence Interval', color="#b9cfe7"
    elif ScalerBoolForPlot == False:
        ax.fill_between(x_more, y_more_option3 + ParameterError, y_more_option3 - ParameterError, color="yellow", edgecolor="", alpha=1.0, label='68% Confidence Interval') # label='67% Confidence Interval', color="#b9cfe7"

    plt.xlabel('|R| / (GV)', horizontalalignment='right', x=1.0)
    if Name == 'CCLevel':
        plt.ylabel('CCLevel ratio (ISS/MC)'     , horizontalalignment='right', y=1.0)
    elif Name == 'CCProtonOverProton':
        plt.ylabel('CCProtonOverProton (ISS/MC)', horizontalalignment='right', y=1.0)

    plt.legend(fontsize=50)

    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_Option3CurveFit_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LinearX_' + 'ScalerBoolForPlot_' + str(ScalerBoolForPlot) + '_' + sigma + '_Sigma' + '.pdf')

    plt.close()


def Plot_Option4_BandWidthFor68PercentPoints(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, BandWidthFor68PercentPoints_more, RigidityBinningCenter, Name, Scaler, ScalerBoolForPlot, y_model_option3, y_more_option3, CCcut_TF):

    fig, ax = plt.subplots()

    plt.xlim(10,550)
    plt.ylim(0, 3)

    ## plot raw data
    if ScalerBoolForPlot == False:
        plt.errorbar( RigidityBinningCenter, RatioUserd, yerr=RatioErrorUserd, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )
    elif ScalerBoolForPlot == True:
        plt.errorbar( RigidityBinningCenter, RatioUserd, yerr=RatioErrorUserd/Scaler, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )

    ## plot Linear Fit 
    plt.hlines(y_model_option3, x_more[0], x_more[-1], linewidth=5, alpha=0.5, linestyle='dashed', color="red")
    ax.fill_between(x_more, y_more_option3 + BandWidthFor68PercentPoints_more[0], y_more_option3 - BandWidthFor68PercentPoints_more[0], color="yellow", edgecolor="", alpha=1.0)

    plt.xlabel('|R| / (GV)', horizontalalignment='right', x=1.0)
    if Name == 'CCLevel':
        plt.ylabel('CCLevel ratio (ISS/MC)'     , horizontalalignment='right', y=1.0)
    elif Name == 'CCProtonOverProton':
        plt.ylabel('CCProtonOverProton (ISS/MC)', horizontalalignment='right', y=1.0)

    plt.legend(fontsize=50)

    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_Option4Manual68PercentPoints_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LinearX_' + 'ScalerBoolForPlot_' + str(ScalerBoolForPlot) + '.pdf')

    plt.close()


def Plot_Option4_2_BandWidth_WithBump(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd, RatioErrorUserd, x_more, BandWidth_UnBumpRange, RigidityBinningCenter, Name, Scaler, ScalerBoolForPlot, y_model_option3, y_more_option3, y_more_Bump, BumpRangeStart, BumpRangeEnd, CCcut_TF):

    fig, ax = plt.subplots()

    plt.xlim(10,550)
    plt.ylim(0, 3)

    ## plot raw data
    if ScalerBoolForPlot == False:
        plt.errorbar( RigidityBinningCenter, RatioUserd, yerr=RatioErrorUserd, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )
    elif ScalerBoolForPlot == True:
        plt.errorbar( RigidityBinningCenter, RatioUserd, yerr=RatioErrorUserd/Scaler, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )

    ## plot Linear Fit
    plt.hlines(y_model_option3, x_more[0], x_more[-1], linewidth=5, alpha=0.5, linestyle='dashed', color="red")
    ax.fill_between(x_more, y_more_option3 + BandWidth_UnBumpRange[0], y_more_option3 - BandWidth_UnBumpRange[0], color="yellow", edgecolor="", alpha=1.0)
    ## Bump
    #plt.plot(x_more[BumpRangeStart:BumpRangeEnd], y_more_Bump[BumpRangeStart:BumpRangeEnd] )
    ax.fill_between(x_more[BumpRangeStart:BumpRangeEnd], y_more_option3 + (y_more_option3-y_more_Bump[BumpRangeStart:BumpRangeEnd]), y_more_Bump[BumpRangeStart:BumpRangeEnd], color="yellow", edgecolor="", alpha=1.0)

    plt.xlabel('|R| / (GV)', horizontalalignment='right', x=1.0)
    if Name == 'CCLevel':
        plt.ylabel('CCLevel ratio (ISS/MC)'     , horizontalalignment='right', y=1.0)
    elif Name == 'CCProtonOverProton':
        plt.ylabel('CCProtonOverProton (ISS/MC)', horizontalalignment='right', y=1.0)

    plt.legend(fontsize=50)

    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_Option4.2Manual100PercentPointsWithBump_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LinearX_' + 'ScalerBoolForPlot_' + str(ScalerBoolForPlot) + '.pdf')

    plt.close()


def Plot_Option4_3_SmoothedBandWidth(NNsuffix, pattern, ISSversion, Binningversion, workpath, suffix, NNpoint, PointsForFit_Begin, PointsForFit_End, RatioUserd_AllPointsForNoFurtherRemoved, RatioErrorUserd_AllPointsForNoFurtherRemoved, x_more, BandWidth_UnBumpRange, RigidityBinningCenter_AllPointsForNoFurtherRemoved, Name, Scaler, ScalerBoolForPlot, y_model_option3, y_more_option3, y_more_Bump, BumpRangeStart_in_x_more, BumpRangeEnd_in_x_more, CCcut_TF, BandWidth_Option4_3_Final, BandWidth_Option4_3_FirstConst_more, BumpRange_PeakValue_in_x_more, FurtherPointRemoved, ShiftValue):

    fig, ax = plt.subplots()

    plt.xlim(10,550)
    plt.ylim(0, 3)

    ## plot raw data
    if ScalerBoolForPlot == False:
        plt.errorbar( RigidityBinningCenter_AllPointsForNoFurtherRemoved[FurtherPointRemoved:], RatioUserd_AllPointsForNoFurtherRemoved[FurtherPointRemoved:], yerr=RatioErrorUserd_AllPointsForNoFurtherRemoved[FurtherPointRemoved:], xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )
    elif ScalerBoolForPlot == True:
        plt.errorbar( RigidityBinningCenter_AllPointsForNoFurtherRemoved[FurtherPointRemoved:], RatioUserd_AllPointsForNoFurtherRemoved[FurtherPointRemoved:], yerr=RatioErrorUserd_AllPointsForNoFurtherRemoved[FurtherPointRemoved:]/Scaler, xerr=0, fmt='o', markerfacecolor="black", ecolor="black", linewidth=5 )

    ## plot Linear Fit
    plt.hlines(y_model_option3, x_more[0], x_more[-1], linewidth=5, alpha=0.5, linestyle='dashed', color="red")
    ax.fill_between(x_more, y_more_option3 + BandWidth_UnBumpRange[0], y_more_option3 - BandWidth_UnBumpRange[0], color="yellow", edgecolor="", alpha=1.0)

    ## Plot Bump Fit
    #plt.plot(x_more[BumpRangeStart_in_x_more:BumpRangeEnd_in_x_more], y_more_Bump[BumpRangeStart_in_x_more:BumpRangeEnd_in_x_more] )

    ## Plot Bump
    ax.fill_between(x_more[BumpRangeStart_in_x_more:BumpRangeEnd_in_x_more], y_more_option3 + (y_more_option3 - y_more_Bump[BumpRangeStart_in_x_more:BumpRangeEnd_in_x_more])-ShiftValue, y_more_Bump[BumpRangeStart_in_x_more:BumpRangeEnd_in_x_more]+ShiftValue, color="yellow", edgecolor="", alpha=1.0)
    print("y_more_option3: " + str(y_more_option3))
    print("BandWidth_Option4_3_FirstConst_more:" + str(BandWidth_Option4_3_FirstConst_more))

    ## Plot Const Part at beginning
    ax.fill_between(x_more[0:BumpRange_PeakValue_in_x_more]                , y_more_option3 - BandWidth_Option4_3_FirstConst_more, y_more_option3 + BandWidth_Option4_3_FirstConst_more, color="yellow", edgecolor="", alpha=1.0)

    plt.xlabel('|R| / (GV)', horizontalalignment='right', x=1.0)
    if Name == 'CCLevel':
        plt.ylabel('CCLevel ratio (ISS/MC)'     , horizontalalignment='right', y=1.0)
    elif Name == 'CCProtonOverProton':
        plt.ylabel('CCProtonOverProton (ISS/MC)', horizontalalignment='right', y=1.0)

    plt.legend(fontsize=50)

    plt.savefig(workpath + '/templatefit/negative/FixedCC/Plot_UncertaintyBand/UncertaintyBand_' + Name + '_Option4.3_Pattern_' + pattern + NNsuffix + "_CCcut_TF_" + str(CCcut_TF) + "_" + Binningversion + suffix + '_LinearX_' + 'ScalerBoolForPlot_' + str(ScalerBoolForPlot) + '.pdf')

    plt.close()



#### others

def ListToVector(rawlist):
    vector = ROOT.vector('double')(len(rawlist))
    vector.clear()
    for i in range(len(rawlist)):
        vector.insert(vector.begin()+i, rawlist[i])
    return vector


