import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
import ROOT
from root_numpy import root2array, tree2array
import PythonPlotDefaultParameters
from scipy.optimize import curve_fit
import binning
from iminuit import Minuit
from iminuit.cost import LeastSquares
from scipy.optimize import minimize
import uncertainties
from uncertainties import unumpy
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import ScalarFormatter


def Plot_RigidityResolution(workpath, BinningCenter, Resolution, ResolutionError, p_ResModel, pattern, ISSversion, Binningversion, NNsuffix):
    plt.figure(figsize=(30,18))
    ax=plt.gca()

    plt.xlabel('Rigidity',fontsize=60)
    plt.ylabel(r'$\sigma$ of the gaussian fit',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)

    plt.xlim(10,550)
    ax.set_xscale('log')
    ax.set_xticks([20, 30, 50, 70, 100, 200, 400])
    ax.get_xaxis().set_major_formatter(ScalarFormatter())

    ## Plot RigidityResolution without Fit
    plt.errorbar(BinningCenter[0:], Resolution, yerr=ResolutionError, fmt='o', markersize=15, color="black", label='Resolution')
    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "Resolution_WithoutFit_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Linear.pdf")
    plt.yscale('log')
    plt.savefig( workpath + '/templatefit/negative/' + 'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "Resolution_WithoutFit_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Log.pdf")

    ## Plot RigidityResolution with LinearFit
    plt.yscale('linear')
    plt.plot( np.linspace(14, 525, 100), 10**p_ResModel(np.log10(np.linspace(14, 525, 100))), '--', label='ResolutionModel')
    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "Resolution_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Linear.pdf")
    plt.yscale('log')
    ax.tick_params(axis='y', direction='in', which='both', labelsize=40)
    plt.savefig( workpath + '/templatefit/negative/' + 'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "Resolution_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Log.pdf")
    plt.close()


def Plot_CCLevelISS(workpath, BinningCenter, CCLevelISS, CCLevelISSError, p_CCLevelModel, pattern, ISSversion, Binningversion, NNsuffix):
    plt.figure(figsize=(30,18))
    ax=plt.gca()

    plt.xlabel('Rigidity',fontsize=60)
    plt.ylabel('CCLevel',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)

    plt.xlim(10,550)
    ax.set_xscale('log')
    ax.set_xticks([20, 30, 50, 70, 100, 200, 400])
    ax.get_xaxis().set_major_formatter(ScalarFormatter())

    ## Plot CCLevel without Fit
    plt.errorbar(BinningCenter[0:], CCLevelISS, yerr=CCLevelISSError, fmt='o', markersize=15, color="black", label='CCLevel')
   
    ylim = ax.get_ylim() 
    minvalue = ylim[0]
    maxvalue = ylim[1]
    if str(pattern) == "4":
        plt.vlines(38.9, minvalue, maxvalue, linestyles = "solid")
    if str(pattern) == "1":
        plt.vlines(147,  minvalue, maxvalue, linestyles = "solid")
    if str(pattern) == "2":
        plt.vlines(175,  minvalue, maxvalue, linestyles = "solid")
    
    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "CCLevelWithoutFit_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Linear.pdf")

    plt.yscale('log')
    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "CCLevelWithoutFit_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Log.pdf")

    ## Plot CLevel with LinearFit in Log
    plt.yscale('linear')
    plt.plot( np.linspace(14, 525, 100), 10**p_CCLevelModel(np.log10(np.linspace(14, 525, 100))), '--', label='CCLevelModel')
    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "CCLevelFit_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Linear.pdf")
    plt.yscale('log')
    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "CCLevelFit_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Log.pdf")
    plt.close()


def Plot_CCLevelISSOverRigidityResolution(workpath, pattern, ISSversion, Binningversion, BinningCenter, CCLevelISS, CCLevelISSError, Resolution_Rescale, ResolutionError_Rescale, ResolutionFitLine_Rescale, p_CCLevelModel, p_ResModel, NNsuffix):
    plt.figure(figsize=(30,18))
    ax=plt.gca()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Rigidity',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.errorbar(BinningCenter[0:], CCLevelISS, yerr=CCLevelISSError, fmt='o', markersize=10, color="black", label='CCLevel')
    plt.errorbar(BinningCenter[0:], Resolution_Rescale, yerr=ResolutionError_Rescale, fmt='o', markersize=10, color="green", label='Resolution_Rescale')
    #plt.errorbar(BinningCenter[0:], ResolutionFitLine_Rescale, yerr=ResolutionError_Rescale, fmt='o', markersize=10, color="red", label='ResolutionFitLine_Rescale')
    plt.plot( np.linspace(14, 525, 200), 10**np.poly1d(np.polyfit(np.log10(BinningCenter[0:]), np.log10(ResolutionFitLine_Rescale), 1) )(np.log10(np.linspace(14, 525, 200))), '-', label='Resolution Linear Model (After Rescale)')
    plt.plot( np.linspace(14, 525, 100), 10**p_CCLevelModel(np.log10(np.linspace(14, 525, 100))), '--', label='CCLevel Linear Model')
    plt.plot( np.linspace(14, 525, 100), 10**p_ResModel(np.log10(np.linspace(14, 525, 100))), '--', label='Resolution Linear Model (Before Rescale)')
    plt.legend()
    plt.ylabel('CCLevel',fontsize=60)
    plt.savefig(workpath + '/templatefit/negative/' + 'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "Both_Pattern_" + pattern + NNsuffix + "_" +ISSversion + "_" + Binningversion + "_Log.pdf")
    plt.close()


def Plot_CCLevelISS_AllPatterns(workpath, BinningCenter, ISSversion, Binningversion, NNsuffix):
    ## Load CCLevel forr all patterns
    f_CCLevelISS_P0    = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + "0" + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")
    CCLevelISS_P0      = np.asarray(f_CCLevelISS_P0.Get("CCLevelISS_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    CCLevelISSError_P0 = np.asarray(f_CCLevelISS_P0.Get("CCLevelISS_Error_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    f_CCLevelISS_P0.Close()

    f_CCLevelISS_P1    = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + "1" + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")
    CCLevelISS_P1      = np.asarray(f_CCLevelISS_P1.Get("CCLevelISS_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    CCLevelISSError_P1 = np.asarray(f_CCLevelISS_P1.Get("CCLevelISS_Error_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    f_CCLevelISS_P1.Close()

    f_CCLevelISS_P2    = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + "2" + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")
    CCLevelISS_P2      = np.asarray(f_CCLevelISS_P2.Get("CCLevelISS_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    CCLevelISSError_P2 = np.asarray(f_CCLevelISS_P2.Get("CCLevelISS_Error_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    f_CCLevelISS_P2.Close()

    f_CCLevelISS_P4    = TFile.Open(workpath + "/templatefit/negative/FixedCC/CCLevel_Pattern_" + "4" + "_" + str(ISSversion) + "_" + str(Binningversion) + ".root", "READ")
    CCLevelISS_P4      = np.asarray(f_CCLevelISS_P4.Get("CCLevelISS_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    CCLevelISSError_P4 = np.asarray(f_CCLevelISS_P4.Get("CCLevelISS_Error_Signal_Efficiency_1.0_CCcut_TF_0.00"))
    f_CCLevelISS_P4.Close()


    #### Fix the first a few points due to limited statistics.
    CCLevelISS_P0[0]      = CCLevelISS_P0[1] * 0.6
    CCLevelISSError_P0[0] = CCLevelISSError_P0[1] * 1.8


    ## Plot
    plt.figure(figsize=(30,18))
    ax=plt.gca()

    plt.xlabel('Rigidity',fontsize=60)
    plt.ylabel('CCLevel',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)

    plt.xlim(10,550)
    ax.set_xscale('log')
    ax.set_xticks([20, 30, 50, 70, 100, 200, 400])
    ax.get_xaxis().set_major_formatter(ScalarFormatter())

    ## Plot CCLevel for Used Range
    plt.errorbar(BinningCenter[0:], CCLevelISS_P0, yerr=CCLevelISSError_P0, fmt='o', markersize=20, color='r', label='CCLevel (P0)')
    plt.errorbar(BinningCenter[0:27], CCLevelISS_P1[0:27], yerr=CCLevelISSError_P1[0:27], fmt='o', markersize=20, color='b', label='CCLevel (P1)')
    plt.errorbar(BinningCenter[0:28], CCLevelISS_P2[0:28], yerr=CCLevelISSError_P2[0:28], fmt='o', markersize=20, color='g', label='CCLevel (P2)')
    plt.errorbar(BinningCenter[0:13], CCLevelISS_P4[0:13], yerr=CCLevelISSError_P4[0:13], fmt='o', markersize=20, color='m', label='CCLevel (P4)')

    ## Plot CCLevel for NOT Used Range
    plt.errorbar(BinningCenter[27:], CCLevelISS_P1[27:], yerr=CCLevelISSError_P1[27:], marker='o', markersize=20, markerfacecolor='none', markeredgecolor='b', linewidth=0, markeredgewidth='6')
    plt.errorbar(BinningCenter[28:], CCLevelISS_P2[28:], yerr=CCLevelISSError_P2[28:], marker='o', markersize=20, markerfacecolor='none', markeredgecolor='g', linewidth=0, markeredgewidth='6')
    plt.errorbar(BinningCenter[13:], CCLevelISS_P4[13:], yerr=CCLevelISSError_P4[13:], marker='o', markersize=20, markerfacecolor='none', markeredgecolor='m', linewidth=0, markeredgewidth='6')

    print("BinningCenter is" )
    print(BinningCenter)
    print(BinningCenter[12]) ## BinningCenter[12] = 37.5 (36.1_38.9), 
    print(BinningCenter[26]) ## BinningCenter[26] = 136.0 (125_147)
    print(BinningCenter[27]) ## BinningCenter[27] = 161.0 (147_175)

    ylim = ax.get_ylim()
    minvalue = ylim[0]
    maxvalue = ylim[1]
    plt.vlines(38.9, minvalue, maxvalue, linestyles = "solid")
    plt.vlines(147,  minvalue, maxvalue, linestyles = "solid")
    plt.vlines(175,  minvalue, maxvalue, linestyles = "solid")

    plt.legend(loc='best', prop={'size': 35})

    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "CCLevelWithoutFit_ALLPatterns_" + ISSversion + "_" + Binningversion + "_Linear.pdf")
    plt.yscale('log')
    plt.savefig(workpath+'/templatefit/negative/'+'FixedCC/Plot_ISSCCLevelOverRigidityResolution/' + "CCLevelWithoutFit_AllPatterns_" + ISSversion + "_" + Binningversion + "_Log.pdf")







def ListToVector(rawlist):
    vector = ROOT.vector('double')(len(rawlist))
    vector.clear()
    for i in range(len(rawlist)):
        vector.insert(vector.begin()+i, rawlist[i])
    return vector
