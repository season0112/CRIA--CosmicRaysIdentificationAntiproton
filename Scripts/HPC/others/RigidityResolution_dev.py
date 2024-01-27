#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import root2array, tree2array
import PythonPlotDefaultParameters
from scipy.optimize import curve_fit
import binning
import ROOT
from matplotlib.ticker import FormatStrFormatter


def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def Plot_RigidityResolution(workpath, BinningCenter, Pattern, Binningversion, SigmaAll, SigmaErrorAll):
    plt.figure(figsize=(30,18))
    ax=plt.gca()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Rigidity',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.errorbar(BinningCenter, SigmaAll, yerr=SigmaErrorAll, label="", fmt='o', markersize=15, color="black")
    plt.savefig(workpath + '/templatefit/negative/' + 'FixedCC/Plot_RigidityResolution/' + "RigidityResolution_Patttern_" + str(Pattern) + "_" + str(Binningversion) + ".pdf")


def Plot_Chi2dof(workpath, BinningCenter, Pattern, Binningversion, Chi2dofAll):
    plt.figure(figsize=(30,18))
    ax=plt.gca()
    #plt.yscale('log')
    ax.axes.set_ylim( [0, 5] )
    plt.xlabel('Rigidity',fontsize=60)
    plt.ylabel('Chi2dof',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)

    plt.tick_params(axis='x', which='minor', bottom=True, top=False, labelbottom=True, direction='in', length=9, width=7, labelsize=27, colors='black')
    plt.tick_params(axis='x', which='major',  bottom=True, top=False, labelbottom=True, direction='in', length=12, width=7, labelsize=33, colors='black')
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))

    plt.errorbar(BinningCenter, Chi2dofAll, yerr=0, label="Chi2/dof", fmt='o', markersize=15, color="black")

    plt.savefig(workpath + '/templatefit/negative/' + 'FixedCC/Plot_RigidityResolution/' + "Chi2dof_Patttern_" + str(Pattern) + "_" + str(Binningversion) + "_LinearX.pdf")
    plt.xscale('log')
    plt.tick_params(axis='x', which='minor', bottom=True, top=False, labelbottom=True, direction='in', length=9, width=7, labelsize=27, colors='black')
    plt.tick_params(axis='x', which='major',  bottom=True, top=False, labelbottom=True, direction='in', length=12, width=7, labelsize=33, colors='black')
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    plt.savefig(workpath + '/templatefit/negative/' + 'FixedCC/Plot_RigidityResolution/' + "Chi2dof_Patttern_" + str(Pattern) + "_" + str(Binningversion) + "_LogX.pdf")


def Plot_FitCore(workpath, BinningEdgeName, Pattern, Binningversion, center, hist, popt, i, reduced_chi_squared):
    plt.figure(figsize=(30,18))
    ax=plt.gca()
    plt.yscale('log')
    plt.xlabel(r'$R_T$/$R_M$-1',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.errorbar(center, hist, yerr=np.sqrt(hist), label="", fmt='o', markersize=15, color="black")
    #plt.plot(center[FitStart:FitEnd], gauss_function(center[FitStart:FitEnd], *popt), 'ro:', label='fit')
    #plt.ylim(0.1, np.max(hist)*1.2)
    plt.text(0.77, 0.9, 'Chi2/dof='+str("{:.2f}".format(reduced_chi_squared)), transform=ax.transAxes, fontsize=40)
    plt.savefig( workpath + '/templatefit/negative/' + 'FixedCC/Plot_RigidityResolution/' + "Fit_Core_Patttern_" + str(Pattern) + "_" + str(BinningEdgeName[i]) + ".pdf")
    plt.close()


def Plot_1overR_M(workpath, BinningEdgeName, Pattern, Binningversion, center_1OverR_M, hist_1OverR_M, i):
    plt.figure(figsize=(30,18))
    ax=plt.gca()
    plt.yscale('log')
    #plt.xlabel(r'$R_T$/$R_M$-1',fontsize=60)
    plt.xlabel(r'$1/R_M$',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    plt.errorbar(center_1OverR_M, hist_1OverR_M, yerr=np.sqrt(hist_1OverR_M), label="", fmt='o', markersize=15, color="black")
    #plt.plot(center_res[FitStart:FitEnd], gauss_function(center_res[FitStart:FitEnd], *popt_res), 'ro:', label='fit')
    #plt.ylim(0.1, np.max(hist_res)*1.2)
    #plt.text(0.77, 0.9, 'Chi2/dof='+str("{:.2f}".format(reduced_chi_squared_res)), transform=ax.transAxes, fontsize=40)
    plt.savefig( workpath + '/templatefit/negative/' + 'FixedCC/Plot_RigidityResolution/' + "1OverR_M_Patttern_" + str(Pattern) + "_" + str(BinningEdgeName[i]) + ".pdf")
    plt.close()


def ListToVector(rawlist):
    vector = ROOT.vector('double')(len(rawlist))
    vector.clear()
    for i in range(len(rawlist)):
        vector.insert(vector.begin()+i, rawlist[i])
    return vector


def main(FitEnd):
    SigmaAll      = []
    SigmaErrorAll = []
    Chi2dofAll    = []

    index = 0

    for i in range(BinningEdgeName.shape[0]):
        ## Load data
        protonMCname = "B1042_pr.pl1.1800_7.6_all"
        #protonMCname = "B1042_pr.pl1.flux.l1a9.2016000_7.6_all"
        protondata_positive = root2array(highpath + "/" + protonMCname + "_Tree_positive_" + str(BinningEdgeName[i]) + ".root", "ExampleAnalysisTree", selection="Pattern==" + str(Pattern))
        #protondata_negative = root2array(highpath + "/" + protonMCname + "_Tree_negative_" + str(BinningEdgeName[i]) + ".root", "ExampleAnalysisTree", selection="Pattern==" + str(Pattern))
        #protondata_all = np.append(protondata_positive, protondata_negative)
        protondata_all = protondata_positive

        ## Define Rigidity Resolution ((1/R_m - 1/R_t) / (1/R_t) = R_t / R_m - 1)
        #resolution = protondata_all['MCPrimaryMomentum']/protondata_all['Rigidity'] - 1
        #resolution = 1/protondata_all['Rigidity'] - 1/protondata_all['MCPrimaryMomentum']

        One_Over_R_M = 1/protondata_all['Rigidity']
        One_Over_R_T = 1/protondata_all['MCPrimaryMomentum']

        ## Fill Histogram
        #hist, bins = np.histogram(resolution, bins=50, range=(-1, 1) )
        #hist, bins = np.histogram(resolution, bins=50,  )
        
        #hist, bins = np.histogram2d(One_Over_R_M, One_Over_R_T, bins=50)     
        H, xedges, yedges = np.histogram2d(One_Over_R_M, One_Over_R_T, bins=50)

        center = (bins[:-1] + bins[1:]) / 2

        ## Fit the core 
        
        if Pattern == '0':
            if index > 21:
                FitEnd = 30
        

        print('\n')
        print("Fit rigidity bin is: " + str(BinningEdgeName[i]))
        print("Actual Fit Range: " + str(center[FitStart]) + " to " + str(center[FitEnd]))
        popt, pcov = curve_fit(gauss_function, center[FitStart:FitEnd], hist[FitStart:FitEnd], p0 = [2.7, 6.3, 1.2], maxfev=50000000  ) 

        chi_squared = np.sum( (gauss_function(center[FitStart:FitEnd], *popt) - hist[FitStart:FitEnd])**2 / gauss_function(center[FitStart:FitEnd]**2, *popt) ) # Check: Chi2
        reduced_chi_squared = chi_squared/( (FitEnd-FitStart) - len(popt) )
        print("chi2/dof: " + str(reduced_chi_squared) )

        ## Plot Fit Core
        Plot_FitCore(workpath, BinningEdgeName, Pattern, Binningversion, center, hist, popt, i, reduced_chi_squared)

        SigmaAll.append(abs(popt[2]))
        SigmaErrorAll.append(np.sqrt(np.diag(pcov))[2])
        Chi2dofAll.append(reduced_chi_squared)

        ## index
        index = index + 1

    '''
    #### Save result
    f_Resolution = TFile.Open(workpath + "/templatefit/negative/FixedCC/RigidityResolution_Pattern_" + str(Pattern) + "_" + str(Binningversion) + ".root", "RECREATE")
    f_Resolution.WriteObject( ListToVector(SigmaAll), "Sigma") 
    f_Resolution.WriteObject( ListToVector(SigmaErrorAll), "SigmaError") 
    f_Resolution.Close()
    '''

    #### Plot RigidityResolution
    Plot_RigidityResolution(workpath, BinningCenter, Pattern, Binningversion, SigmaAll, SigmaErrorAll)
    Plot_Chi2dof(workpath, BinningCenter, Pattern, Binningversion, Chi2dofAll)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--pattern', help='Which tracker patterns you choose')
    parser.add_argument('--binningversion', help='binningversion, 450version or 525version')
    arguments = parser.parse_args()

    if (arguments.pattern):
        Pattern = arguments.pattern
    else:
        print("You need to choose a tracker pattern !")
        os._exit(0)

    if (arguments.binningversion):
        Binningversion = arguments.binningversion
    else:
        print("You need to choose a binning version!")
        os._exit(0)

    highpath = os.getenv('HPCHIGHENERGYDATADIR')
    workpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis/data'

    if Binningversion == "525version":
        BinningEdgeName = binning.bins_525_zhili
        BinningCenter = binning.Newbinnings_525_center_zhili[26:]
    elif Binningversion == "450version":
        BinningEdgeName = binning.bins_450
        BinningCenter = binning.published2016binnings_center[26:]

    BinningEdgeName = BinningEdgeName[-9:]
    BinningCenter   = BinningCenter[-9:]

    if Pattern == "0":
        ## Old setting for v2.0
        #FitStart = 10
        #FitEnd = 35
        ## Test
        FitStart = 10
        FitEnd = 32
    elif Pattern == "1":
        FitStart = 12
        FitEnd = 28
    elif Pattern == "2":
        FitStart = 10
        FitEnd = 35
    elif Pattern == "4":
        FitStart = 8
        FitEnd = 40
    elif Pattern == "5":
        FitStart = 10
        FitEnd = 35
    elif Pattern == "3":
        FitStart = 10
        FitEnd = 35
    else:
        print("Please check the tracker patterns!")

    main(FitEnd)




