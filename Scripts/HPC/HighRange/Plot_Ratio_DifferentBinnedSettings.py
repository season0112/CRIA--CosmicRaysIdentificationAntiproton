#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TTree 
from root_numpy import fill_hist, array2tree
import binning
import PythonPlotDefaultParameters


def main():
    #### publised result
    #from 1 to 450 GV, 57 bins.
    rigidity_bin_number_published = 57 ## 1-450 GV
    antiproton_number_published = binning.antiproton_number_published_full
    proton_number_published = binning.proton_number_published_full
    rigiditybin_published = binning.published2016binnings


    #### Load ROOT File Result

    f_FitResult_negative_CCFree  = TFile.Open(workpath + "/templatefit/negative/FitResult/FitResult_Pattern_" + str(pattern) + "_" + str(ISSversion) + "_" + str(Binningversion) + "_CCFree.root", "READ")
    f_FitResult_positive_CCFree  = TFile.Open(workpath + "/templatefit/positive/FitResult/FitResult_Pattern_" + str(pattern) + "_" + str(ISSversion) + "_" + str(Binningversion) + "_CCFree.root", "READ")

    ## Old
    #CCNumber_All = [20, 20, 9, 20]
    #TRDNumber_All =[12, 20, 11, 16]
    #CCcut_TF_All = [0.00]
    #CCcut_TF_All = np.append(CCcut_TF_All, np.arange(0.20, 0.50, 0.05))  # 0.50 in not included
    ## New
    #CCNumber_All = np.array([9, 11, 13, 19, 20, 21, 25, 27, 29])
    #TRDNumber_All = np.array([9, 11, 13, 19, 20, 21, 25, 27, 29])
    #CCcut_TF_All = np.array([0.00, 0.20, 0.35, 0.40, 0.65, 0.70, 0.80, 0.90])
    ## New2
    #CCNumber_All = np.array([9, 13, 20, 27])
    #TRDNumber_All = np.array([9, 13, 20, 27])
    #CCcut_TF_All = np.array([ 0.20, 0.35, 0.70])
    ## Current Using:
    CCcut_TF_All = np.array([0.00, 0.20, 0.35, 0.40, 0.65, 0.70, 0.80, 0.90])
    CCNumber     = np.array([30, 20, 20, 20, 9])
    TRDNumber    = np.array([20, 20, 16, 12, 11])


    #### template fit with different settings.
    MinMeanFitchi2 = 999

    for CCcut_TF in CCcut_TF_All:
        for ccnumber, trdnumber in zip(CCNumber, TRDNumber): # Previous Goal: last two points different.

        #for ccnumber in CCNumber_All:
            #for trdnumber in TRDNumber_All:

                antiproton_number    = np.asarray(f_FitResult_negative_CCFree.Get("AntiprotonNumber_cccut_" + str(format(CCcut_TF, '.2f')) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
                proton_number_TF     = np.asarray(f_FitResult_positive_CCFree.Get("AntiprotonNumber_cccut_" + str(format(CCcut_TF, '.2f')) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
                delta_antiproton     = np.asarray(f_FitResult_negative_CCFree.Get("AntiprotonNumberError_cccut_" + str(format(CCcut_TF, '.2f')) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
                fitchi2              = np.asarray(f_FitResult_negative_CCFree.Get("chi2dof_cccut_" + str(format(CCcut_TF, '.2f')) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))

                if abs(MinMeanFitchi2-1) > abs(np.mean(fitchi2)-1):
                    print("Updating best fit parameters:")
                    MinMeanFitchi2 = np.mean(fitchi2)
                    print("MinMeanFitchi2: " + str(MinMeanFitchi2))
                    print("CCcut_TF:"  + str(CCcut_TF))
                    print("ccnumber:"  + str(ccnumber))
                    print("trdnumber:" + str(trdnumber))
                    #print("statistic_error: " + str(statistic_error[-1:]))


                ratio = antiproton_number/proton_number_TF
                ratio_published = antiproton_number_published/proton_number_published
                publishedPointCenter = binning.published2016binnings_center
                statistic_error = delta_antiproton/proton_number_TF


                #### Plot
                plt.figure(figsize=(32,18))
                ax = plt.gca()

                plt.errorbar(rigiditybin_center  , ratio          , yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="blue" , ecolor="blue", label='This analysis')
                plt.errorbar(publishedPointCenter, ratio_published, yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="red", ecolor="red", label='2016PRLpaper' )

                plt.grid(True)
                ax.set_xlim([10, 600])
                ax.set_ylim([0.00001, 0.0003])
                plt.xticks(fontsize=50)
                plt.yticks(fontsize=50)
                plt.xscale("log")
                #plt.yscale("log")
                plt.xlabel('Rigidity (GV)', fontsize=60)
                plt.ylabel('Ratio'        , fontsize=60)
                plt.legend(loc='lower left', prop={'size': 40})
                plt.savefig("Ratio_Pattern" + pattern + "_cccut_" + str(format(CCcut_TF, '.2f')) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + suffix + ".pdf")
                plt.close()

                plt.figure(figsize=(32,18))
                plt.errorbar(rigiditybin_center, statistic_error, yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="blue", ecolor="blue", label='This analysis')
                plt.xticks(fontsize=50)
                plt.yticks(fontsize=50)
                plt.grid(True)
                plt.legend(loc='upper left', prop={'size': 40})
                plt.savefig("StaError_Pattern" + pattern + "_cccut_" + str(format(CCcut_TF, '.2f')) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + suffix + ".pdf")
                plt.close()

     
                plt.figure(figsize=(32,18))
                plt.errorbar(rigiditybin_center, fitchi2, yerr=0, marker='o', linestyle="None", markersize=25, markerfacecolor="blue", ecolor="blue", label='This analysis')
                plt.xlabel('Rigidity (GV)', fontsize=60)
                plt.ylabel('Chi2'        , fontsize=60)
                plt.legend(loc='lower right', prop={'size': 40})
                plt.xticks(fontsize=50)
                plt.yticks(fontsize=50)
                plt.xscale('log')
                plt.yscale('log')
                plt.savefig("Chi2_Pattern" + pattern + "_cccut_" + str(format(CCcut_TF, '.2f')) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + suffix + ".pdf")
                plt.close()

    f_FitResult_negative_CCFree.Close()
    f_FitResult_positive_CCFree.Close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--issversion', help='ISS data version, pass7.8 or published2016')
    parser.add_argument('--binningversion', help='binningversion, 450version or 525version')
    parser.add_argument('--pattern', help='which tracker patterns you choose')
    arguments = parser.parse_args()

    if (arguments.issversion):
        ISSversion = arguments.issversion
    else:
        print("You need to choose a ISS data cersion!")
        os._exit(0)

    if (arguments.binningversion):
        Binningversion = arguments.binningversion
    else:
        print("You need to choose a binning version!")
        os._exit(0)

    if (arguments.pattern):
        pattern = arguments.pattern
    else:
        print("You need to choose a tracker pattern !")
        os._exit(0)

    if ISSversion == "pass7.8":
        suffix = ""
    elif ISSversion == "published2016":
        suffix = "_May2015"
    elif ISSversion == "PhyRep2021":
        suffix = "_Nov2017"

    if Binningversion == "525version":
        resultdir = "results_525version"
        resultdir_uncertainty = "results_525version_uncertainty"
    elif Binningversion == "450version":
        resultdir = "results_450version"
        resultdir_uncertainty = "results_450version_uncertainty"

    if Binningversion == "450version":
        rigidity_bin_number = 31 ## 14.1-450 GV
        rigiditybin = binning.published2016binnings[26:58] # 14.1-450
        rigiditybin_center = binning.published2016binnings_center[26:58]
    elif Binningversion == "525version":
        rigidity_bin_number = 32 ## 14.1-525 GV
        #rigiditybin = binning.Newbinnings_525[26:59] # 14.1-525
        #rigiditybin_center = binning.Newbinnings_525_center[26:59]
        rigiditybin = binning.Newbinnings_525_zhili[26:59]
        rigiditybin_center = binning.Newbinnings_525_center_zhili[26:59]


    workpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis/data'
    resultpath = os.getenv('HPCHIGHENERGYDATADIR') + '/ISS_anylsis'

    main()


