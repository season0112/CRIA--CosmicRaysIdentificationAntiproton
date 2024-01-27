#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TTree 
from root_numpy import fill_hist, array2tree
import binning

def main():
    #### Load template fit result
    ## all but not last 2 
    #CCNumber_general = "20"
    #TRDNumber_general = "12"
    #CCNumber_general = "20"
    #TRDNumber_general = "20"
    ## last 2 points
    #CCNumber_last2 = "9"
    #TRDNumber_last2 = "11"
    #CCNumber_last2 = "20"
    #TRDNumber_last2 = "16"

    #CCNumber = [20, 20, 9, 20]
    #TRDNumber =[12, 20, 11, 16]

    #### Setting: (for trdrange:0-1.6) 
    if pattern == "0":
        CCNumber  = [9]
        TRDNumber = [11]
        CCcut_TF  = "0.65"
    if pattern == "1":
        CCNumber  = [9]
        TRDNumber = [11]
        CCcut_TF  = "0.65"
    elif pattern == "2":
        CCNumber  = [20]
        TRDNumber = [12]
        CCcut_TF  = "0.65"
    elif pattern == "3":
        CCNumber  = [20]
        TRDNumber = [20]
        CCcut_TF  = "0.20"
    elif pattern == "4":
        CCNumber  = [20]
        TRDNumber = [20]
        CCcut_TF  = "0.20"
    elif pattern == "5":
        CCNumber  = [20]
        TRDNumber = [20]
        CCcut_TF  = "0.20"
    elif pattern == "-1":
        CCNumber  = [20]
        TRDNumber = [20]
        CCcut_TF  = "0.20"

    f_FitResult_negative_CCFree  = TFile.Open(workpath + "/templatefit/negative/FitResult/FitResult_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + "_CCFree.root", "READ")
    f_FitResult_positive_CCFree  = TFile.Open(workpath + "/templatefit/positive/FitResult/FitResult_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + "_CCFree.root", "READ")
    f_FitResult_negative_CCFixed = TFile.Open(workpath + "/templatefit/negative/FitResult/FitResult_Pattern_" + str(pattern) + NNsuffix + "_" + str(ISSversion) + "_" + str(Binningversion) + "_CCFixed.root", "READ")

    for ccnumber, trdnumber in zip(CCNumber, TRDNumber): # Previous Goal: last two points different.
        antiproton_number    = np.asarray(f_FitResult_negative_CCFree.Get("AntiprotonNumber_cccut_" + str(CCcut_TF) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
        a_number_uncertainty = np.asarray(f_FitResult_negative_CCFixed.Get("AntiprotonNumber_cccut_" + str(CCcut_TF) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
        proton_number_TF     = np.asarray(f_FitResult_positive_CCFree.Get("AntiprotonNumber_cccut_" + str(CCcut_TF) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
        delta_antiproton     = np.asarray(f_FitResult_negative_CCFree.Get("AntiprotonNumberError_cccut_" + str(CCcut_TF) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
        fitchi2              = np.asarray(f_FitResult_negative_CCFree.Get("chi2dof_cccut_" + str(CCcut_TF) + "_CCN_" + str(ccnumber) + "_TRDN_" + str(trdnumber) + "_ISSVersion" + suffix))
        delta_antiproton_system_CC = np.abs(a_number_uncertainty - antiproton_number)

    f_FitResult_negative_CCFree.Close()
    f_FitResult_positive_CCFree.Close()
    f_FitResult_negative_CCFixed.Close()


    #### produce tree of statistic_error and delta antiproton.
    statistic_error = delta_antiproton/proton_number_TF
    statistic_error.dtype = [('statistic_error','double')]
    t_statistic_error = array2tree(statistic_error, name='tstatistic_error')

    delta_antiproton.dtype = [('delta_antiproton','double')]
    t_delta_antiproton = array2tree(delta_antiproton, name='tdelta_antiproton')

    system_CC = delta_antiproton_system_CC/proton_number_TF
    system_CC.dtype = [('system_CC','double')]
    t_system_CC = array2tree(system_CC, name='tsystem_CC')


    #### publised result
    #from 1 to 450 GV, 57 bins. 
    rigidity_bin_number_published = 57 ## 1-450 GV
    antiproton_number_published = binning.antiproton_number_published_full
    proton_number_published = binning.proton_number_published_full
    rigiditybin_published = binning.published2016binnings
    errorbar_published = binning.publishederrorbar
    rigiditypoint_published = np.array( np.zeros(rigiditybin_published.shape[0]-1) )

    errorbar_published.dtype = [('error','float32')]
    t_errorbar_published = array2tree(errorbar_published, name='terrorbar_published')


    #### Make root histogram for unfolding
    rigiditypoint = np.array( np.zeros(rigiditybin.shape[0]-1) )

    g_fitchi2 = TGraph(rigidity_bin_number, rigiditybin_center, fitchi2);

    h_proton_number = TH1D("proton_number","", rigidity_bin_number, rigiditybin)
    for i in range(rigiditybin.shape[0]-1):
        rigiditypoint[i] = (rigiditybin[i] + rigiditybin[i+1])/2
        b = np.ones((1,int(proton_number_TF[i]))) * rigiditypoint[i]
        fill_hist(h_proton_number, b[0])

    h_antiproton_number = TH1D("antiproton_number","", rigidity_bin_number, rigiditybin)
    for i in range(rigiditybin.shape[0]-1):
        rigiditypoint[i] = (rigiditybin[i] + rigiditybin[i+1])/2
        b = np.ones((1,int(antiproton_number[i]))) * rigiditypoint[i]
        fill_hist(h_antiproton_number, b[0])

    h_MeasuringTime = TH1D("MeasuringTime","", rigidity_bin_number, rigiditybin)
    for i in range(rigiditybin.shape[0]-1):
        rigiditypoint[i] = (rigiditybin[i] + rigiditybin[i+1])/2
        b = np.ones((1,1)) * rigiditypoint[i]
        fill_hist(h_MeasuringTime, b[0])

    h_Acceptance = TH1D("Acceptance","", rigidity_bin_number, rigiditybin)
    for i in range(rigiditybin.shape[0]-1):
        rigiditypoint[i] = (rigiditybin[i] + rigiditybin[i+1])/2
        b = np.ones((1,1)) * rigiditypoint[i]
        fill_hist(h_Acceptance, b[0])

    h_TriggerEfficiency = TH1D("TriggerEfficiency","", rigidity_bin_number, rigiditybin)
    for i in range(rigiditybin.shape[0]-1):
        rigiditypoint[i] = (rigiditybin[i] + rigiditybin[i+1])/2
        b = np.ones((1,1)) * rigiditypoint[i]
        fill_hist(h_TriggerEfficiency, b[0])

    h_proton_number_published = TH1D("proton_number_published","", rigidity_bin_number_published, rigiditybin_published)
    for i in range(rigiditybin_published.shape[0]-1):
        rigiditypoint_published[i] = (rigiditybin_published[i] + rigiditybin_published[i+1])/2
        b = np.ones((1,int(proton_number_published[i]))) * rigiditypoint_published[i]
        fill_hist(h_proton_number_published, b[0])

    h_antiproton_number_published = TH1D("antiproton_number_published","", rigidity_bin_number_published, rigiditybin_published)
    for i in range(rigiditybin_published.shape[0]-1):
        rigiditypoint_published[i] = (rigiditybin_published[i] + rigiditybin_published[i+1])/2
        b = np.ones((1,int(antiproton_number_published[i]))) * rigiditypoint_published[i]
        fill_hist(h_antiproton_number_published, b[0])


    #### Save result in ROOT file
    ROOTFile  = TFile(resultpath + "/unfolding/RawRatio_Pattern_" + str(pattern) + NNsuffix + "_" + str(Binningversion) + suffix + ".root","RECREATE")
    g_fitchi2.Write("g_fitchi2")

    h_proton_number.Write()
    h_antiproton_number.Write()
    #h_Sys_CC.Write()

    h_MeasuringTime.Write()
    h_Acceptance.Write()
    h_TriggerEfficiency.Write()

    h_proton_number_published.Write()
    h_antiproton_number_published.Write()

    t_errorbar_published.Write()
    t_statistic_error.Write()
    t_delta_antiproton.Write()
    t_system_CC.Write()

    ROOTFile.Close()

    

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--issversion', help='ISS data version, pass7.8 or published2016')
    parser.add_argument('--binningversion', help='binningversion, 450version or 525version')
    parser.add_argument('--pattern', help='which tracker patterns you choose')
    parser.add_argument('--ifVGGNN', help='if you want to turn to VGGNN, please turn to Yes.')
    arguments = parser.parse_args()

    if (arguments.ifVGGNN):
        IfVGGNN = arguments.ifVGGNN
    else:
        print("You need to choose which CC estimator you want to use!")
        os._exit(0)
    if IfVGGNN == "No":
        NNsuffix = ""
    elif IfVGGNN == "Yes":
        NNsuffix = "_VGG16NN"

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


