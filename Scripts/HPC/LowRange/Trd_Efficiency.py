#!/usr/bin/env python
import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
from root_numpy import array2tree, array2root, tree2array
import argparse
import os
import binning

parser = argparse.ArgumentParser()
parser.add_argument('-C','--cluster', help='which cluster you choose')
parser.add_argument('-E','--efficiency', type=float, help='what antiproton signal efficiency you choose.')
arguments = parser.parse_args()

efficiency = arguments.efficiency

if arguments.cluster == "JUAMS":
    print("in progress...")
elif arguments.cluster == "HPC":
    lowpath = os.getenv('HPCLOWENERGYDIR')

lowpath = '/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v7.0'

PositiveCut_exceptTRD = 'TrdLogLikelihoodRatioProtonHeliumTracker < 0.1 && ProtonCCMVABDT > 0.9 && EcalBDT_EnergyD < -0.9 && TrdLogLikelihoodRatioElectronProtonTracker > -1.5' # Richbeta is also not included.

cutvalue_all = []
for i in range (2,16):
    print("Now:" + str(binning.published2016binnings[i]) + "_" + str(binning.published2016binnings[i+1]))
    positive = TFile(lowpath + "/totalall/" + "B1130_pass7_7.8_all_Tree_positive_May2015_" + str(binning.published2016binnings[i]) + "_" + str(binning.published2016binnings[i+1]) + ".root")
    tree = positive.Get("AntiprotonLowEnergyTree")
    array = tree2array(tree, selection = PositiveCut_exceptTRD)
    cutvalue = np.sort(array['TrdLogLikelihoodRatioElectronProtonTracker'])[int(array.shape[0]*(1-efficiency))]
    print("Cut value is " + str(cutvalue))
    cutvalue_all.append(cutvalue)

with open(lowpath + "/totalall/Trd_RichBeta_Cut/" + "TrdlikelihoodCutValue_eff_" + str(efficiency) +".txt","w") as f:
    for item in cutvalue_all:
        f.write("%s\n" % item)




