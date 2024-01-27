#!/usr/bin/env python
import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors, TChain
from root_numpy import array2tree, array2root, tree2array
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-C','--cluster', help='which cluster you choose')
parser.add_argument('-E','--efficiency', type=float, help='what antiproton signal efficiency you choose.')
parser.add_argument('-B','--binmerge', type=str, help='how many bins you used for each fit.')
arguments = parser.parse_args()

binmerge = int(arguments.binmerge)
efficiency = arguments.efficiency

if arguments.cluster == "JUAMS":
    print("in progress...")
elif arguments.cluster == "HPC":
    lowpath = os.getenv('HPCLOWENERGYDIR')

PositiveCut_RichBetaall = 'TrdLogLikelihoodRatioProtonHeliumTracker < 0.1 && ProtonCCMVABDT > 0.9 && EcalBDT_EnergyD < -0.9 && RichBeta>0' # TRDlikelihood is also not included.
PositiveCut_RichBeta = 'RichBeta>0'


lowbinning = np.array([0.8, 1, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11, 12, 13, 14.1, 15.3, 16.6, 18, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93, 108, 125, 147, 175, 211, 259, 450])

if arguments.binmerge == "1":
#    rigidity_start = 3 # correspondent to 1.33
#    rigidity_end = 17 # correspondent to 5.90
#    plotrange = range(rigidity_start,rigidity_end, int(binmerge))
    rigidity_start = 0 # correspondent to 0.8
    rigidity_end = 30 # correspondent to 18.0
    plotrange = range(rigidity_start,rigidity_end, int(binmerge))
if arguments.binmerge == "2":
    rigidity_start = 4 # correspondent to 1.51
    rigidity_end = 16 # correspondent to 5.37
    plotrange = range(rigidity_start,rigidity_end, int(binmerge))

cutvalue_all = []

chain = TChain('AntiprotonLowEnergyTree','')

for i in plotrange:
    if binmerge == 1:
        print("Now:" + str(lowbinning[i]) + "_" + str(lowbinning[i+1]))
        chain.AddFile( lowpath + "/totalall/" + "B1130_pass7_7.8_all_Tree_positive_May2015_" + "{:g}".format(lowbinning[i]) + "_" + "{:g}".format(lowbinning[i+1])  + ".root" )
    elif binmerge == 2:
        print("Now:" + str(lowbinning[i]) + "_" + str(lowbinning[i+2]))
        chain.AddFile( lowpath + "/totalall/" + "B1130_pass7_7.8_all_Tree_positive_May2015_" + "{:g}".format(lowbinning[i]) + "_" + "{:g}".format(lowbinning[i+1]) + ".root" )
        chain.AddFile( lowpath + "/totalall/" + "B1130_pass7_7.8_all_Tree_positive_May2015_" + "{:g}".format(lowbinning[i+1]) + "_" + "{:g}".format(lowbinning[i+2]) + ".root" )
    print(chain.GetEntries()) ## Do not delete! important! otherwise will not load this chain.
    tree = chain.GetTree()
    #array = tree2array(tree)
    #array = tree2array(tree, selection = PositiveCut_RichBetaall)
    array = tree2array(tree, selection = PositiveCut_RichBeta)
    cutvalue = np.sort(array['RichBeta'])[int(array.shape[0]*efficiency)]
    print("Cut value is " + str(cutvalue))
    cutvalue_all.append(cutvalue)

with open(lowpath + "/totalall/Trd_RichBeta_Cut/" + "RichBetaCutValue_eff_" + str(efficiency) + "_binmerge_"+ str(binmerge) + ".txt","w") as f:
    for item in cutvalue_all:
        f.write("%s\n" % item)




