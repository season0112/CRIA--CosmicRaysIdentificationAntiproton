#!/usr/bin/env python
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
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors, TChain
import numpy.polynomial.polynomial as poly
from scipy import optimize
from root_numpy import root2array, tree2array
from root_numpy import testdata

def main():

    for Efficiency in Efficiency_all:

        for binmerge in Binmerge:

            cutvalue_all = []
            chain = TChain(lowtreename,'')

            for index in plotrange:

                if binmerge == 1:
                    chain.AddFile(lowworkpath + "/B1042_antipr.pl1.1800_7.6_all_Tree_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+1]) + ".root")
                elif binmerge == 2:
                    chain.AddFile(lowworkpath + "/B1042_antipr.pl1.1800_7.6_all_Tree_" + "{:g}".format(binning.published2016binnings[index]) + "_" + "{:g}".format(binning.published2016binnings[index+1]) + ".root")
                    chain.AddFile(lowworkpath + "/B1042_antipr.pl1.1800_7.6_all_Tree_" + "{:g}".format(binning.published2016binnings[index+1]) + "_" + "{:g}".format(binning.published2016binnings[index+2]) + ".root")

                print("index is " + str(index))
                print(chain.GetEntries()) ## Do not delete! important! otherwise will not load this chain.

                antiprotontree = chain.GetTree()
                antiprotonarray = tree2array(antiprotontree)

                antiprotonarray_large_0 = antiprotonarray['RichBeta'][np.where(antiprotonarray['RichBeta']>0)[0]]
                cutvalue = np.sort(antiprotonarray_large_0)[int(antiprotonarray_large_0.shape[0]*Efficiency)]
                cutvalue_all.append(cutvalue)

            with open(intermediateworkpath + "/RichBetaCutValue_eff_" + str(Efficiency) + "_" + str(binmerge)+".txt","w") as f:
                for item in cutvalue_all:
                    f.write("%s\n" % item)


if __name__ == '__main__':

    #### Argument
    parser = argparse.ArgumentParser()
    parser.add_argument('-C','--cluster', help='which cluster you choose')
    arguments = parser.parse_args()

    Efficiency_all = [0.8, 0.9] 
    Binmerge       = [1, 2]

    rigidity_start = 9 # correspondent to 2.97
    rigidity_end   = 29 # correspondent to 18.0
    plotrange      = range(rigidity_start,rigidity_end, int(Binmerge))

    if arguments.cluster == "JUAMS":
        intermediateworkpath = os.getenv('JUAMSINTERMEDIATEENERGYDIR')
    elif arguments.cluster == "HPC":
        intermediateworkpath = os.getenv('HPCINTERMEDIATEDIR') + '/total'
        intermediatreename   = 'AntiprotonIntermediateEnergyTree'
        lowworkpath = os.getenv('HPCLOWENERGYDIR') + '/totalall'
        lowtreename  = 'AntiprotonLowEnergyTree'

    main()
