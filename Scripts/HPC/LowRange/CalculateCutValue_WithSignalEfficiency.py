#!/usr/bin/env python
import numpy as np
import os
import argparse
import binning
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors, TChain
from root_numpy import root2array, tree2array


def main():

    #### Loop in Trd Efficiency
    for efficiency in UsedEfficiencyAll:

        #### Loop in Binmerge
        for binmerge in BinMergeAll:

            print("Now:" + "EfficiencyType is " + str(EfficiencyType) + ", efficiency is " + str(efficiency) + ", " + "binmerge is " + str(binmerge))
            cutvalue_all = []
            RigidityRange        = range(rigidity_start, rigidity_end, int(binmerge))
            chain = TChain('AntiprotonLowEnergyTree','')

            #### Loop in Rigidity
            for index in RigidityRange:

                print("Now we are in:" + str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+1]))

                #### Load Root File
                if binmerge == 1:
                    chain.AddFile(lowworkpath + "/totalall" + "/" + ProtonDataSet + "_Tree_positive_" + str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+1]) + ".root")
                elif binmerge == 2:
                    chain.AddFile(lowworkpath + "/totalall" + "/" + ProtonDataSet + "_Tree_positive_" + str(binning.published2016binnings[index]) + "_" + str(binning.published2016binnings[index+1]) + ".root")
                    chain.AddFile(lowworkpath + "/totalall" + "/" + ProtonDataSet + "_Tree_positive_" + str(binning.published2016binnings[index+1]) + "_" + str(binning.published2016binnings[index+2]) + ".root")
                print(chain.GetEntries()) ## Do not delete! important! otherwise will not load this chain.

                #### Get Tree after Cuts  
                antiprotontree = chain.GetTree()
                antiprotonarray = tree2array(antiprotontree)

                if EfficiencyType == "TRD":
                    antiprotonarray = antiprotonarray['TrdLogLikelihoodRatioElectronProtonTracker'][np.where(antiprotonarray['TrdLogLikelihoodRatioElectronProtonTracker']>-1.5)[0]]
                elif EfficiencyType == "TOF":
                    antiprotonarray = 1.0/antiprotonarray['TofBeta'] - 1.0 / np.sqrt(antiprotonarray['Rigidity']**2 / ( 0.938**2 + antiprotonarray['Rigidity']**2 ))
                cutvalue = np.flipud(np.sort(antiprotonarray))[int(antiprotonarray.shape[0]*efficiency)]
                cutvalue_all.append(cutvalue)


            #### Save in txt File
            if EfficiencyType == "TRD":
                with open(lowworkpath + "/totalall/" + "TRDLogLikelihood_CutValue_eff_" + str("{:.2f}".format(efficiency)) + "_" + str(binmerge)+".txt","w") as f:
                    for item in cutvalue_all:
                        f.write("%s\n" % item)
            elif EfficiencyType == "TOF":
                with open(lowworkpath + "/totalall/" + "TOFBeta_CutValue_eff_" + str("{:.2f}".format(efficiency)) + "_" + str(binmerge)+".txt","w") as f:
                    for item in cutvalue_all:
                        f.write("%s\n" % item)


            chain.Reset()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--efficiencyType', help='which efficiency you calculate')
    arguments = parser.parse_args()

    if (arguments.efficiencyType):
        EfficiencyType = arguments.efficiencyType
    else:
        print("You need to provide all parameters!")
        os._exit(0)

    rigidity_start = 0  # correspondent to 1.0
    rigidity_end   = 16 # correspondent to 5.9

    lowworkpath = os.getenv('HPCLOWENERGYDIR')

    ProtonDataSet = 'B1042_pr.pl1.1800_7.6_all'
    #ProtonDataSet = 'B1130_pass7_7.8_all'

    #TrdEfficiencyAll = np.arange(0.7, 0.96, 0.02)
    #TofEfficiencyAll = np.arange(0.7, 0.96, 0.02)
    TrdEfficiencyAll = np.arange(0.6, 1.00, 0.01)
    TofEfficiencyAll = np.arange(0.6, 1.00, 0.01)

    if EfficiencyType == "TRD":
        UsedEfficiencyAll = TrdEfficiencyAll
    elif EfficiencyType == "TOF":
        UsedEfficiencyAll = TofEfficiencyAll

    BinMergeAll      = [1, 2]

    main()



