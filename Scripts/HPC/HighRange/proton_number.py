#!/usr/bin/env python 
from __future__ import division
import numpy as np
import math
import json
import collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import os
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, gPad, TGraphErrors
import ROOT
plt.switch_backend('agg')


def ListToVector(rawlist):
    vector = ROOT.vector('double')(len(rawlist))
    vector.clear()
    for i in range(len(rawlist)):
        vector.insert(vector.begin()+i, rawlist[i])
    return vector


def main():

    print("IfVGGP0 is " + str(IfVGGP0))

    CCcutValue_All = np.arange(0.0, 1.0, 0.05)
    ProtonRawNumber_All = [[] for j in range(CCcutValue_All.shape[0])]

    if IfVGGP0 == "No":
        f_ProtonRawNumber = TFile.Open(workpath + "/templatefit/positive/ProtonRawNumber/ProtonRawNumber_Pattern_" + pattern + "_" + ISSversion + "_" + binning_version + ".root", "RECREATE")
    elif IfVGGP0 == "Yes":
        f_ProtonRawNumber = TFile.Open(workpath + "/templatefit/positive/ProtonRawNumber/ProtonRawNumber_Pattern_" + pattern + "_VGG16NN" + "_" + ISSversion + "_" + binning_version + ".root", "RECREATE")


    for index_R, i in enumerate (bins):
        print("Now is " + str(i))
        if IfVGGP0 == "No":
            data = np.load(workpath + "/ISS_positive" + suffix + str(i) + "_Pattern_" + pattern + ".npy")
        elif IfVGGP0 == "Yes":
            data = np.load(workpath + "/ISS_positive" + suffix + str(i) + "_Pattern_" + pattern + "_VGG16NN.npy")
        else:
            data = np.array([])
            print("Please check!")

        for index_CCcut, cccutvalue in enumerate (CCcutValue_All):
            cccutvalue = "{:.2f}".format(cccutvalue)  # keep 2 digits.
            number = data[np.where(data[:,0]>float(cccutvalue))[0],0].shape[0]
            ProtonRawNumber_All[index_CCcut].append(number)

    for index_CCcut, cccutvalue in enumerate (CCcutValue_All):
        cccutvalue = "{:.2f}".format(cccutvalue)  # keep 2 digits.
        f_ProtonRawNumber.WriteObject(ListToVector(ProtonRawNumber_All[index_CCcut]), "ProtonRawNumber_CCcut_"+cccutvalue)
    f_ProtonRawNumber.Close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--issversion', help='ISS data version, pass7.8 or publisehd2016')
    parser.add_argument('--binning', choices=['450version','525version'], help='publisehd(up to 450 GV) or new one (up to 525GV)')
    parser.add_argument('--pattern', help='Which tracker patterns you choose')
    parser.add_argument('--ifVGGP0', type=str, help='For pattern 0, if turn to Neural Network.')
    arguments = parser.parse_args()

    if (arguments.issversion):
        ISSversion = arguments.issversion
    else:
        print("You need to choose a ISS data cersion !")
        os._exit(0)

    if (arguments.binning):
        binning_version = arguments.binning
    else:
        print("You need to choose a binning version !")
        os._exit(0)

    if (arguments.pattern):
        pattern = arguments.pattern
    else:
        print("You need to choose a tracker pattern !")
        os._exit(0)

    if (arguments.ifVGGP0):
        IfVGGP0 = arguments.ifVGGP0
    else:
        print("You need to choose TMVA or Neural Network")
        os._exit(0)

    workpath = str(os.getenv('HPCHIGHENERGYDATADIR')) + str('/ISS_anylsis/data')

    if binning_version == "450version":
        bins = ["14.1_15.3", "15.3_16.6", "16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175", "175_211", "211_259", "259_450"]
    elif binning_version == "525version":
        bins = ["14.1_15.3", "15.3_16.6", "16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175", "175_211", "211_250", "250_330", "330_525"] # "211_259", "259_330"

    if ISSversion == "pass7.8":
        suffix = "_"
    elif ISSversion == "PhyRep2021":
        suffix = "_Nov2017_"
    elif ISSversion == "published2016":
        suffix = "_May2015_"

    main()




