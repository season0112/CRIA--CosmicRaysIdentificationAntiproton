#!/usr/bin/env python
import numpy as np
import math
import json
import collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import os


def main():
    CCProton_Merged = np.array([])
    for i in bins:
        if IfVGGP0 == "No":
            CCProton = np.load(datapath + "/ISS_anylsis/data/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_" + str(i)+ "_Pattern_" + str(pattern) + ".npy")
        elif IfVGGP0 == "Yes":
            CCProton = np.load(datapath + "/ISS_anylsis/data/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_" + str(i)+ "_Pattern_" + str(pattern) + "_VGG16NN.npy")

        if CCProton_Merged.shape[0] == 0:
            CCProton_Merged = CCProton
        else:
            CCProton_Merged = np.row_stack((CCProton_Merged, CCProton))
        if IfVGGP0 == "No":
            np.save( datapath + "/ISS_anylsis/data/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_" + "14.1_38.9" + "_Pattern_" + str(pattern) + ".npy", CCProton_Merged)
        elif IfVGGP0 == "Yes":
            np.save( datapath + "/ISS_anylsis/data/ChargeConfusedProtomTemplate_MC_B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_negative_" + "14.1_38.9" + "_Pattern_" + str(pattern) + "_VGG16NN.npy", CCProton_Merged) 


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pattern', help='which tracker pattern you choose')
    parser.add_argument('--ifVGGP0', type=str, help='For pattern 0, if turn to Neural Network.')
    arguments = parser.parse_args()

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

    datapath = os.getenv('HPCHIGHENERGYDATADIR')

    bins = ["14.1_15.3", "15.3_16.6", "16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9"]

    main()

