#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
from ROOT import TFile, TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import array2tree, array2root, tree2array, root2array
import binning


#### Parse Argument
parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='which cluster you choose')
parser.add_argument('--rigidityrange', help='which rigidityrange you choose, high? intermediate?')
parser.add_argument('--data', help='which data you choose to split')
parser.add_argument('--datamode', help='the dataset is MC or ISS data')
parser.add_argument('-S','--Rigidity_start', type=int, help='Rigidity_start')
parser.add_argument('-E','--Rigidity_end', type=int, help='Rigidity_end')
arguments = parser.parse_args()


#### Load Path
if arguments.rigidityrange == "intermediaterange": 
    if arguments.cluster=="JUAMS":
        workpath = os.getenv('JUAMSINTERMEDIATEENERGYDIR')
    elif arguments.cluster=="HPC":
        workpath = os.getenv('HPCINTERMEDIATEDIR')
    Treename = "AntiprotonIntermediateEnergyTree"
    binningused = binning.published2016binnings
elif arguments.rigidityrange == "highrange":
    if arguments.cluster=="JUAMS":
        print("work in progress......")
        os._exit(0)
    elif arguments.cluster=="HPC":
        workpath = os.getenv('HPCHIGHENERGYDATADIR')
    Treename = "ExampleAnalysisTree"
    binningused = binning.published2016binnings
elif arguments.rigidityrange == "lowrange":
    if arguments.cluster=="JUAMS":
        print("work in progress......")
        os._exit(0)
    elif arguments.cluster=="HPC":
        workpath = os.getenv('HPCLOWENERGYDIR')
    Treename = "AntiprotonLowEnergyTree"
    binningused = np.array([ 0.8, 1, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11, 12, 13, 14.1, 15.3, 16.6, 18,]) ## Naive fix to add 0.8-1.0 GV

else:
    print("you must choose a rigidity range correctly!")
    os._exit(0)

rootfile = arguments.data


#### Load root file and cut ont unixtime for time split
negative_array = root2array( workpath + "/" + str(rootfile)+".root", Treename)
negative_May2015 = negative_array[np.where(negative_array['TimeStamp']<1432598400)[0]] #May 26, 2015 12:00:00 AM
negative_Nov2017 = negative_array[np.where(negative_array['TimeStamp']<1510484400)[0]] #November 12, 2017 11:00:00 AM


#### Split in Rigidity
for rigidity_edge in range(arguments.Rigidity_start,arguments.Rigidity_end):  # for intermediate range:9-29: 2.97-18.0. for high range:28-57: 16.6-450.

    rigidity = "{:g}".format(binningused[rigidity_edge]) + "_" + "{:g}".format(binningused[rigidity_edge+1])
    
    # Split for "MC" or "full time ISS negative data"
    negative_array_tem = negative_array[np.where(negative_array['Rigidity']>-binningused[rigidity_edge+1])[0]]
    negative_array_tem = negative_array_tem[np.where(negative_array_tem['Rigidity']<-binningused[rigidity_edge])[0]]
    array2root(negative_array_tem, workpath + '/' + arguments.data + '_'+str(rigidity)+'.root', Treename, 'recreate')

    # Split for PRL time ISS negative data
    if arguments.datamode == "MC":
        print("the dataset is MC. No Timestamp process will be  done.")
    elif arguments.datamode == "ISSdata":
        negative_array_tem2 = negative_May2015[np.where(negative_May2015['Rigidity']>-binningused[rigidity_edge+1])[0]]
        negative_array_tem2 = negative_array_tem2[np.where(negative_array_tem2['Rigidity']<-binningused[rigidity_edge])[0]]
        array2root(negative_array_tem2, workpath + '/' + arguments.data + '_'+ 'May2015_' + str(rigidity)+'.root', Treename, 'recreate')
    else:
        os._exit(0)
    print(rigidity+"finished.")
    
    # Split for PhysRep time ISS negative data
    if arguments.datamode == "MC":
        print("the dataset is MC. No Timestamp process will be  done.")
    elif arguments.datamode == "ISSdata":
        negative_array_tem2 = negative_Nov2017[np.where(negative_Nov2017['Rigidity']>-binningused[rigidity_edge+1])[0]]
        negative_array_tem2 = negative_array_tem2[np.where(negative_array_tem2['Rigidity']<-binningused[rigidity_edge])[0]]
        array2root(negative_array_tem2, workpath + '/' + arguments.data + '_'+ 'Nov2017_' + str(rigidity)+'.root', Treename, 'recreate')
    else:
        os._exit(0)
    print(rigidity+"finished.")


#### Special step for high range: other rigidity bins split.
if arguments.rigidityrange == "highrange":
    print("Since the rigidity range is High Rigidity Rarnge, the special dealing for 450 and 525 version differences starts. ")

    left  = [259, 330, 211, 250] 
    right = [330, 525, 250, 330]

    for (i, j) in zip(left, right):
        # Split for full time ISS negative data
        negative_array_tem = negative_array[np.where(negative_array['Rigidity']>-j)[0]]
        negative_array_tem = negative_array_tem[np.where(negative_array_tem['Rigidity']<-i)[0]]
        array2root(negative_array_tem, workpath + '/' + arguments.data + '_' + str(i) + '_' + str(j) +'.root', Treename, 'recreate')

        # Split for PRL time ISS negative data
        if arguments.datamode == "MC":
            print("the dataset is MC. No Timestamp process will be  done.")
        elif arguments.datamode == "ISSdata":
            negative_array_tem2 = negative_May2015[np.where(negative_May2015['Rigidity']>-j)[0]]
            negative_array_tem2 = negative_array_tem2[np.where(negative_array_tem2['Rigidity']<-i)[0]]
            array2root(negative_array_tem2, workpath + '/' + arguments.data + '_'+ 'May2015_' + str(i) + '_' + str(j) + '.root', Treename, 'recreate')
        else:
            os._exit(0)

        # Split for PhysRep time ISS negative data
        if arguments.datamode == "MC":
            print("the dataset is MC. No Timestamp process will be  done.")
        elif arguments.datamode == "ISSdata":
            negative_array_tem2 = negative_Nov2017[np.where(negative_Nov2017['Rigidity']>-j)[0]]
            negative_array_tem2 = negative_array_tem2[np.where(negative_array_tem2['Rigidity']<-i)[0]]
            array2root(negative_array_tem2, workpath + '/' + arguments.data + '_'+ 'Nov2017_' + str(i) + '_' + str(j) + '.root', Treename, 'recreate')
        else:
            os._exit(0)

        print( str(i) + "_" + str(j) + "finished.")


