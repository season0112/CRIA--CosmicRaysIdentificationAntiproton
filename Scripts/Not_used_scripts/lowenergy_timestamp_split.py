#!/usr/bin/env python
from root_numpy import root2array, tree2array, array2root
from ROOT import TFile
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', help='current path')
parser.add_argument('--rigiditystart', type = int, help='which rigidity')
parser.add_argument('--rigidityend', type = int, help='which rigidity')
arguments = parser.parse_args()

rigiditystart = arguments.rigiditystart
rigidityend = arguments.rigidityend

os.chdir(arguments.path)

binning = np.array([ 0.8, 1, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11, 12, 13, 14.1, 15.3, 16.6, 18,])
PRL_splitdata = 1432598400 #May 26, 2015 12:00:00 AM
PhyRep_splitdata = 1510484400 #November 12, 2017 11:00:00 AM

negative = root2array("B1130_pass7_7.8_all_Tree_negative.root","AntiprotonLowEnergyTree")
antiproton = root2array("B1042_antipr.pl1.1800_7.6_all_Tree.root","AntiprotonLowEnergyTree")
electron = root2array("B1091_el.pl1.0_25200_7.6_all_Tree.root","AntiprotonLowEnergyTree")
#deuteron = root2array("B1128_d.pl1ph.021000_7.6_all_Tree.root","AntiprotonLowEnergyTree")
proton = root2array("B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative.root","AntiprotonLowEnergyTree")


# Rigidity split for MC and full time ISS negative data
for i in range(rigiditystart,rigidityend):
    print ("Full time R split begin:"+ str(binning[i]) + " to " + str(binning[i+1]))
    # Positive R
    deuteron_tem = deuteron[np.where((binning[i]<deuteron['Rigidity']) & (deuteron['Rigidity']<binning[i+1]))[0]]
    array2root(deuteron_tem, 'B1128_d.pl1ph.021000_7.6_all_Tree'+'_'+str(binning[i])+'_'+str(binning[i+1])+'.root', 'AntiprotonLowEnergyTree')
    # Negative R
    negative_tem = negative[np.where((-binning[i+1]<negative['Rigidity']) & (negative['Rigidity']<-binning[i]))[0]]
    antiproton_tem = antiproton[np.where((-binning[i+1]<antiproton['Rigidity']) & (antiproton['Rigidity']<-binning[i]))[0]]
    #electron_tem = electron[np.where((-binning[i+1]<electron['Rigidity']) & (electron['Rigidity']<-binning[i]))[0]]
    proton_tem = proton[np.where((-binning[i+1]<proton['Rigidity']) & (proton['Rigidity']<-binning[i]))[0]]
    array2root(negative_tem, 'B1130_pass7_7.8_all_Tree_negative'+'_'+str(binning[i])+'_'+str(binning[i+1])+'.root', 'AntiprotonLowEnergyTree')
    array2root(antiproton_tem, 'B1042_antipr.pl1.1800_7.6_all_Tree'+'_'+str(binning[i])+'_'+str(binning[i+1])+'.root', 'AntiprotonLowEnergyTree')
    #array2root(electron_tem, 'B1091_el.pl1.0_25200_7.6_all_Tree'+'_'+str(binning[i])+'_'+str(binning[i+1])+'.root', 'AntiprotonLowEnergyTree')
    array2root(proton_tem, 'B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative'+'_'+str(binning[i])+'_'+str(binning[i+1])+'.root', 'AntiprotonLowEnergyTree') 
print ("Full time R split finished.")

# Rigidity split for PRL time ISS negative data
negative_May2015 = negative[np.where(negative['TimeStamp']<PRL_splitdata)[0]]
array2root(negative_May2015, 'B1130_pass7_7.8_all_Tree_negative_May2015.root', 'AntiprotonLowEnergyTree')
for i in range(rigiditystart,rigidityend):
    print ("PRL paper time R split begin:"+ str(binning[i]) + " to " + str(binning[i+1]))
    negative_tem2 = negative_May2015[np.where((-binning[i+1]<negative_May2015['Rigidity']) & (negative_May2015['Rigidity']<-binning[i]))[0]]
    array2root(negative_tem2, 'B1130_pass7_7.8_all_Tree_negative_May2015'+'_'+str(binning[i])+'_'+str(binning[i+1])+'.root', 'AntiprotonLowEnergyTree')

# Rigidity split for PhysRep time ISS negative data
negative_Nov2017 = negative[np.where(negative['TimeStamp']<PhyRep_splitdata)[0]]
array2root(negative_Nov2017, 'B1130_pass7_7.8_all_Tree_negative_Nov2017.root', 'AntiprotonLowEnergyTree')
for i in range(rigiditystart,rigidityend):
    print ("2021 PhyeReport time R split begin:"+ str(binning[i]) + " to " + str(binning[i+1]))
    negative_tem2 = negative_Nov2017[np.where((-binning[i+1]<negative_Nov2017['Rigidity']) & (negative_Nov2017['Rigidity']<-binning[i]))[0]]
    array2root(negative_tem2, 'B1130_pass7_7.8_all_Tree_negative_Nov2017'+'_'+str(binning[i])+'_'+str(binning[i+1])+'.root', 'AntiprotonLowEnergyTree')


