#!/usr/bin/env python
from __future__ import division
from ROOT import TFile
from root_numpy import root2array, tree2array
import numpy as np
import sys
import os
import heapq
import time
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--rootfilenumber', help='rootfilenumber.')
parser.add_argument('--jobnumber', help='jobnumber.')
arguments = parser.parse_args()

begintime = time.time()

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path, exist_ok = True)
#########################################################################################

for treelistname in range(int(arguments.jobnumber)):
    treelistname = '{:05d}'.format(treelistname)
    rootfile = arguments.rootfilenumber
    result = None
    with open ('lists/'+'x'+treelistname,  "rb") as f:
        lines=f.readlines()
    for treecount in range(0,int(rootfile)):
        filename=os.path.basename(lines[treecount].strip()).decode()
        myfile = TFile(filename)
        intree = myfile.Get('AntiprotonIntermediateEnergyTree')
        array = tree2array(intree)
        if result is None:
            result = array
        else:
            result = np.append(result,array)


sign_Rigidity = result['Rigidity']/np.abs(result['Rigidity'])

plt.figure(figsize=(18,9))
plt.hist2d(result['TrdLogLikelihoodRatioElectronProtonTracker']*sign_Rigidity, result['RichBeta'],bins=200, range=[[-2.5,2.5],[0.95,1.02]] )
plt.savefig('test.png')

np.save( "test.npy", result )

#########################################################################################
endtime = time.time()
print ((endtime - begintime)/60)

