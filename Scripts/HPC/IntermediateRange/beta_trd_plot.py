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
import matplotlib as mpl

proton = np.load("B1042_pr.pl1.flux.l1o9.2016000_7.6_all.npy")
electron = np.load("B1091_el.pl1.0_25200_7.6_all.npy")

proton_sign_Rigidity = proton['Rigidity']/np.abs(proton['Rigidity'])
electron_sign_Rigidity = electron['Rigidity']/np.abs(electron['Rigidity'])


plt.figure(figsize=(18,9))
#plt.hist2d(proton['TrdLogLikelihoodRatioElectronProtonTracker']*proton_sign_Rigidity, proton['RichBeta'],bins=200, range=[[-2.5,2.5],[0.95,1.02]],norm=mpl.colors.LogNorm())
plt.hist2d(electron['TrdLogLikelihoodRatioElectronProtonTracker']*electron_sign_Rigidity, electron['RichBeta'],bins=200, range=[[-2.5,2.5],[0.95,1.02]],norm=mpl.colors.LogNorm() )
plt.colorbar()
plt.savefig('beta_trd_plot.png')








