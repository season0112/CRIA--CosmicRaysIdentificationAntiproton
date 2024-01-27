import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist

correct_proton = np.load('positive_250_330_pattern_0.npy')
correct_antiproton = np.load('negative_250_330_pattern_0.npy')

plt.figure(figsize=(18,18))
plt.hist(correct_proton[:,9],80,  alpha=0.5, range=(-10,10), density=True,  label='correct_proton',facecolor='blue',edgecolor='black')
plt.hist(correct_antiproton[:,9],80,  alpha=0.5, range=(-10,10), density=True,  label='correct_antiproton',facecolor='red',edgecolor='green')
plt.xlabel('L24L58RigidityMatching',fontsize=30)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.legend(loc='upper right',fontsize=30)
plt.savefig('9'+'.png')

plt.figure(figsize=(18,18))
plt.hist(correct_proton[:,8],80,  alpha=0.5, range=(-2.5,2.5), density=True,  label='correct_proton',facecolor='blue',edgecolor='black')
plt.hist(correct_antiproton[:,8],80,  alpha=0.5, range=(-2.5,2.5), density=True,  label='correct_antiproton',facecolor='red',edgecolor='green')
plt.xlabel('L1L9RigidityMatching',fontsize=30)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.legend(loc='upper right',fontsize=30)
plt.savefig('8'+'.png')

plt.figure(figsize=(18,18))
plt.hist(correct_proton[:,7],80,  alpha=0.5, range=(-2.5,2.5), density=True,  label='correct_proton',facecolor='blue',edgecolor='black')
plt.hist(correct_antiproton[:,7],80,  alpha=0.5, range=(-2.5,2.5), density=True,  label='correct_antiproton',facecolor='red',edgecolor='green')
plt.xlabel('InnerMaxSpanRigidityMatching',fontsize=30)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.legend(loc='upper right',fontsize=30)
plt.savefig('7'+'.png')


