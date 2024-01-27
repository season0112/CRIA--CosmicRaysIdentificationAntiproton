#!/usr/bin/env python
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile
from root_numpy import fill_hist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import argparse
plt.switch_backend('agg')

fitresult = np.ones((16,30))
antiprotonratio = np.ones((16,30))
j=0

for cccutvalue in np.arange(0.2,1,0.05):
#for cccutvalue in np.arange(0.2,0.25,0.05):
    cccutvalue = "{:.2f}".format(cccutvalue)
    proton_number_TF = open("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_pass6/negative/proton_number_pass6_"+cccutvalue+".txt","r")
    f = open("/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_pass6/negative/results/fit_resultscccut"+cccutvalue+"CCN20TRDN12.txt","r") 
    proton_number = proton_number_TF.readlines()
    lines = f.readlines()
    for i in range(len(lines)):
        fitresult[j,i] = float(lines[i].strip('\n').split()[0:1][0].strip(','))
        antiprotonratio[j,i] = fitresult[j,i]/float(proton_number[i].strip("\n"))
    j=j+1
print(fitresult)
print(fitresult.shape)

NNbinnings = np.array([16.6, 18.0, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93, 108, 125, 147, 175, 211, 259, 330, 525])
NNpoint=np.array( np.zeros(NNbinnings.shape[0]-1) )
for i in range(NNpoint.shape[0]):
    NNpoint[i]=(NNbinnings[i]+NNbinnings[i+1])/2

fig, ax = plt.subplots(figsize=(22,18))
ax.set_yscale('log')
for i in range(16):
    plt.plot(NNpoint, antiprotonratio[i,:], '-o',lw=3, label='cccut' + "{:.2f}".format(np.arange(0.2,1,0.05)[i]))
#plt.plot(NNpoint, antiprotonratio[0,:], '-o',lw=3, label='cccut'+str(0.2) )
plt.xlabel('Rigidity',fontsize=20)
plt.ylabel('Antiproton to Proton Ratio',fontsize=20)
#plt.grid(True, which='both', axis='y')
plt.grid(True, which='both', axis='both')
my_y_ticks = np.array([0.00016,0.00018,0.0002,0.00022,0.00030])
plt.yticks(my_y_ticks)
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='minor',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.legend(loc='best',fontsize=20)
plt.vlines(16.6, min(antiprotonratio[0,:]), max(antiprotonratio[0,:]), colors = "c", linestyles = "dashed")
plt.savefig('/home/bo791269/v7.0_05.03.2019/VGG16_v7.0/VGG16_v7_mva_only_yescccut_yesecalcut_yestrd_v3.5_withmvaonly_ele_posi_same_2electronMCused_protonfromTF_trdccbinTuned/pattern_0/ISS_anylsis/data/templatefit_pass6/negative/results/result_cccut.png')

