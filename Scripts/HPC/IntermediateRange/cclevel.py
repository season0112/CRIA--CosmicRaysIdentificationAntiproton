#!/usr/bin/env python
import numpy as np
import os
from root_numpy import array2root, root2array
import argparse
import binning
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TMultiGraph, TGraphErrors

data_beforecccut = root2array("$HPCINTERMEDIATEDIR/total/cclevel/B1042_pr.pl1.1800_7.6_all_Tree.root", "AntiprotonIntermediateEnergyTree")
data_aftercccut = data_beforecccut[np.where(data_beforecccut["ProtonCCMVABDT"]>0.9)[0]]
binning = np.array([2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11, 12, 13, 14.1, 15.3, 16.6, 18])

binningcenter  = np.array([
    3.13,  3.465,  3.83,  4.225,  4.655,  5.125,  5.635,  6.185,
    6.78,  7.425, 8.12,   8.87,  9.68,   10.55,  11.5,   12.5,   13.55,
    14.7,  15.95, 17.3])

ratio_beforecccut=np.array([])
ratio_aftercccut=np.array([])
error_beforecccut=np.array([])
error_aftercccut=np.array([])

for index in range(binning.shape[0]-1):
    positive_beforecccut = data_beforecccut[np.where( (data_beforecccut["Rigidity"]>binning[index]) & (data_beforecccut["Rigidity"]<binning[index+1]))[0]]
    negative_beforecccut = data_beforecccut[np.where( (data_beforecccut["Rigidity"]>-binning[index+1]) & (data_beforecccut["Rigidity"]<-binning[index]))[0]]
    p_count_before = positive_beforecccut.shape[0]
    n_count_before = negative_beforecccut.shape[0]
    ratio_tem_beforecccut = n_count_before/(p_count_before+n_count_before)
    error_tem_beforecccut = np.sqrt( p_count_before**2*n_count_before/(n_count_before+p_count_before)**4 + n_count_before**2*p_count_before/(n_count_before+p_count_before)**4 )
    ratio_beforecccut = np.append(ratio_beforecccut, ratio_tem_beforecccut)
    error_beforecccut = np.append(error_beforecccut, error_tem_beforecccut)

for index in range(binning.shape[0]-1):
    positive_aftercccut = data_aftercccut[np.where( (data_aftercccut["Rigidity"]>binning[index]) & (data_aftercccut["Rigidity"]<binning[index+1]))[0]]
    negative_aftercccut = data_aftercccut[np.where( (data_aftercccut["Rigidity"]>-binning[index+1]) & (data_aftercccut["Rigidity"]<-binning[index]))[0]]
    p_count_after = positive_aftercccut.shape[0]
    n_count_after = negative_aftercccut.shape[0]
    ratio_tem_aftercccut = n_count_after/(p_count_after+n_count_after)
    error_tem_aftercccut = np.sqrt( p_count_after**2*n_count_after/(n_count_after+p_count_after)**4 + n_count_after**2*p_count_after/(n_count_after+p_count_after)**4 )
    ratio_aftercccut = np.append(ratio_aftercccut,ratio_tem_aftercccut)
    error_aftercccut = np.append(error_aftercccut, error_tem_aftercccut)





ex = np.zeros(binningcenter.shape[0])
ccbefore = TGraphErrors (binningcenter.shape[0], binningcenter , ratio_beforecccut*1000, ex, error_beforecccut*1000 )
ccafter = TGraphErrors (binningcenter.shape[0], binningcenter , ratio_aftercccut*1000, ex, error_aftercccut*1000 )

c2 = TCanvas("c1","c1",1000,500)
#c2.SetLogy()
ccbefore.Draw("A*")
ccafter.Draw("same*")
ccbefore.SetMarkerStyle(53)
ccafter.SetMarkerStyle(54)
#ccbefore.SetMarkerSize(20)
#ccafter.SetMarkerSize(20)
ccbefore.SetMarkerColor(4)
ccafter.SetMarkerColor(2)
ccbefore.SetLineColor(4)
ccafter.SetLineColor(2)
ccbefore.GetXaxis().SetTitle("Rigidity (GV)")
ccbefore.GetYaxis().SetTitle("CCLevel")
ccbefore.SetTitle("CCLevel")
gPad.SetGrid(0,0)
#ccbefore.GetYaxis().SetMoreLogLabels();
leg =TLegend(.7,.7,.9,.9,)
leg.SetFillColor(0)
leg.AddEntry(ccbefore,"Before CCcut")
leg.AddEntry(ccafter,"After CCcut")
leg.Draw()
c2.Update()
c2.SaveAs("$HPCINTERMEDIATEDIR/total/cclevel/test2.pdf")

'''
plt.figure(figsize=(18,9))
plt.plot(binningcenter, ratio_beforecccut , 'bo-',lw=3, label='CClevel before cccut')
plt.plot(binningcenter, ratio_aftercccut , 'rs-',lw=3, label='CClevel after cccut')
plt.xlabel('Rigidity',fontsize=22)
plt.ylabel('CClevel',fontsize=22)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#yaxis.SetMoreLogLabels()
plt.yscale('log')
#yaxis.set_minor_formatter(ticker.LogFormatterMathtext(labelOnlyBase=False))
plt.legend( loc='best',fontsize=20)
plt.savefig('$HPCINTERMEDIATEDIR/total/cclevel/test.png')
'''





