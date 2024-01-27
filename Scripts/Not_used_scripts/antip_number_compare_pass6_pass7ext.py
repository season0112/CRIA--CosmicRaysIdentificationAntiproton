#!/usr/bin/env python
from __future__ import division
import numpy as np
import math
import json
import collections
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
from root_numpy import array2tree, array2root, fill_hist
from ROOT import TFile, TH1D, TCanvas, gPad, gStyle, TPad, TLegend


#### NN result pass6
NNbinnings_pass6 = np.array([38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93, 108, 125, 147, 175, 211, 259, 330, 525])
antiproton_number_pass6 = np.array([ 732,613,575,512,435,361,332,283,262,211 ,369, 280, 200, 168, 175, 104, 101, 87, 70]) #### bins: 330-525:TRD:11 CC:9, cccut=0.7 259-330:TRD:11 CC:9.cccut=0.7 others:TRD:12 CC:20 cccut=0.2
proton_number_pass6 = np.array([4211803,3668867,3184018,2835609,2447051,2161489,1898563,1695342,1472061,1298441,2146489,1714722,1289892,1087198,865939,672067,520330,367076,384926])  #### raw proton number bins: 330-525 and 259-330:cccut=0.7 others:cccut=0.2

#### NN result pass7
NNbinnings_pass7ext = np.array([16.6, 18, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93, 108, 125, 147, 175, 211, 259, 330, 525])
antiproton_number_pass7ext = np.array([3655, 3216, 2847, 2502, 2309, 1978, 1717, 1520, 1261, 1120, 1051, 841, 686, 666, 579, 496, 422, 406, 347, 302, 251, 431, 320, 231, 204, 188, 114, 111, 84, 79]) #### bins: 330-525:TRD:11 CC:9, cccut=0.7 259-330:TRD:11 CC:9.cccut=0.7 others:TRD:12 CC:20 cccut=0.2
proton_number_pass7ext = np.array([21979364, 19263347, 16689065, 14389772, 13011251, 11052934, 9415821, 8369767, 7102510, 6265889, 5489010, 4701078, 4093107, 3546336, 3153585, 2717903, 2396250, 2102873, 1874303, 1624258, 1431628, 2361218, 1855990, 1408431, 1172795, 965490, 745050, 584025, 423701, 443781]) #### raw proton number bins: 330-525 and 259-330:cccut=0.7 others:cccut:0.2

#### Published2016
published2016binnings = np.array([1, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11, 12, 13, 14.1, 15.3, 16.6, 18, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1, 38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9, 80.5, 93, 108, 125, 147, 175, 211, 259, 450])
antiproton_number_published = np.array([15816, 15049, 14426, 13511, 12943, 11723, 10411, 9508, 7876, 7212, 6127, 2697, 2353, 1962, 1772, 1528, 1300, 1143, 1002, 916, 841, 1270, 980, 733, 573, 233, 83, 72, 66]) ## 4.5 years data
proton_number_published = np.array([79080000., 75623115., 72492462., 66886138., 65040201., 58034653., 51539603., 47069306., 41020833., 36060000., 30944444., 14046875., 12515957., 10548387.,  9277486.,  8000000.,  6701030.,  6211956., 5661016.,  4978260.,  4335051.,  7055555.,  5600000.,  4261627., 3237288.,  1339080.,   446236.,   327272.,   461538.]) ## 4.5 years data
antiproton_number_expect = antiproton_number_published/4.5*7.5
proton_number_expect = proton_number_published/4.5*7.5

NNpoint=np.array( np.zeros(NNbinnings_pass6.shape[0]-1) )
int=np.array( np.zeros(NNbinnings_pass6.shape[0]-1) )
for i in range(NNpoint.shape[0]):
    NNpoint[i]=(NNbinnings_pass6[i]+NNbinnings_pass6[i+1])/2
NNbinningswidth = np.ones(NNbinnings_pass6.shape[0]-1)
for i in range(NNbinnings_pass6.shape[0]-1):
    NNbinningswidth[i] = NNbinnings_pass6[i+1] - NNbinnings_pass6[i]

published2016point=np.array( np.zeros(published2016binnings.shape[0]-1) )
for i in range(published2016point.shape[0]):
    published2016point[i]=(published2016binnings[i]+published2016binnings[i+1])/2


TH_antiproton_pass6 = TH1D("","", 19 ,NNbinnings_pass6)
TH_antiproton_pass7ext = TH1D("","", 19 ,NNbinnings_pass6)
TH_proton_pass6 = TH1D("","", 19 ,NNbinnings_pass6)
TH_proton_pass7ext = TH1D("","", 19 ,NNbinnings_pass6)
TH_antiproton_expect = TH1D("","", 18 ,published2016binnings[-19:])
TH_proton_expect = TH1D("","", 18 ,published2016binnings[-19:])


for i in range(1,20):
    TH_antiproton_pass6.SetBinContent(i,antiproton_number_pass6[i-1])
    TH_antiproton_pass7ext.SetBinContent(i,antiproton_number_pass7ext[i-1+11])
    TH_proton_pass6.SetBinContent(i,proton_number_pass6[i-1])
    TH_proton_pass7ext.SetBinContent(i,proton_number_pass7ext[i-1+11])
for i in range(1,19):
    TH_antiproton_expect.SetBinContent(i,antiproton_number_expect[i-1+11])
    TH_proton_expect.SetBinContent(i,proton_number_expect[i-1+11])

c1 = TCanvas()
gStyle.SetOptStat(0)
TH_antiproton_pass6.Draw("")
TH_antiproton_pass7ext.Draw("same")
xaxis = TH_antiproton_pass6.GetXaxis()
xaxis.SetTitle("|Rigidity| (GV)")
xaxis.SetMoreLogLabels()
xaxis.SetLabelSize(0.05)
xaxis.SetTitleSize(0.05)
#xaxis.SetTickSize(0.1)
yaxis = TH_antiproton_pass6.GetYaxis()
#yaxis.SetTitle("MC / ISS")
yaxis.SetLabelSize(0.05)
yaxis.SetTitleSize(0.06)
yaxis.SetMoreLogLabels()
gPad.SetLogx()
gPad.SetLogy()
TH_antiproton_pass6.SetLineColor(632)
TH_antiproton_pass7ext.SetLineColor(1)
leg = TLegend(0.7,0.7,0.9,0.85)
leg.AddEntry(TH_antiproton_pass6,"pass6","lp")
leg.AddEntry(TH_antiproton_pass7ext,"pass7ext","lp")
leg.Draw()
c1.SaveAs("/home/bo791269/v7.0_05.03.2019/v7.0_results/v7.0_yesecalcut_yescccut_v3.5/antiproton_to_proton_ratio_plot/antiproton.pdf")


c2 = TCanvas()
gStyle.SetOptStat(0)
TH_proton_pass6.Draw("")
TH_proton_pass7ext.Draw("same")
xaxis = TH_proton_pass6.GetXaxis()
xaxis.SetTitle("|Rigidity| (GV)")
xaxis.SetMoreLogLabels()
xaxis.SetLabelSize(0.05)
xaxis.SetTitleSize(0.05)
#xaxis.SetTickSize(0.1)
yaxis = TH_proton_pass6.GetYaxis()
#yaxis.SetTitle("MC / ISS")
yaxis.SetLabelSize(0.05)
yaxis.SetTitleSize(0.06)
yaxis.SetMoreLogLabels()
gPad.SetLogx()
gPad.SetLogy()
TH_proton_pass6.SetLineColor(632)
TH_proton_pass7ext.SetLineColor(1)
leg = TLegend(0.7,0.7,0.9,0.85)
leg.AddEntry(TH_proton_pass6,"pass6","lp")
leg.AddEntry(TH_proton_pass7ext,"pass7ext","lp")
leg.Draw()
c2.SaveAs("/home/bo791269/v7.0_05.03.2019/v7.0_results/v7.0_yesecalcut_yescccut_v3.5/antiproton_to_proton_ratio_plot/proton.pdf")


c3 = TCanvas()
gStyle.SetOptStat(0)
TH_antiproton_expect.Draw("")
TH_antiproton_pass7ext.Draw("same")
xaxis = TH_antiproton_expect.GetXaxis()
xaxis.SetTitle("|Rigidity| (GV)")
xaxis.SetMoreLogLabels()
xaxis.SetLabelSize(0.05)
xaxis.SetTitleSize(0.05)
#xaxis.SetTickSize(0.1)
yaxis = TH_antiproton_expect.GetYaxis()
#yaxis.SetTitle("MC / ISS")
yaxis.SetLabelSize(0.05)
yaxis.SetTitleSize(0.06)
yaxis.SetMoreLogLabels()
gPad.SetLogx()
gPad.SetLogy()
TH_antiproton_expect.SetLineColor(632)
TH_antiproton_pass7ext.SetLineColor(1)
leg = TLegend(0.7,0.7,0.9,0.85)
leg.AddEntry(TH_antiproton_expect,"Expected from 2016paper","lp")
leg.AddEntry(TH_antiproton_pass7ext,"pass7ext","lp")
leg.Draw()
c3.SaveAs("/home/bo791269/v7.0_05.03.2019/v7.0_results/v7.0_yesecalcut_yescccut_v3.5/antiproton_to_proton_ratio_plot/antiproton_comparewithexpect.pdf")



c4 = TCanvas()
gStyle.SetOptStat(0)
TH_proton_expect.Draw("")
TH_proton_pass7ext.Draw("same")
xaxis = TH_proton_expect.GetXaxis()
xaxis.SetTitle("|Rigidity| (GV)")
xaxis.SetMoreLogLabels()
xaxis.SetLabelSize(0.05)
xaxis.SetTitleSize(0.05)
#xaxis.SetTickSize(0.1)
yaxis = TH_proton_expect.GetYaxis()
#yaxis.SetTitle("MC / ISS")
yaxis.SetLabelSize(0.05)
yaxis.SetTitleSize(0.06)
yaxis.SetMoreLogLabels()
gPad.SetLogx()
gPad.SetLogy()
TH_proton_expect.SetLineColor(632)
TH_proton_pass7ext.SetLineColor(1)
leg = TLegend(0.7,0.7,0.9,0.85)
leg.AddEntry(TH_proton_expect,"Expected from 2016paper","lp")
leg.AddEntry(TH_proton_pass7ext,"pass7ext","lp")
leg.Draw()
c4.SaveAs("/home/bo791269/v7.0_05.03.2019/v7.0_results/v7.0_yesecalcut_yescccut_v3.5/antiproton_to_proton_ratio_plot/proton_comparewithexpect.pdf")











