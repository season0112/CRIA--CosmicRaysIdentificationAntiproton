from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TMultiGraph
from root_numpy import fill_hist

bin = [80.5, 93.0, 108.0, 125.0, 147.0, 175.0, 211.0, 250.0, 330.0, 525.0]
rigidity2 = np.array([])
for i in range(len(bin)-1):
    rigidity2 = np.append(rigidity2, np.array((bin[i]+bin[i+1])/2))
rigidity = np.append( np.array([33.6, 39.0, 48.65, 56.25, 65.0, 75.1]) , rigidity2 )
rigidity.astype(np.float)

ccnumber = np.array([1490, 1990, 2544, 3085, 3676, 4736, 5843, 7984, 9951, 13566, 18153, 23759, 25410, 47545, 86928]).astype(np.float)
correctnumber = np.array([92187167, 72002335, 48641517, 37497075, 29533378, 28661684, 17866111, 14404357, 10865174, 9209412, 7380493, 5768877, 3810445, 4223271, 3732153]).astype(np.float)
ccratio = ccnumber/(correctnumber+ccnumber)

ISSpositive = np.array([16165878, 12455567, 8304814, 6352927, 4953535, 3823826, 2961983, 2367438, 1780154, 1500592, 1197313, 930211, 614647, 683018, 610281]).astype(np.float)
ISSnegative = np.array([39986, 28817, 17690, 12754, 9576, 7191, 5572, 4624, 3988, 3950, 4469, 4905, 5066, 9397, 16391]).astype(np.float)
ISS_chargeconfused = ISSpositive*ccratio


#######################################################################################
c1 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
c1.SetLogy()

zeroline = [0]*rigidity.shape[0]

mg = TMultiGraph("mg","mg")
grshade = TGraph(2*rigidity.shape[0])
grshade2 = TGraph(2*rigidity.shape[0])
for i in range(rigidity.shape[0]):
  grshade.SetPoint(i,rigidity[i],ISSnegative[i])
  grshade.SetPoint(rigidity.shape[0]+i,rigidity[rigidity.shape[0]-i-1],ISS_chargeconfused[rigidity.shape[0]-i-1])
for i in range(rigidity.shape[0]):
  grshade2.SetPoint(i,rigidity[i],ISS_chargeconfused[i])
  grshade2.SetPoint(rigidity.shape[0]+i,rigidity[rigidity.shape[0]-i-1],zeroline[rigidity.shape[0]-i-1])

ccratio_plot = TGraph (rigidity.shape[0], rigidity,ccratio)
ISS_negative_plot = TGraph (rigidity.shape[0], rigidity, ISSnegative)
ISS_chargeconfused_plot = TGraph (rigidity.shape[0], rigidity, ISS_chargeconfused)
ISS_positive_plot = TGraph (rigidity.shape[0], rigidity, ISSpositive)

ccratio_plot.SetLineColor(1)
ISS_negative_plot.SetLineColor(4)
ISS_chargeconfused_plot.SetLineColor(6)
ISS_positive_plot.SetLineColor(8)
ccratio_plot.SetMarkerStyle(20)
ISS_negative_plot.SetMarkerStyle(21)
ISS_chargeconfused_plot.SetMarkerStyle(22)
ISS_positive_plot.SetMarkerStyle(23)

mg.SetTitle("; Rigidity (GV); ")
#mg.Add(ccratio_plot)
mg.Add(ISS_negative_plot)
mg.Add(ISS_positive_plot)
mg.Add(ISS_chargeconfused_plot)
mg.Draw("ACP")
grshade.SetFillStyle(3004)
grshade.SetFillColor(2)
grshade.Draw("F")
grshade2.SetFillStyle(3016)
grshade2.SetFillColor(1)
grshade2.Draw("F")


leg =TLegend(.7,.7,.9,.9,)
leg.SetFillColor(0)
#leg.AddEntry(ccratio_plot,"Charge Confusion Level","LP")
leg.AddEntry(ISS_negative_plot,"ISS negative data","LP")
leg.AddEntry(ISS_positive_plot,"ISS positive data","LP")
leg.AddEntry(ISS_chargeconfused_plot,"Charge Confused Event","LP")
leg.Draw()

c1.Update()
c1.SaveAs("ccratio.pdf")

#######################################################################################

c2 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
#c2.SetLogy()
cc_percentage = ISS_chargeconfused/ISSnegative
cc_percentage_plot = TGraph (rigidity.shape[0], rigidity, cc_percentage)
cc_percentage_plot.GetXaxis().SetTitle("Rigidity (GV)")
cc_percentage_plot.SetTitle("Percentage of charge confused proton in ISS Negative data")
#cc_percentage_plot.GetYaxis().SetLimits(0,1000);
#cc_percentage_plot.GetYaxis().SetRange(0,1000)
#cc_percentage_plot.GetYaxis().SetRange(0,1000)

cc_percentage_plot.Draw("AC*")
c2.Update()
c2.SaveAs("cc_percentage.pdf")



