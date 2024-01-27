import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist

correct_proton = np.load('positive_250_330_pattern_0.npy')
correct_antiproton = np.load('negative_250_330_pattern_0.npy')

TH_Proton = TH1D("correct_proton","", 100,-10,10)
TH_Antiproton = TH1D("correct_antiproton","", 100,-10,10)
TH_Proton.SetFillColor(6)
TH_Proton.SetFillStyle(3004)
TH_Proton.SetLineColor(6)
TH_Antiproton.SetFillColor(4)
TH_Antiproton.SetFillStyle(3005)
TH_Antiproton.SetLineColor(4)

fill_hist(TH_Proton, correct_proton[:,9])
fill_hist(TH_Antiproton, correct_antiproton[:,9] )

scale = 150.0/TH_Proton.Integral()
TH_Proton.Scale(scale)
scale = 150.0/TH_Antiproton.Integral()
TH_Antiproton.Scale(scale)

c1 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
TH_Proton.Draw()
TH_Antiproton.Draw("same")

leg =TLegend(.1,.7,.3,.9,)
leg.SetFillColor(0)
leg.AddEntry(TH_Proton,"correct_proton")
leg.AddEntry(TH_Antiproton,"correct_antiproton")
leg.Draw()

TH_Proton.GetXaxis().SetTitle("L24L58RigidityMatching");

c1.Update()
c1.SaveAs("proton"+str(9)+".pdf")

'''
plt.figure(figsize=(18,18))
plt.hist(correct_proton[:,9],80,  alpha=0.5, range=(-10,10), density=True,  label='positive',facecolor='blue',edgecolor='black')
plt.hist(correct_antiproton[:,9],80,  alpha=0.5, range=(-10,10), density=True,  label='negative',facecolor='red',edgecolor='green')
plt.savefig('9'+'.png')
'''
