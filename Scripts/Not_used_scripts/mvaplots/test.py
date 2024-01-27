import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist

Correct = np.load('negative_147_175_pattern_0.npy')

TH_Correct = TH1D("Correct","", 100,-5,5)
TH_Correct.SetFillColor(6)
TH_Correct.SetFillStyle(3004)
TH_Correct.SetLineColor(6)

fill_hist(TH_Correct, Correct[:,0])

scale = 150.0/TH_Correct.Integral()
TH_Correct.Scale(scale)

c1 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
TH_Correct.Draw()

leg =TLegend(.1,.7,.3,.9,)
leg.SetFillColor(0)
leg.AddEntry(TH_Correct,"ChargeCorrect")
leg.Draw()

TH_Correct.GetXaxis().SetTitle("L1L9RigidityMatching");

c1.Update()
c1.SaveAs("test.pdf")


