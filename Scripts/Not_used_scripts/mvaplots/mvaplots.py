import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist

Correct = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1o9.2016000_7.6_all/147-1000GeV/transferdata/positive/147_175GeV/pattern0/positive_147_175_pattern_0.npy',encoding="latin1")
Confused = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/147-1000GeV/transferdata/negative/147_175GeV/pattern0/negative_147_175_pattern_0.npy',encoding="latin1")
'''
Correct = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/negative/175_211GeV/pattern0/negative_175_211_pattern_0.npy',encoding="latin1")
Confused = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/positive/175_211GeV/pattern0/positive_175_211_pattern_0.npy',encoding="latin1")
Correct2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/negative/211_250GeV/pattern0/negative_211_250_pattern_0.npy',encoding="latin1")
Confused2 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/positive/211_250GeV/pattern0/positive_211_250_pattern_0.npy',encoding="latin1")
Correct3 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/negative/250_330GeV/pattern0/negative_250_330_pattern_0.npy',encoding="latin1")
Confused3 = np.load('/hpcwork/jara0052/sichen/analysis_7.0/B1042_antipr.pl1.1800_7.6_all/147-1000GeV/transferdata/positive/250_330GeV/pattern0/positive_250_330_pattern_0.npy',encoding="latin1")
Correct = np.row_stack((Correct,Correct2))
Correct = np.row_stack((Correct,Correct3))
Confused = np.row_stack((Confused,Confused2))
Confused = np.row_stack((Confused,Confused3))
'''


for i in range(16):
    TH_Correct = TH1D("Correct","", 100,-5,5)
    TH_Confused = TH1D("Confused","", 100,-5,5)
    TH_Correct.SetFillColor(6)
    TH_Correct.SetFillStyle(3004)
    TH_Correct.SetLineColor(6)
    TH_Confused.SetFillColor(4)
    TH_Confused.SetFillStyle(3005)
    TH_Confused.SetLineColor(4)

    fill_hist(TH_Correct, Correct[:,i])
    fill_hist(TH_Confused, Confused[:,i] )

    scale = 150.0/TH_Correct.Integral()
    TH_Correct.Scale(scale)
    scale = 150.0/TH_Confused.Integral()
    TH_Confused.Scale(scale)

    c1 = TCanvas()
    gPad.SetGrid()
    gPad.SetFrameFillColor(0)
    gStyle.SetOptStat("00000000")
    TH_Correct.Draw()
    TH_Confused.Draw("same")

    leg =TLegend(.1,.7,.3,.9,)
    leg.SetFillColor(0)
    leg.AddEntry(TH_Correct,"ChargeCorrect")
    leg.AddEntry(TH_Confused,"ChargeConfused")
    leg.Draw()

    TH_Correct.GetXaxis().SetTitle("L1L9RigidityMatching");

    c1.Update()
    c1.SaveAs("proton"+str(i)+".pdf")


