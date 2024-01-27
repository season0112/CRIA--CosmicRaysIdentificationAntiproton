import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist, root2array, tree2array
import os

bins = "330_525"

highpath = os.getenv('HPCHIGHENERGYDATADIR')

############### Load data ################################################################################
proton = root2array(highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_" + bins + ".root")
antiproton = root2array(highpath + "/B1042_antipr.pl1.1800_7.6_all_Tree_" + bins + ".root" )


'''
for i in range(16):
	TH_proton = TH1D("proton","", 100,-5,5)
	TH_antiproton = TH1D("antiproton","", 100,-5,5)
	TH_proton.SetFillColor(6)
	TH_proton.SetFillStyle(3004)
	TH_proton.SetLineColor(6)
	TH_antiproton.SetFillColor(4)
	TH_antiproton.SetFillStyle(3005)
	TH_antiproton.SetLineColor(4)

	fill_hist(TH_proton, proton[:,i])
	fill_hist(TH_antiproton, antiproton[:,i] )

	scale = 150.0/TH_proton.Integral()
	TH_proton.Scale(scale)
	scale = 150.0/TH_antiproton.Integral()
	TH_antiproton.Scale(scale)

	c1 = TCanvas()
	gPad.SetGrid()
	gPad.SetFrameFillColor(0)
	gStyle.SetOptStat("00000000")
	TH_proton.Draw()
	TH_antiproton.Draw("same")

	leg =TLegend(.1,.7,.3,.9,)
	leg.SetFillColor(0)
	leg.AddEntry(TH_proton,"Chargeproton")
	leg.AddEntry(TH_antiproton,"Chargeantiproton")
	leg.Draw()

	TH_proton.GetXaxis().SetTitle("L1L9RigidityMatching");

	c1.Update()
	c1.SaveAs("proton_antiproton_correct"+str(i)+".pdf")
'''

#####################################################################
plt.figure(figsize=(18,18))
plt.hist(proton['InnerMaxSpanRigidityMatching'], 80, range=(-3,3), log=True, alpha=0.5, label='Proton', facecolor='blue', edgecolor='black'  )
plt.hist(antiproton['InnerMaxSpanRigidityMatching'], 80, range=(-3,3), log=True, alpha=0.5, label='Antiproton', facecolor='blue', edgecolor='black'  )
plt.savefig(highpath + '/mvaplots/mva_4.png')


