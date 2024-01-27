import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend
from root_numpy import fill_hist

MC1042=np.load('positive_147_175_pattern_0.npy')
ISS=np.load('positive_147_175_pattern_0_ISS.npy')


TH_MC1042 = TH1D("MC1042","", 100,0.4,1.4)
TH_ISS = TH1D("ISS","ISS", 100,0.4,1.4)
TH_MC1042.SetFillColor(6)
TH_MC1042.SetFillStyle(3004)
TH_MC1042.SetLineColor(6)
TH_ISS.SetFillColor(4)
TH_ISS.SetFillStyle(3005)
TH_ISS.SetLineColor(4)

fill_hist(TH_MC1042, MC1042[:,0])
fill_hist(TH_ISS, ISS[:,0] )

scale = 150/TH_MC1042.Integral()
TH_MC1042.Scale(scale)
scale = 150/TH_ISS.Integral()
TH_ISS.Scale(scale)

c1 = TCanvas()
gPad.SetGrid()
gPad.SetFrameFillColor(0)
gStyle.SetOptStat("00000000")
TH_MC1042.Draw()
TH_ISS.Draw("same")

leg =TLegend(.1,.7,.3,.9,)
leg.SetFillColor(0)
leg.AddEntry(TH_MC1042,"MCBproton1042")
leg.AddEntry(TH_ISS,"ISS")
leg.Draw()

TH_MC1042.GetXaxis().SetTitle("TrdLogLikelihood(e/p)");

c1.Update()
c1.SaveAs("plot_root.emf")


'''
plt.hist(MC1042[:,0],bins=80,range=(0,1.5),log=True,normed=True,alpha=0.5,edgecolor='black',facecolor='blue',label='MC:B1024proton' )
plt.hist(ISS[:,0],bins=80,range=(0,1.5),log=True,normed=True,alpha=0.5,edgecolor='black',facecolor='green',label='ISS' )
plt.legend(loc='upper left',fontsize=10)
plt.savefig('plot.png')
'''
