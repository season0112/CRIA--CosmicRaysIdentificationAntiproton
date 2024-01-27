import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F, gStyle
import binning
import argparse
import os
import matplotlib.pyplot as plt
import PythonPlotDefaultParameters
import CalculateSysErrorFromShape_tool

def Plot_ratio_histogram(hist_ratio, BinningAll, RIndex, RigidityRange):
    c2 = TCanvas("c2", "c2", 1000, 500);
    gStyle.SetPadBorderMode(0);
    gStyle.SetFrameBorderMode(0);
    gStyle.SetTitleFontSize(0.1);
    gStyle.SetTitleFont(62);
    hist_ratio.Draw("")
    hist_ratio.SetTitle('Antiproton to proton ratio in ' + str(BinningAll[RIndex]) + ' GV');
    xaxis = hist_ratio.GetXaxis();
    yaxis = hist_ratio.GetYaxis();
    xaxis.SetTitle("Antiproton to proton ratio (*10^-5)");
    xaxis.SetTitleFont(62);
    yaxis.SetTitleFont(62);
    xaxis.SetTitleSize(0.05);
    yaxis.SetTitleSize(0.04);
    xaxis.SetLabelFont(62);
    xaxis.SetLabelSize(0.05);
    yaxis.SetLabelFont(62);
    yaxis.SetLabelSize(0.05);
    gPad.SetBottomMargin(0.15);
    xaxis.SetTitleOffset(1.0);
    c2.SaveAs("RMSPlots/RMS_" + str(BinningAll[RIndex]) + 'GV_' + RigidityRange + "_LinearY.pdf")
    gPad.SetLogy();
    c2.SaveAs("RMSPlots/RMS_" + str(BinningAll[RIndex]) + 'GV_' + RigidityRange + "_LogY.pdf")
    c2.Close();


def Plot_RMS(Covnumber, BinningCenterAll, SysErrorShape_original, RigidityRange):
    plt.figure(figsize=(18,9))
    ax=plt.gca()
    #ax.axes.set_ylim(0, 0.70)
    plt.plot(np.array(BinningCenterAll), np.array(SysErrorShape_original)*100000, "o", markersize=15)
    plt.plot(BinningCenterAll, CalculateSysErrorFromShape_tool.smooth(SysErrorShape_original, Covnumber)*100000, "*", markersize=10) ##Parametrization (should be updated with Parametrization of template)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.xlabel('Rigidity (GV)',fontsize=30)
    plt.ylabel('RMS (*10^-5)',fontsize=30)
    plt.savefig('RMSPlots/RMSPlots_' + 'Covnumber'+ str(Covnumber) + "_" + RigidityRange + '.pdf')
    plt.close()


