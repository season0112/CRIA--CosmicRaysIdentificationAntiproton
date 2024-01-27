import numpy as np
import uproot
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, TH1F, TH2F, gStyle
import matplotlib.pyplot as plt
import PythonPlotDefaultParameters

def Plot_RatiovsSignalEfficiency(totaleff, ratioall, binname, RigidityRange, TimeMode, TimeSplitMode, TimeIndex):

    plt.figure(figsize=(30,18))
    ax=plt.gca()

    plt.plot(totaleff, ratioall, marker='o', color="blue", lw=0, markersize=20, label='')

    yaxisfactor = 1.5
    ax.set_ylim( (plt.ylim()[1]+plt.ylim()[0])/2 - ((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor, (plt.ylim()[1]+plt.ylim()[0])/2 + ((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor )

    plt.xlabel('Signal Efficiency'                          , fontsize=60, horizontalalignment='right', x=1.0)
    plt.ylabel(r'$\rm{\overline{p}/p}$ ($\times$ $10^{5}$)', fontsize=60, horizontalalignment='right', y=1.0)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    ax.tick_params(axis='both', which='both', direction='in', length=30, width=10)

    plt.savefig('Ratio_vs_SignalEff_' + RigidityRange + '_' + binname + '_' + TimeMode + '_' + str(TimeSplitMode) + '_' + str(TimeIndex) + '.pdf')
    plt.close()


def Plot_RMSvsRigidity(BinningCenterAll, SystemError_SigEff, RigidityRange, TimeMode, TimeSplitMode, TimeIndex):

    plt.figure(figsize=(18,9))
    ax=plt.gca()

    plt.plot(BinningCenterAll, SystemError_SigEff, marker='o', color="blue", lw=0, markersize=15, label='')

    yaxisfactor = 1.5
    ax.set_ylim( (plt.ylim()[1]+plt.ylim()[0])/2 - ((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor, (plt.ylim()[1]+plt.ylim()[0])/2 + ((plt.ylim()[1]-plt.ylim()[0])/2)*yaxisfactor )

    plt.xlabel('Rigidity (GV)', fontsize=30)
    plt.ylabel('RMS'          , fontsize=35)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    #plt.legend( loc='best',fontsize=20)
    plt.savefig('RMS_vs_Rigidity_' + RigidityRange + '_' + TimeMode + '_' + str(TimeSplitMode) + '_' + str(TimeIndex) + '.pdf')
    plt.close()


def Plot_RMS(hist_ratio, binname, RigidityRange, TimeMode, TimeSplitMode, TimeIndex):

    c2 = TCanvas("c2", "c2", 1000, 800);

    gStyle.SetPadBorderMode(0);
    gStyle.SetFrameBorderMode(0);
    gStyle.SetTitleFontSize(0.1);
    gStyle.SetTitleFont(62);
    gStyle.SetOptFit(1)

    hist_ratio.Draw("")
    hist_ratio.SetLineColor(4)
    hist_ratio.SetLineWidth(2)

    gPad.Update()

    st = hist_ratio.FindObject("stats")
    #gStyle.SetOptStat(111111111)
    st.SetX1NDC(0.1) # new x start position
    st.SetX2NDC(0.4) # new x end position
    st.SetY1NDC(0.7)
    st.SetY2NDC(0.9)
    st.SetTextSize(0.045);

    xaxis = hist_ratio.GetXaxis();
    yaxis = hist_ratio.GetYaxis();
    xaxis.SetTitle("#bar{p}/p (#times 10^{5})");
    xaxis.SetTitleFont(62);
    yaxis.SetTitleFont(62);
    xaxis.SetTitleSize(0.045);
    yaxis.SetTitleSize(0.045);
    xaxis.SetLabelFont(62);
    xaxis.SetLabelSize(0.05);
    yaxis.SetLabelFont(62);
    yaxis.SetLabelSize(0.05);

    gPad.SetBottomMargin(0.15);
    xaxis.SetTitleOffset(1.0);

    c2.SaveAs("RMS_" + RigidityRange + '_' + binname + '_' + TimeMode + '_' + str(TimeSplitMode) + '_' + str(TimeIndex) + "_LinearY.pdf");
    gPad.SetLogy();
    c2.SaveAs("RMS_" + RigidityRange + '_' + binname + '_' + TimeMode + '_' + str(TimeSplitMode) + '_' + str(TimeIndex) + "_LogY.pdf");
    c2.Close();





