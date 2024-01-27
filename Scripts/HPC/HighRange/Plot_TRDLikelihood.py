import numpy as np
import os
import matplotlib.pyplot as plt
import binning
import uproot
from root_numpy import fill_hist
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TGraphErrors

def RemoveXError(grapherror):
    for i in range(grapherror.GetN()):
        grapherror.SetPointError(i, 0, grapherror.GetErrorY(i))


def Plot(h_proton, h_electron, h_proton_ISS, h_electron_ISS, g_proton_ISS, g_electron_ISS, rigiditybin, pattern, protonMCName, binnumber):

    c2 = TCanvas("c2", "c2", 800, 600);

    gStyle.SetPadBorderMode(0);
    gStyle.SetFrameBorderMode(0);
    gStyle.SetOptStat(0)

    h_electron    .Draw("HIST")
    h_proton      .Draw("HIST same")
    #h_proton_ISS  .Draw("HIST same")
    #h_electron_ISS.Draw("HIST same")
    g_proton_ISS  .Draw("P")
    g_electron_ISS.Draw("P")

    h_electron    .SetLineColor(2) # red
    h_proton      .SetLineColor(4) # blue
    #h_electron_ISS.SetLineColor(3) # green
    #h_proton_ISS  .SetLineColor(6) # pink
    g_electron_ISS.SetLineColor(2)
    g_proton_ISS  .SetLineColor(4)

    #g_proton_ISS  .
    #g_electron_ISS.
    g_proton_ISS  .SetMarkerStyle(15)
    g_electron_ISS.SetMarkerStyle(15)
    g_electron_ISS.SetMarkerColor(2)
    g_proton_ISS  .SetMarkerColor(4)
    g_electron_ISS.SetMarkerSize(0.9)
    g_proton_ISS  .SetMarkerSize(0.9)

    xaxis = h_electron.GetXaxis();
    yaxis = h_electron.GetYaxis();
    xaxis.SetTitle("#Lambda_{TRD}"); # $\Lambda_\mathrm{TRD}$
    yaxis.SetTitle("Probability");
    xaxis.SetTitleFont(62);
    yaxis.SetTitleFont(62);
    xaxis.SetTitleSize(0.045);
    yaxis.SetTitleSize(0.045);
    xaxis.SetLabelFont(62);
    xaxis.SetLabelSize(0.05);
    yaxis.SetLabelFont(62);
    yaxis.SetLabelSize(0.05);

    #gPad.SetLeftMargin(0.16);
    #gPad.SetBottomMargin(0.15);
    #xaxis.SetTitleOffset(1.0);
    yaxis.SetTitleOffset(1.1);

    legend1 = TLegend(0.7, 0.65, 0.85, 0.85);
    legend1.AddEntry(h_electron    , "Electron (MC)" , "lpf");
    legend1.AddEntry(h_proton      , "Proton (MC)"   , "lpf");
    legend1.AddEntry(g_electron_ISS, "Electron (ISS)", "p");
    legend1.AddEntry(g_proton_ISS  , "Proton (ISS)"  , "p");
    legend1.SetTextSize(0.04);
    legend1.SetTextFont(62);
    legend1.SetBorderSize(0);
    legend1.Draw();

    #c2.SaveAs("TRD_" + protonMCName + "rigidity_" + rigiditybin + "_pattern_" + pattern + ".pdf")
    gPad.SetLogy();
    c2.SaveAs("TRD_" + protonMCName + "rigidity_" + rigiditybin + "_pattern_" + pattern + "_bin_" + str(binnumber) + "_LogY.pdf")
    c2.Close();

def main():
    #rigiditybin = '21.1_22.8'
    #rigiditybin = '93_108'
    #rigiditybin = '175_211'
    #rigiditybin_all = ['21.1_22.8', '93_108', '175_211']
    #rigiditybin_all = ['175_211']
    #rigiditybin_all = ['21.1_22.8']
    rigiditybin_all = ['14.1_15.3']
    #bins_525_zhili
    #binslow

    #pattern_all = ['-1', '0', '1', '2', '3', '4', '5']
    pattern_all = ['0']

    #protonMCName_all = ['ChargeCorrectProtonTemplate_MC_B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_Tree_positive_', 'ChargeCorrectProtonTemplate_MC_B1042_pr.pl1.1800_7.6_all_Tree_positive_']
    protonMCName_all = ['ChargeCorrectProtonTemplate_MC_B1042_pr.pl1.1800_7.6_all_Tree_positive_']

    binnumber = 100
  
    for protonMCName in protonMCName_all:
        for rigiditybin in rigiditybin_all:
            for pattern in pattern_all:
                #protonMCName  = "ChargeCorrectProtonTemplate_MC_B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_Tree_positive_"
                #protonMCName   = "ChargeCorrectProtonTemplate_MC_B1042_pr.pl1.1800_7.6_all_Tree_positive_"
                electronMCName = "ElectronTemplate_MC_B1091_el.pl1.0_25_2000_7.6_all_Tree_"

                proton_MC    = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/" + protonMCName   + rigiditybin + "_Pattern_" + pattern + ".npy")
                electron_MC  = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/" + electronMCName + rigiditybin + "_Pattern_" + pattern + ".npy")
                proton_ISS   = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/ISS_positive_" + rigiditybin + "_Pattern_" + pattern + ".npy")
                electron_ISS = np.load("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v3.0/ISS_anylsis/data/ElectronTemplate_Data__" + rigiditybin + "_Pattern_" + pattern + ".npy")
                electron_MC  = electron_MC[np.where(electron_MC[:,0]>0.80)[0]]
                electron_ISS = electron_ISS[np.where(electron_ISS[:,0]>0.80)[0]]

                h_proton_MC    = TH1D("Proton (MC)"   , "", binnumber, 0, 2.0)
                h_electron_MC  = TH1D("Electron (MC)" , "", binnumber, 0, 2.0)
                h_proton_ISS   = TH1D("Proton (ISS)"  , "", binnumber, 0, 2.0)
                h_electron_ISS = TH1D("Electron (ISS)", "", binnumber, 0, 2.0)

                fill_hist(h_proton_MC   , proton_MC[:, 1])
                fill_hist(h_electron_MC , electron_MC[:, 1])
                fill_hist(h_proton_ISS  , proton_ISS[:, 1])
                fill_hist(h_electron_ISS, electron_ISS[:, 1])

                scale_proton_MC = 1/h_proton_MC.Integral()
                h_proton_MC.Scale(scale_proton_MC)

                scale_electron_MC = 1/h_electron_MC.Integral()
                h_electron_MC.Scale(scale_electron_MC)

                scale_proton_ISS = 1/h_proton_ISS.Integral()
                h_proton_ISS.Scale(scale_proton_ISS)

                scale_electron_ISS = 1/h_electron_ISS.Integral()
                h_electron_ISS.Scale(scale_electron_ISS)

                g_proton_ISS   = TGraphErrors(h_proton_ISS)
                g_electron_ISS = TGraphErrors(h_electron_ISS)

                RemoveXError(g_proton_ISS)
                RemoveXError(g_electron_ISS)

                Plot(h_proton_MC, h_electron_MC, h_proton_ISS, h_electron_ISS, g_proton_ISS, g_electron_ISS, rigiditybin, pattern, protonMCName, binnumber)

if __name__ == "__main__":
    main()


