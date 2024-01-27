import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TGraphErrors, TAxis, TPad
import binning

def GetMVAvariables(FullNumpyData, sign):

    MVANumpy = np.array([FullNumpyData['TrdLogLikelihoodRatioElectronProtonTracker'], FullNumpyData['RigidityAsymmetry'], FullNumpyData['RigidityAsymmetryL9'], FullNumpyData['Chi2TrackerYAsymmetry'], FullNumpyData['InnerMaxSpanRigidityMatching']*sign, FullNumpyData['L1L9RigidityMatching']*sign, FullNumpyData['L24L58RigidityMatching']*sign, FullNumpyData['Log10Chi2TrackerXInner'], FullNumpyData['Log10Chi2TrackerYInner'], FullNumpyData['Log10Chi2TrackerX'], FullNumpyData['Log10Chi2TrackerY'], FullNumpyData['TrackerL58L24ChargeAsymmetry'], FullNumpyData['TrackerL9Charge'], FullNumpyData['TrackerL78Charge'], FullNumpyData['UpperTofCharge'], FullNumpyData['LowerTofCharge'], FullNumpyData['Weight']])

    return MVANumpy




def PlotMVAvariables(highpath, i, binmerge, namelabel, namelabel_symbol, index, TH_Correct_MC, TH_Confused_MC, g_Correct_data, g_Confused_data, MC_name, legend_xlow, legend_xhigh, legend_ylow, legend_yhigh, Pattern):

    c1 = TCanvas()

    p1 = TPad("", "", 0.0, 0.0, 1.0, 1.0, 0) #(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1.Draw()
    p1.cd()

    gPad.SetFrameFillColor(0)
    gStyle.SetOptStat("00000000")

    TH_Correct_MC.SetFillColor(4)
    TH_Correct_MC.SetFillStyle(3004)
    TH_Correct_MC.SetLineColor(4)

    TH_Confused_MC.SetFillColor(2)
    TH_Confused_MC.SetFillStyle(3005)
    TH_Confused_MC.SetLineColor(2)

    g_Correct_data.SetMarkerStyle(15)
    g_Correct_data.SetMarkerColor(4)
    g_Correct_data.SetLineColor(4)

    g_Confused_data.SetMarkerStyle(15)
    g_Confused_data.SetMarkerColor(2)
    g_Confused_data.SetLineColor(2)
    gStyle.SetErrorX(0)

    if index == 5 or index == 1 or index == 2:  
        TH_Correct_MC.Draw("HIST")
        TH_Confused_MC.Draw("HIST same")
    else:
        TH_Confused_MC.Draw("HIST")
        TH_Correct_MC.Draw("HIST same")

    TH_Correct_MC.SetLineWidth(1)
    TH_Confused_MC.SetLineWidth(1)

    #TH_Correct_data.Draw("HIST same")
    #TH_Confused_data.Draw("HIST same")
    g_Correct_data.Draw("E1 P")
    g_Confused_data.Draw("E1 P")

    leg =TLegend(legend_xlow, legend_ylow, legend_xhigh, legend_yhigh)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.SetTextFont(62)
    leg.AddEntry(TH_Correct_MC    , "Charge Correct (MC)")
    leg.AddEntry(TH_Confused_MC   , "Charge Confused (MC)")
    #leg.AddEntry(TH_Correct_data , "Charge Correct_data")
    #leg.AddEntry(TH_Confused_data, "Charge Confused_data")
    leg.AddEntry(g_Correct_data   , "Charge Correct (ISS)")
    leg.AddEntry(g_Confused_data  , "Charge Confused (ISS)")
    leg.Draw()

    if index == 5 or index == 1 or index == 2:
        xaxis = TH_Correct_MC.GetXaxis()
        yaxis = TH_Correct_MC.GetYaxis()
    else:
        xaxis = TH_Confused_MC.GetXaxis()
        yaxis = TH_Confused_MC.GetYaxis()
    xaxis.SetTitle(namelabel_symbol[index])
    yaxis.SetTitle("Normalized events")
    xaxis.SetTitleFont(62)
    xaxis.SetTitleSize(0.045) # 0.045
    xaxis.SetLabelFont(62)
    xaxis.SetLabelSize(0.05)
    yaxis.SetTitleFont(62)
    yaxis.SetTitleSize(0.045)
    yaxis.SetLabelFont(62)
    yaxis.SetLabelSize(0.05)

    p1.SetLeftMargin(0.16);
    p1.SetBottomMargin(0.17);
    xaxis.SetTitleOffset(1.2)

    #c1.SaveAs(highpath + "/mvaplots/proton_" + binning.bins[i] + "_to_" + binning.bins[i+binmerge-1] + str("_") + namelabel[index] + str("_") + str(index) + str("_") + MC_name + ".pdf")
    gPad.SetLogy();
    c1.Update()
    c1.SaveAs(highpath + "/mvaplots/proton_" + binning.bins[i] + "_to_" + binning.bins[i+binmerge-1] + str("_") + namelabel[index] + str("_") + str(index) + "_Logy" + str("_") + MC_name + "_Pattern_" + Pattern + ".pdf")

    return c1













