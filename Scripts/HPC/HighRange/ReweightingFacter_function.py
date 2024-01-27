import numpy as np
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TChain, gDirectory, std, gROOT
from root_numpy import root2array, tree2array, fill_hist
import binning
import matplotlib.pyplot as plt
from ctypes import *
import uproot
from scipy.optimize import curve_fit

def LoadMCHist(highpath, binningWidth, LowEdge, HighEdge, BinNumber, pattern, MCName):

    PatternCut = "Pattern==" + str(pattern)

    print("load MC Root File.")

    if MCName == "B1042_pr.pl1.1800_7.6_all":
        MCROOTFILE = uproot.concatenate( [highpath + "/" + MCName + "_Tree_positive.root:ExampleAnalysisTree", highpath + "/" + MCName + "_Tree_negative.root:ExampleAnalysisTree"], ["MCPrimaryMomentum", "Weight"], cut=PatternCut, library='np')
        #MCROOTFILE = uproot.concatenate( [highpath + "/" + MCName + "_Tree_negative.root:ExampleAnalysisTree"], ["MCPrimaryMomentum", "Weight"], cut=PatternCut, library='np')

    elif MCName == "B1042_pr.pl1.flux.l1a9.2016000_7.6_all":
        #MCROOTFILE = uproot.concatenate( [highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/positive/results/reduced_root/" + "ExampleAnalysis_Tree_00000.root:ExampleAnalysisTree", highpath + "/" + MCName + "_Tree_negative.root:ExampleAnalysisTree"], ["MCPrimaryMomentum", "Weight"], cut=PatternCut, library='np')    
        #MCROOTFILE = uproot.concatenate( [ \
        #             highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/positive/results/reduced_root/" + "ExampleAnalysis_Tree_00000.root:ExampleAnalysisTree", \
        #             highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/positive/results/reduced_root/" + "ExampleAnalysis_Tree_00001.root:ExampleAnalysisTree", \
        #             highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/positive/results/reduced_root/" + "ExampleAnalysis_Tree_00002.root:ExampleAnalysisTree",], ["MCPrimaryMomentum", "Weight"], cut=PatternCut, library='np')
        MCROOTFILE = uproot.concatenate( [highpath + "/B1042_pr.pl1.flux.l1a9.2016000_7.6_all/positive/results/ExampleAnalysis_Tree_*_00020.root:ExampleAnalysisTree",], ["MCPrimaryMomentum", "Weight"], cut=PatternCut, library='np')

    elif MCName == "B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all":
        #MCROOTFILE = uproot.concatenate( [ \
        #             highpath + "/B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_525version/positive/results/reduced_root/" + "ExampleAnalysis_Tree_00000.root:ExampleAnalysisTree", \
        #             highpath + "/B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_525version/positive/results/reduced_root/" + "ExampleAnalysis_Tree_00001.root:ExampleAnalysisTree", \
        #             highpath + "/B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_525version/positive/results/reduced_root/" + "ExampleAnalysis_Tree_00002.root:ExampleAnalysisTree",], ["MCPrimaryMomentum", "Weight"], cut=PatternCut, library='np')
        MCROOTFILE = uproot.concatenate( [ highpath + "/B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_525version/positive/results/ExampleAnalysis_Tree_*_00020.root:ExampleAnalysisTree", highpath + "/B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_525version/positive/results/ExampleAnalysis_Tree_*_00050.root:ExampleAnalysisTree"], ["MCPrimaryMomentum","Weight"], cut=PatternCut, library='np')


    print("loading Complete.")

    #### Create Histogram with Merged Rigidity Bins and in 14.1-525GeV range
    h_MCPrimaryMomentum = TH1D("", "", 31, binning.Newbinnings_525_zhili[26:-1])
    #fill_hist(h_MCPrimaryMomentum, MCROOTFILE['MCPrimaryMomentum'], MCROOTFILE['Weight'])
    fill_hist(h_MCPrimaryMomentum, MCROOTFILE['MCPrimaryMomentum'])
    # Divide numbers over BinningWidth
    #print("h_MCPrimaryMomentum.GetNbinsX() is " + str(h_MCPrimaryMomentum.GetNbinsX()))
    for i in range(h_MCPrimaryMomentum.GetNbinsX()):
        h_MCPrimaryMomentum.SetBinContent(i+1, h_MCPrimaryMomentum.GetBinContent(i+1) / binningWidth[i])

    #### Create Histogram with fine Rigidity Bins and in large range (1 GV to 1.8 TV, maybe extent to higher.)
    h_MCPrimaryMomentum_more = TH1D("", "", BinNumber, LowEdge, HighEdge)
    fill_hist(h_MCPrimaryMomentum_more, MCROOTFILE['MCPrimaryMomentum'])
    #fill_hist(h_MCPrimaryMomentum_more, MCROOTFILE['MCPrimaryMomentum'], MCROOTFILE['Weight'])
    # Divide numbers over BinningWidth
    for i in range(h_MCPrimaryMomentum_more.GetNbinsX()):
        h_MCPrimaryMomentum_more.SetBinContent(i+1, h_MCPrimaryMomentum_more.GetBinContent(i+1) / ((HighEdge-LowEdge)/BinNumber) )

    return h_MCPrimaryMomentum, h_MCPrimaryMomentum_more, MCROOTFILE['Weight'], MCROOTFILE['MCPrimaryMomentum']


def LoadProtonRPL2015(resultpath, LowEdge, HighEdge, BinNumber):
    ProtonPRL2015  = TFile( "/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/ProtonFlux/AMSProtonflux2015.root" , "read")
    g_ProtonPRL2015 = ProtonPRL2015.Get("graph1")

    R_Proton    = []
    Flux_Proton = []
    for i in range(g_ProtonPRL2015.GetN()):
        R_Proton.append(g_ProtonPRL2015.GetX()[i])
        Flux_Proton.append(g_ProtonPRL2015.GetY()[i])

    R_Proton    = np.array(R_Proton)
    Flux_Proton = np.array(Flux_Proton)


    #### Create Histogram with Merged Rigidity Bins and in 14.1-525GeV range
    R_Proton_MergedPart    = []
    Flux_Proton_MergedPart = []
    for i in range(0, R_Proton[49:].shape[0]-1, 2):
        R_Proton_MergedPart.append( (R_Proton[49:][i] + R_Proton[49:][i+1]) / 2 )
        Flux_Proton_MergedPart.append( (Flux_Proton[49:][i] + Flux_Proton[49:][i+1]) / 2 )
    R_Proton_matched    = np.append(R_Proton[0:49]   , R_Proton_MergedPart[0:8])
    Flux_Proton_matched = np.append(Flux_Proton[0:49], Flux_Proton_MergedPart[0:8])
    R_Proton_matched    = R_Proton_matched[26:]
    Flux_Proton_matched = Flux_Proton_matched[26:]
    gROOT.cd();
    h_ProtonPRL = TH1D("", "", 31, binning.Newbinnings_525_zhili[26:-1])
    for i in range(h_ProtonPRL.GetNbinsX()):
        h_ProtonPRL.SetBinContent(i+1, Flux_Proton_matched[i])

    #### Create Histogram with fine Rigidity Bins and in large range (1 GV to 1.8 TV, maybe extent to higher.)
    def powerlaw(E, C, gamma):
        return C * np.power(E, gamma) 

    popt, pcov = curve_fit(powerlaw, R_Proton[26:], Flux_Proton[26:], maxfev=10000) # R_Proton[26]=14.7

    #R_Proton_More = binning.ProtonBinCenter
    R_Proton_More = np.arange( LowEdge, HighEdge, (HighEdge-LowEdge)/BinNumber ) # Here R_Proton_More could be higher, but depends on MC Generate Momentum range.
    Y_Fit = powerlaw(R_Proton_More, *popt) 

    '''
    ## Check Fit
    fig, ax = plt.subplots(figsize=(30,18))
    plt.errorbar(R_Proton,  Flux_Proton, yerr=0, marker='o',linestyle="None",markersize=12, markerfacecolor="blue",ecolor="blue" )
    plt.errorbar(R_Proton_More, Y_Fit, yerr=0, marker='o',linestyle="None",markersize=12, markerfacecolor="red",ecolor="red" )
    ax.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    plt.ylabel("",fontsize=50)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + "check" + ".pdf")
    '''

    h_ProtonPRL_PowerLawFit = TH1D("", "", BinNumber, LowEdge, HighEdge)
    for i in range(h_ProtonPRL_PowerLawFit.GetNbinsX()):
        h_ProtonPRL_PowerLawFit.SetBinContent(i+1, Y_Fit[i])

    h_ProtonPRL_AllRange = TH1D("", "", Flux_Proton.shape[0], binning.ProtonBinEdge)
    for i in range(h_ProtonPRL_AllRange.GetNbinsX()):
        h_ProtonPRL_AllRange.SetBinContent(i+1, Flux_Proton[i])

    return h_ProtonPRL, h_ProtonPRL_PowerLawFit, h_ProtonPRL_AllRange


def Plot_MCPrimaryMomentum(resultpath, MCName, pattern, h_MCPrimaryMomentum, h_MCPrimaryMomentum_more):
    c_acc = TCanvas("c_acc","c_acc",1000,500)

    h_MCPrimaryMomentum.Draw("HIST")
    h_MCPrimaryMomentum_more.Draw("same HIST")

    h_MCPrimaryMomentum.SetLineColor(2)
    h_MCPrimaryMomentum_more.SetLineColor(3)

    gPad.SetLogx();
    gPad.SetLogy();

    xaxis = h_MCPrimaryMomentum.GetXaxis()
    xaxis.SetMoreLogLabels();
    xaxis.SetMoreLogLabels();
    c_acc.SaveAs(resultpath + '/antiproton_to_proton_ratio_plot/' + "MCPrimaryMomentum_" + MCName + "_Pattern_" + str(pattern) + ".pdf");

def Plot_ProtonReference(resultpath, MCName, pattern, h_ProtonPRL, h_ProtonPRL_PowerLawFit, h_ProtonPRL_AllRange):
    c_acc = TCanvas("c_acc","c_acc",1000,500)

    h_ProtonPRL_AllRange.Draw("HIST")
    h_ProtonPRL_PowerLawFit.Draw("HIST same")
    h_ProtonPRL.Draw("HIST same ")
    #h_ProtonPRL_AllRange.Draw("HIST same")

    h_ProtonPRL.SetLineColor(2)
    h_ProtonPRL_PowerLawFit.SetLineColor(3)
    h_ProtonPRL_AllRange.SetLineColor(4)

    gPad.SetLogx();
    gPad.SetLogy();
    xaxis = h_ProtonPRL.GetXaxis()
    xaxis.SetMoreLogLabels();
    xaxis.SetMoreLogLabels();
    c_acc.SaveAs(resultpath + '/antiproton_to_proton_ratio_plot/' + "ProtonReference" + ".pdf");



def Plot_ReweightFactor(resultpath, MCName, pattern, ReweightingFactor, binname):
    c_acc = TCanvas("c_acc","c_acc",1000,500)
    #h_MCPrimaryMomentum.Draw("")
    #g_ProtonPRL2015.Draw("same")

    ReweightingFactor.Draw("")

    gPad.SetLogy();

    xaxis = ReweightingFactor.GetXaxis()
    xaxis.SetMoreLogLabels();
    c_acc.SaveAs(resultpath + '/antiproton_to_proton_ratio_plot/' + "Reweight_Factor_" + binname + "_" + MCName + "_Pattern_" + str(pattern) + ".pdf");


def Plot_MCWeight(resultpath, MCName, pattern, MCWeight, MCPrimaryMomentum, LowEdge, HighEdge, BinNumber):

    WeightNumber = 500
    H, xedges, yedges = np.histogram2d(MCPrimaryMomentum, MCWeight, bins=[BinNumber, WeightNumber], range=[ [LowEdge,HighEdge], [min(MCWeight), max(MCWeight)] ])
    #print(H.shape)
    #print(H)

    width = (HighEdge-LowEdge)/BinNumber
    a=[] 
    b=[]
    for i in range(BinNumber):
        a.append( sum(H[i,:]) )
        b.append( width/2 + width*i)

    fig, ax = plt.subplots(figsize=(30,18))
    #plt.hist(MCWeight, bins=BinNumber, range=(LowEdge,HighEdge), density=False, facecolor='g', alpha=0.75)
    plt.errorbar(b, a, yerr=0, fmt='o', markersize=20, color='r', label='')
    ax.tick_params(direction='in',length=10, width=3)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=80)
    plt.ylabel("",fontsize=50)
    #plt.xscale('log')
    plt.yscale('log')
    plt.savefig(resultpath + '/antiproton_to_proton_ratio_plot/' + "MCWeight_" + MCName + "_Pattern_" + str(pattern) + ".pdf")


















