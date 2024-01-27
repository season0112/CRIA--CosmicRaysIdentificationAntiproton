import numpy as np
import matplotlib.pyplot as plt
from ROOT import TH2D, TH1D, TCanvas, gPad, TF1, TGraph, gEnv, TFile, gStyle, TLegend, TGraphErrors, TAxis, TPad
from root_numpy import fill_hist, root2array, tree2array


def LoadFitResult(highpath, CCCut, CCN, TRDN):
    with open (highpath + '/ISS_anylsis/data/templatefit_pass7.8/negative/results_525version' + '/fit_resultscccut' + CCCut + 'CCN' + CCN + 'TRDN' + TRDN + '.txt') as fileresult:
        a_number = fileresult.readlines()
    a_number = map(lambda s: float(s.split(',')[0]) ,a_number)
    a_number = list(a_number)
    a_number = np.array(a_number)
    return a_number


def draw_Correct_1D(Correct_data, binnumber, Correct_MC, MarkerSize):
    hist, bins = np.histogram(Correct_data[:,0], bins=binnumber, range=(0, 1) )
    center     = (bins[:-1] + bins[1:]) / 2
    factor     = 1 / Correct_data[:,0].shape[0] * Correct_MC[:,0].shape[0]

    plt.errorbar(center, hist * factor, yerr=np.sqrt(hist) * factor, label="Charge Correct (ISS)", fmt='o', markersize=MarkerSize, color="blue")
    plt.hist( Correct_MC[:,0], bins=binnumber, alpha=1.0, label='Charge Correct (MC)', facecolor='None', edgecolor='blue', linewidth=10, histtype='step')


def draw_Confused_1D(Confused_data, binnumber, Correct_MC, Confused_MC, ccpercentage, MarkerSize, factor_CorrectConfusedCountRatio):

    # ISS Data
    hist, bins                     = np.histogram(Confused_data[:,0], bins=binnumber, range=(0, 1) )
    center                         = (bins[:-1] + bins[1:]) / 2
    hist_CorrectMC, bins_CorrectMC = np.histogram( Correct_MC[:,0], bins=binnumber, range=(0, 1) )
    hist_realCC                    = ( hist/Confused_data[:,0].shape[0] - ccpercentage * hist_CorrectMC/Correct_MC[:,0].shape[0] ) / (1-ccpercentage)
    factor_realCC                  = 1 * Confused_MC[:,0].shape[0]

    plt.errorbar(center, hist_realCC * factor_realCC * factor_CorrectConfusedCountRatio, yerr=np.sqrt(hist_realCC * factor_realCC * factor_CorrectConfusedCountRatio), label="Charge Confused (ISS)", fmt='o', markersize=MarkerSize, color="red")

    # MC
    plt.hist( Confused_MC[:,0].repeat(factor_CorrectConfusedCountRatio), bins=binnumber, alpha=1.0, label='Charge Confused (MC)', facecolor='None', edgecolor='red', linewidth=10, histtype='step')


def plot_Correct_1D(Correct_data, binnumber, Correct_MC, highpath, Rigiditybin, i, MarkerSize):
    plt.figure(figsize=(30,18))
    ax=plt.gca()

    draw_Correct_1D(Correct_data, binnumber, Correct_MC, MarkerSize)

    plt.legend(loc='best',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    ax.tick_params(direction='in',length=10, width=3, axis='y', which='minor')
    ax.tick_params(direction='in',length=20, width=6, axis='y', which='major')
    plt.xlabel(r"$\Lambda_{CC}$",fontsize=40)

    plt.savefig(highpath + "/CCEstimatorplots/Correct_MC_1D_" + str(Rigiditybin[i]) + ".pdf")
    plt.yscale('log')
    plt.savefig(highpath + "/CCEstimatorplots/Correct_MC_1D_Logy_" + str(Rigiditybin[i]) + ".pdf")
    plt.close()


def plot_Confused_1D(Confused_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, MarkerSize, factor_CorrectConfusedCountRatio):
    plt.figure(figsize=(30,18))
    ax=plt.gca()

    draw_Confused_1D(Confused_data, binnumber, Correct_MC, Confused_MC, ccpercentage, MarkerSize, factor_CorrectConfusedCountRatio)

    plt.legend(loc='best',fontsize=60)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    ax.tick_params(direction='in',length=10, width=3, axis='y', which='minor')
    ax.tick_params(direction='in',length=20, width=6, axis='y', which='major')
    plt.xlabel(r"$\Lambda_{CC}$",fontsize=40)

    plt.savefig(highpath + "/CCEstimatorplots/Confused_MC_1D_" + str(Rigiditybin[i]) + ".pdf")
    plt.yscale('log')
    plt.savefig(highpath + "/CCEstimatorplots/Confused_MC_1D_Logy_" + str(Rigiditybin[i]) + ".pdf")
    plt.close()


def plot_CorrectAndConfused_1D(Confused_data, Correct_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, MarkerSize, factor_CorrectConfusedCountRatio):
    plt.figure(figsize=(30,18))
    ax=plt.gca()

    ## Correct
    draw_Correct_1D(Correct_data, binnumber, Correct_MC, MarkerSize)
    ## Confused
    draw_Confused_1D(Confused_data, binnumber, Correct_MC, Confused_MC, ccpercentage, MarkerSize, factor_CorrectConfusedCountRatio)

    plt.legend(loc='best',fontsize=40)
    plt.xticks(fontsize=50)
    plt.yticks(fontsize=50)
    ax.tick_params(direction='in',length=10, width=3, axis='y', which='minor')
    ax.tick_params(direction='in',length=20, width=6, axis='y', which='major')
    plt.xlabel(r"$\Lambda_{CC}$",fontsize=40)

    plt.savefig(highpath + "/CCEstimatorplots/CorrectAndConfused_1D_" + str(Rigiditybin[i]) + ".pdf")
    plt.yscale('log')
    plt.savefig(highpath + "/CCEstimatorplots/CorrectAndConfused_1D_Logy_" + str(Rigiditybin[i]) + ".pdf")
    plt.close()


def plotRoot_CorrectAndConfused_1D(Confused_data, Correct_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, factor_CorrectConfusedCountRatio, CCCut, CCN, TRDN):

    c1 = TCanvas("", "", 800, 600)

    p1 = TPad("", "", 0.0, 0.0, 1.0, 1.0, 0) #(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1.Draw()
    p1.cd()

    gPad.SetFrameFillColor(0)
    gStyle.SetOptStat("00000000")

    TH_Correct_MC    = TH1D("Correct_MC"   , "", binnumber, 0, 1)
    TH_Confused_MC   = TH1D("Confused_MC"  , "", binnumber, 0, 1)

    fill_hist(TH_Correct_MC   , Correct_MC[:,0]   )
    fill_hist(TH_Confused_MC  , Confused_MC[:,0].repeat(factor_CorrectConfusedCountRatio)  )

    ## Correct Data
    hist, bins = np.histogram(Correct_data[:,0], bins=binnumber, range=(0, 1) )
    center     = (bins[:-1] + bins[1:]) / 2
    factor     = 1 / Correct_data[:,0].shape[0] * Correct_MC[:,0].shape[0]
    g_Correct_data  = TGraphErrors(len(center), np.array(center.tolist()), np.array(hist)*factor, np.array([0]*len(center)), np.sqrt(hist)*factor)

    ## Confused Data
    hist, bins                     = np.histogram(Confused_data[:,0], bins=binnumber, range=(0, 1) )
    center                         = (bins[:-1] + bins[1:]) / 2
    hist_CorrectMC, bins_CorrectMC = np.histogram( Correct_MC[:,0], bins=binnumber, range=(0, 1) )
    hist_realCC                    = ( hist/Confused_data[:,0].shape[0] - ccpercentage * hist_CorrectMC/Correct_MC[:,0].shape[0] ) / (1-ccpercentage)
    factor_realCC                  = 1 * Confused_MC[:,0].shape[0]
    g_Confused_data = TGraphErrors(len(center), np.array(center.tolist()), np.array(hist_realCC)* factor_realCC * factor_CorrectConfusedCountRatio, np.array([0]*len(center)), np.sqrt(hist_realCC*factor_realCC)*factor_CorrectConfusedCountRatio )

    TH_Correct_MC   .Sumw2();
    TH_Confused_MC  .Sumw2();

    TH_Correct_MC.SetFillColor(4)
    TH_Correct_MC.SetFillStyle(3004)
    TH_Correct_MC.SetLineColor(4)
    TH_Confused_MC.SetFillColor(2)
    TH_Confused_MC.SetFillStyle(3005)
    TH_Confused_MC.SetLineColor(2)

    g_Correct_data.SetMarkerStyle(15)
    g_Correct_data.SetMarkerColor(4)
    g_Correct_data.SetLineColor(4)
    g_Correct_data.SetMarkerSize(0.9)
    g_Confused_data.SetMarkerStyle(15)
    g_Confused_data.SetMarkerColor(2)
    g_Confused_data.SetLineColor(2)
    g_Confused_data.SetMarkerSize(0.9)

    TH_Correct_MC  .Draw("HIST")
    TH_Confused_MC .Draw("HIST same")
    g_Correct_data .Draw("E1 P")
    g_Confused_data.Draw("E1 P")

    leg =TLegend(0.35, 0.65, 0.55, 0.85)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)
    leg.SetTextFont(62)
    leg.AddEntry(TH_Correct_MC    , "Charge Correct (MC)")
    leg.AddEntry(TH_Confused_MC   , "Charge Confused (MC)")
    leg.AddEntry(g_Correct_data   , "Charge Correct (ISS)")
    leg.AddEntry(g_Confused_data  , "Charge Confused (ISS)")
    leg.Draw()

    xaxis = TH_Correct_MC.GetXaxis()
    yaxis = TH_Correct_MC.GetYaxis()
    xaxis.SetTitle('#Lambda_{CC}')
    yaxis.SetTitle('Normalized events')
    xaxis.SetTitleFont(62)
    xaxis.SetTitleSize(0.045) 
    xaxis.SetLabelFont(62)
    xaxis.SetLabelSize(0.05)
    yaxis.SetTitleFont(62)
    yaxis.SetTitleSize(0.045)
    yaxis.SetLabelFont(62)
    yaxis.SetLabelSize(0.05)

    yaxis.SetTitleOffset(1.1);
    #p1.SetLeftMargin(0.16);
    #p1.SetBottomMargin(0.2);

    gPad.SetLogy();
    c1.Update()
    c1.SaveAs(highpath + "/CCEstimatorplots/CorrectAndConfused_1D_Logy_" + str(Rigiditybin[i]) + "_binnumber_" + str(binnumber) + "_CCCut_" + str(CCCut) + "_CCN_" + str(CCN) + "_TRDN_" + str(TRDN) + "_InROOT.pdf")


def Plot_All(Confused_data, Correct_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, MarkerSize, factor_CorrectConfusedCountRatio, CCCut, CCN, TRDN):
    '''
    ## Correct 1D
    plot_Correct_1D(Correct_data, binnumber, Correct_MC, highpath, Rigiditybin, i, MarkerSize)

    ## Confused 1D
    plot_Confused_1D(Confused_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, MarkerSize, factor_CorrectConfusedCountRatio)

    ## Correct And Confused 1D
    plot_CorrectAndConfused_1D(Confused_data, Correct_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, MarkerSize, factor_CorrectConfusedCountRatio)
    '''

    ## Correct And Confused 1D (Plot in ROOT)
    plotRoot_CorrectAndConfused_1D(Confused_data, Correct_data, binnumber, Correct_MC, Confused_MC, ccpercentage, highpath, Rigiditybin, i, factor_CorrectConfusedCountRatio, CCCut, CCN, TRDN)

    '''
    #### 2D
    plt.figure(figsize=(30,18))
    plt.hist2d( Correct_MC[:,0], Correct_MC[:,1],  bins=100, norm=mpl.colors.LogNorm())
    plt.savefig(highpath + "/CCEstimatorplots/Correct_MC_2D_" + str(Rigiditybin[i]) + ".pdf")
    plt.close()

    plt.figure(figsize=(30,18))
    plt.hist2d( Confused_MC[:,0], Confused_MC[:,1],  bins=100, norm=mpl.colors.LogNorm())
    plt.savefig(highpath + "/CCEstimatorplots/Confused_MC_2D_" + str(Rigiditybin[i]) + ".pdf")
    plt.close()

    plt.figure(figsize=(30,18))
    plt.hist2d( Correct_data[:,0], Correct_data[:,1],  bins=100, norm=mpl.colors.LogNorm())
    plt.savefig(highpath + "/CCEstimatorplots/Correct_data_2D_" + str(Rigiditybin[i]) + ".pdf")
    plt.close()

    plt.figure(figsize=(30,18))
    plt.hist2d( Confused_data[:,0], Confused_data[:,1],  bins=100, norm=mpl.colors.LogNorm())
    plt.savefig(highpath + "/CCEstimatorplots/Confused_data_2D_" + str(Rigiditybin[i]) + ".pdf")
    plt.close()
    '''



def CalculateRejectionPower(binnumber, SignalArray, BackgroundArray):

    SignalAccumulatecount = 0
    FractionSignalSet     = np.array([])

    for SignalIterate in range(binnumber):
        SignalAccumulatecount = SignalAccumulatecount + np.histogram(SignalArray, binnumber, range=(0,1))[0][binnumber - SignalIterate - 1]
        FractionSignal        = SignalAccumulatecount / SignalArray.shape[0]
        FractionSignalSet     = np.append(FractionSignalSet, FractionSignal)

    BackgroundAccumulatecount = 0
    RejectionBackgroundSet     = np.array([])

    for BackgroundIterate in range(binnumber):
        BackgroundAccumulatecount = BackgroundAccumulatecount + np.histogram(BackgroundArray, binnumber, range=(0,1))[0][binnumber - BackgroundIterate - 1]   
        RejectionBackground       = 1.0/(BackgroundAccumulatecount/BackgroundArray.shape[0])
        RejectionBackgroundSet    = np.append(RejectionBackgroundSet, RejectionBackground) 

    return FractionSignalSet, RejectionBackgroundSet

    
def Plot_RejectionPower(FractionSignalSet, RejectionBackgroundSet, highpath, Rigiditybin, i, binnumber, CCCut, CCN, TRDN):

    c1 = TCanvas("", "", 800, 600)

    p1 = TPad("", "", 0.0, 0.0, 1.0, 1.0, 0) #(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1.Draw()
    p1.cd()

    gPad.SetFrameFillColor(0)
    gStyle.SetOptStat("00000000")

    g_RejectionPower  = TGraphErrors(len(FractionSignalSet), FractionSignalSet, RejectionBackgroundSet, np.array([0]*len(FractionSignalSet)), np.array([0]*len(FractionSignalSet)))

    g_RejectionPower.SetMarkerStyle(15)
    g_RejectionPower.SetMarkerColor(4)
    g_RejectionPower.SetLineColor(4)
    g_RejectionPower.SetMarkerSize(0.9)

    g_RejectionPower.Draw("")

    g_RejectionPower.SetTitle("")

    xaxis = g_RejectionPower.GetXaxis()
    yaxis = g_RejectionPower.GetYaxis()
    xaxis.SetTitle('Signal efficiency')
    yaxis.SetTitle('Background rejection power')
    xaxis.SetTitleFont(62)
    xaxis.SetTitleSize(0.045) 
    xaxis.SetLabelFont(62)
    xaxis.SetLabelSize(0.05)
    yaxis.SetTitleFont(62)
    yaxis.SetTitleSize(0.045)
    yaxis.SetLabelFont(62)
    yaxis.SetLabelSize(0.05)

    yaxis.SetTitleOffset(1.1);

    #p1.SetBottomMargin(0.17);
    #p1.SetLeftMargin(0.14);

    gPad.SetLogy();
    c1.Update()
    c1.SaveAs(highpath + "/CCEstimatorplots/RejectionPower" + str(Rigiditybin[i]) + "_binnumber_" + str(binnumber) + "_CCCut_" + str(CCCut) + "_CCN_" + str(CCN) + "_TRDN_" + str(TRDN) + "_InROOT.pdf")





