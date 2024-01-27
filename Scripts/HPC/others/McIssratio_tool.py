from ROOT import TFile, TH1D, TCanvas, gPad, gStyle, TPad, TLegend

def plot(TH_ratio, TH_iss, TH_mc, mcdataset):

    c1 = TCanvas()

    gStyle.SetOptStat(0)

    p2 = TPad("p2","p2",0.,0.,1.,0.3)
    p2.Draw()
    p2.SetTopMargin(0.001)
    p2.SetBottomMargin(0.3)
    p2.SetLogx()
    p2.SetLogy()

    p1 = TPad("p1","p1",0.,0.3,1.,1.0)
    p1.Draw()
    p1.SetBottomMargin(0.001)
    p1.cd()
    p1.SetLogx()
    p1.SetLogy()

    # Plot Event Numbers
    TH_mc.Draw("")
    TH_iss.Draw("SAME")

    TH_mc .SetLineColor(4)
    TH_iss.SetLineColor(2)
    TH_mc .SetLineWidth(3)
    TH_iss.SetLineWidth(3)

    xaxis_up = TH_mc.GetXaxis()
    yaxis_up = TH_mc.GetYaxis()
    if mcdataset == "B1042_pr.pl1.flux.l1a9.2016000_7.6_all":
        xaxis_up.SetRange(22,30)
    xaxis_up.SetTitleFont(62);
    yaxis_up.SetTitleFont(62);
    xaxis_up.SetTitleSize(0.05);
    yaxis_up.SetTitleSize(0.05);
    xaxis_up.SetLabelFont(62);
    yaxis_up.SetLabelFont(62);
    xaxis_up.SetLabelSize(0.05);
    yaxis_up.SetLabelSize(0.05);
    #xaxis_up.SetTickSize(0.1)
    yaxis_up.SetTitle("N / #Delta R")
    yaxis_up.SetRangeUser(0.05, 4000000)

    leg = TLegend(0.75, 0.7, 0.89, 0.85)
    leg.AddEntry(TH_mc , "MC" , "l")
    leg.AddEntry(TH_iss, "ISS", "l")
    leg.SetTextSize(0.04);
    leg.SetTextFont(62);
    leg.SetBorderSize(0);
    leg.Draw()

    # Plot ratio
    p2.cd()

    TH_ratio.Draw("")

    TH_ratio.SetLineWidth(3)
    TH_ratio.SetLineColor(1)

    xaxis_down = TH_ratio.GetXaxis()
    yaxis_down = TH_ratio.GetYaxis()
    xaxis_down.SetTitle("|R| / (GV)")
    yaxis_down.SetTitle("MC / ISS")

    xaxis_down.SetMoreLogLabels()

    xaxis_down.SetTitleFont(62);
    yaxis_down.SetTitleFont(62);
    xaxis_down.SetTitleSize(0.118);
    yaxis_down.SetTitleSize(0.118);
    xaxis_down.SetLabelFont(62);
    yaxis_down.SetLabelFont(62);
    xaxis_down.SetLabelSize(0.12);
    yaxis_down.SetLabelSize(0.12);
    #xaxis_down.SetTickSize(0.1)
    xaxis_down.SetTitleOffset(1.29);
    yaxis_down.SetTitleOffset(0.45);

    if mcdataset == "B1042_pr.pl1.flux.l1a9.2016000_7.6_all":
        xaxis_down.SetRange(22,30)
        yaxis_down.SetMoreLogLabels()

    c1.SaveAs(str(mcdataset)+"_mc_iss_ratio.pdf")




