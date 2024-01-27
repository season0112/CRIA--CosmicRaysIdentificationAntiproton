

void Plot_EffectiveAcceptance(TGraphAsymmErrors *Acceptance_low, TGraphAsymmErrors *g_EffectiveAcceptance_high_P0, TGraphAsymmErrors *g_EffectiveAcceptance_high_P1, TGraphAsymmErrors *g_EffectiveAcceptance_high_P2, TGraphAsymmErrors *g_EffectiveAcceptance_high_P4, string specie){

    TCanvas c1("c1","c1",1000,1000);

    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0, 0);
    p1->Draw();
    p1->cd();

    Acceptance_low    ->Draw("AP");

    g_EffectiveAcceptance_high_P0         ->Draw("same P");
    g_EffectiveAcceptance_high_P1         ->Draw("same P");
    g_EffectiveAcceptance_high_P2         ->Draw("same P");
    g_EffectiveAcceptance_high_P4         ->Draw("same P");

    Acceptance_low->SetTitle("");

    Acceptance_low->SetMarkerStyle(15);
    Acceptance_low->SetMarkerColor(2);
    Acceptance_low->SetMarkerSize(0.9);
    g_EffectiveAcceptance_high_P0->SetMarkerStyle(15);
    g_EffectiveAcceptance_high_P0->SetMarkerColor(8);
    g_EffectiveAcceptance_high_P0->SetMarkerSize(0.9);
    g_EffectiveAcceptance_high_P1->SetMarkerStyle(15);
    g_EffectiveAcceptance_high_P1->SetMarkerColor(9);
    g_EffectiveAcceptance_high_P1->SetMarkerSize(0.9);
    g_EffectiveAcceptance_high_P2->SetMarkerStyle(15);
    g_EffectiveAcceptance_high_P2->SetMarkerColor(12);
    g_EffectiveAcceptance_high_P2->SetMarkerSize(0.9);
    g_EffectiveAcceptance_high_P4->SetMarkerStyle(15);
    g_EffectiveAcceptance_high_P4->SetMarkerColor(30);
    g_EffectiveAcceptance_high_P4->SetMarkerSize(0.9);

    TAxis * xaxis = Acceptance_low->GetXaxis();
    TAxis * yaxis = Acceptance_low->GetYaxis();
    xaxis->SetTitle("|R|_{true} / (GV)");
    xaxis->SetMoreLogLabels();
    double ymax;
    if (specie == std::string("Proton")){
        yaxis->SetTitle("A_{p} / (cm^{2} sr)");
        ymax = 2500;
    }
    else if (specie == std::string("Antiproton")){
        yaxis->SetTitle("A_{p} / (cm^{2} sr)");

        TLatex latex;
        latex.SetNDC(1);
        latex.SetTextSize(0.055);
        latex.SetTextAngle(90);
        latex.SetTextAlign(13);  //align at top
        latex.DrawLatex(0.054, 0.718, "#minus");
        ymax = 2500;
    }
    TLine line1(15.3, 0, 15.3, ymax);
    line1.Draw("same");

    xaxis->SetLimits(1.0,600);
    yaxis->SetRangeUser(0, ymax);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.04);
    yaxis->SetTitleSize(0.04);
    xaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.04);
    yaxis->SetLabelFont(62);
    yaxis->SetLabelSize(0.04);

    //xaxis->SetLabelOffset(-0.01);
    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(2.0);
    p1->SetBottomMargin(0.2);
    p1->SetLeftMargin(0.2);

    gPad->SetLogx();

    TLegend *legend = new TLegend(0.22, 0.25, 0.37, 0.50);
    legend->AddEntry(Acceptance_low               , "InnerCentral", "p");
    legend->AddEntry(g_EffectiveAcceptance_high_P0, "FullSpan"    , "p");
    legend->AddEntry(g_EffectiveAcceptance_high_P1, "Inner + L1"  , "p");
    legend->AddEntry(g_EffectiveAcceptance_high_P2, "Inner + L9"  , "p");
    legend->AddEntry(g_EffectiveAcceptance_high_P4, "Inner Only"  , "p");
    legend->SetTextSize(0.04);
    legend->SetTextFont(62);
    legend->SetBorderSize(0);
    legend->Draw();

    c1.SaveAs( (string("Acceptance_") + specie + string(".pdf")).c_str());

}





