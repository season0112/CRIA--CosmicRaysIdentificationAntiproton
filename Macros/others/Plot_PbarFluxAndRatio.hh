
void RemoveXaxis(TGraphAsymmErrors *result){
    for (int i=0; i < result->GetN(); i++){
        result->SetPointEXhigh(i,0);
        result->SetPointEXlow (i,0);
    }
}

void Plot_PbarOverProton_CompareInLowEnergy(TGraphAsymmErrors *ratioPAMELA, TGraphErrors *ratio_pass78){

    TCanvas c4("c4","c4",1000,500);

    ratioPAMELA->SetLineColor(4);
    ratioPAMELA->SetMarkerColor(4);
    ratioPAMELA->SetMarkerStyle(16);
    ratioPAMELA->SetMarkerSize(0.8);
    ratio_pass78->SetLineColor(2);
    ratio_pass78->SetMarkerColor(2);
    ratio_pass78->SetMarkerStyle(18);
    ratio_pass78->SetMarkerSize(0.8);

    ratioPAMELA->Draw("AP");
    ratio_pass78->Draw("same P");

    //ratioPAMELA->GetXaxis()->SetRangeUser(0.5, 5);
    ratioPAMELA->GetXaxis()->SetLimits(0, 5);
    ratioPAMELA->GetYaxis()->SetRangeUser(0, 0.00014);    

    ratioPAMELA->GetXaxis()->SetTitle("Rigidity (GV)");
    ratioPAMELA->GetYaxis()->SetTitle("Antiproton to Proton Ratio");
    ratioPAMELA->SetTitle("");
    ratioPAMELA->GetXaxis()->SetTitleFont(62);
    ratioPAMELA->GetYaxis()->SetTitleFont(62);
    ratioPAMELA->GetXaxis()->SetLabelFont(62);
    ratioPAMELA->GetXaxis()->SetLabelSize(0.05);
    ratioPAMELA->GetYaxis()->SetLabelFont(62);
    ratioPAMELA->GetYaxis()->SetLabelSize(0.05);

    gPad->SetLeftMargin(0.15);

    TLegend *legend4 = new TLegend(0.25,0.68,0.70,0.88);
    legend4->AddEntry(ratioPAMELA,"PAMELA(2006/06-2010/01)","lpf");
    legend4->AddEntry(ratio_pass78,"AMS(2011/05-2020/06)(Statistical Error Only)","lpf");
    legend4->SetTextFont(62);
    legend4->SetTextSize(0.027);
    legend4->Draw();
    c4.SaveAs("PbarOverProton_LowRatioCompare.pdf");

}


void Plot_PbarOverProton_Compare_above10(TGraphAsymmErrors *ratioCAPRICE98, TGraphAsymmErrors *ratioHEATpbar, TGraphAsymmErrors *ratioPAMELA3){

    TCanvas *c3 = new TCanvas;
    gPad->SetLogx();
    gPad->SetLogy();

    ratioCAPRICE98->GetXaxis()->SetRangeUser(8, 525);
    ratioCAPRICE98->GetXaxis()->SetLimits(8, 525);
    ratioCAPRICE98->SetLineColor(2);
    ratioCAPRICE98->SetMarkerColor(2);

    ratioHEATpbar->SetLineColor(4);
    ratioHEATpbar->SetMarkerColor(4);

    ratioCAPRICE98->Draw("AP*");
    ratioHEATpbar->Draw("same *");
    ratioPAMELA3->Draw("same *");

    ratioCAPRICE98->SetTitle("");
    ratioCAPRICE98->GetXaxis()->SetMoreLogLabels();
    ratioCAPRICE98->GetYaxis()->SetRangeUser(0.00003, 0.004);
    ratioCAPRICE98->GetXaxis()->SetRangeUser(8, 525);
    ratioCAPRICE98->GetXaxis()->SetLimits(8, 525);
    ratioCAPRICE98->GetXaxis()->SetTitle("Ekn (GeV/n)");
    ratioCAPRICE98->GetYaxis()->SetTitle("#bar{p}/p");
    ratioCAPRICE98->GetXaxis()->SetTitleSize(0.05);
    ratioCAPRICE98->GetYaxis()->SetTitleSize(0.05);

    TLegend *legend1 = new TLegend(0.45,0.68,0.78,0.88);
    legend1->AddEntry(ratioCAPRICE98,"CAPRICE98","lpf");
    legend1->AddEntry(ratioHEATpbar,"HEAT-pbar","lpf");
    legend1->AddEntry(ratioPAMELA3,"PAMELA(2006/07-2009/12)","lpf");
    legend1->Draw();

    c3->SaveAs("PbarOverProton_Compare_above10.pdf");
}


void Plot_PbarFluxMulR_Compare_above10(TGraphAsymmErrors *ams, TGraphAsymmErrors *CAPRICE98, TGraphAsymmErrors *MASS91, TGraphAsymmErrors *PAMELA2){

    TCanvas *c = new TCanvas;
    gPad->SetLogx();
    gPad->SetLogy();

    ams->GetXaxis()->SetRangeUser(8, 525);
    ams->GetXaxis()->SetLimits(8, 525);

    ams      ->Draw("AP*");
    CAPRICE98->Draw("same *");
    MASS91   ->Draw("same *");
    PAMELA2  ->Draw("same *");

    ams->SetTitle("");
    ams->GetXaxis()->SetMoreLogLabels();
    //ams->GetYaxis()->SetMoreLogLabels();
    ams->GetYaxis()->SetRangeUser(0.4, 14);
    ams->GetXaxis()->SetRangeUser(8, 525);
    ams->GetXaxis()->SetLimits(8, 525);
    ams->GetXaxis()->SetTitle("|R| / (GV)"); //???
    ams->GetYaxis()->SetTitle("#Phi#upointR^{2.7}");
    ams->GetXaxis()->SetTitleSize(0.05);
    ams->GetYaxis()->SetTitleSize(0.05);
    ams->GetXaxis()->SetLabelSize(0.05);
    ams->GetYaxis()->SetLabelSize(0.05);

    TLegend *legend2 = new TLegend(0.45, 0.68, 0.78, 0.88);
    legend2->AddEntry(ams      , "AMS02"                  , "lpf");
    legend2->AddEntry(CAPRICE98, "CAPRICE98"              , "lpf");
    legend2->AddEntry(MASS91   , "MASS91"                 , "lpf");
    legend2->AddEntry(PAMELA2  , "PAMELA(2006/07-2009/12)", "lpf");
    legend2->Draw();

    c->SaveAs("Pbar_Compare_above10.pdf");
}



void Plot_PbarFlux_CRDB_Compare(TGraphAsymmErrors *AMS02, TGraphAsymmErrors *BESSPolarII, TGraphAsymmErrors *PAMELA, TGraphAsymmErrors *BESSPolarI, TGraphAsymmErrors *BESSTeV, TGraphAsymmErrors *BESS00, TGraphAsymmErrors *BESS93, TGraphAsymmErrors *BESS95, TGraphAsymmErrors *BESS97, TGraphAsymmErrors *BESS98, TGraphAsymmErrors *BESS99, TGraphAsymmErrors *CAPRICE94, TGraphAsymmErrors *CAPRICE98, TGraphAsymmErrors *IMAX92, TGraphAsymmErrors *MASS91){

    TCanvas *c3 = new TCanvas("c3","c3",800,600);

    TPad *p1 = new TPad("p1"   , "p1"   , 0.05, 0.0  , 1.0, 1.0, 0); //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1->Draw();
    p1->cd();

    gPad->SetLogx();
    gPad->SetLogy();

    AMS02      ->Draw("AP E ");
    BESSPolarII->Draw("same P");  // ->Draw("same P X0");
    PAMELA     ->Draw("same P");

    BESSPolarI ->Draw("same P");
    BESSTeV    ->Draw("same P");
    BESS00     ->Draw("same P");
    BESS93     ->Draw("same P");
    BESS95     ->Draw("same P");
    BESS97     ->Draw("same P");
    BESS98     ->Draw("same P");
    BESS99     ->Draw("same P");
    CAPRICE94  ->Draw("same P");
    CAPRICE98  ->Draw("same P");
    IMAX92     ->Draw("same P");
    MASS91     ->Draw("same P");
    //AMS02->Draw("A P");

    AMS02->GetYaxis()->SetRangeUser(0.0000001, 0.2);
    AMS02->GetXaxis()->SetLimits(0.5,700);
    AMS02->SetTitle("");

    gStyle->SetErrorX(0);
    c3->Update();

    AMS02      ->SetMarkerColor(1);
    AMS02      ->SetLineColor(1);
    PAMELA     ->SetMarkerColor(2);
    PAMELA     ->SetLineColor(2);
    BESSPolarII->SetMarkerColor(4);
    BESSPolarII->SetLineColor(4);
    BESSPolarI ->SetMarkerColor(3);
    BESSPolarI ->SetLineColor(3);
    BESSTeV    ->SetMarkerColor(91);
    BESSTeV    ->SetLineColor(91);
    BESS00     ->SetMarkerColor(6);
    BESS00     ->SetLineColor(6);
    BESS93     ->SetMarkerColor(7);
    BESS93     ->SetLineColor(7);
    BESS95     ->SetMarkerColor(8);
    BESS95     ->SetLineColor(8);
    BESS97     ->SetMarkerColor(9);
    BESS97     ->SetLineColor(9);
    BESS98     ->SetMarkerColor(11);
    BESS98     ->SetLineColor(11);
    BESS99     ->SetMarkerColor(38);
    BESS99     ->SetLineColor(38);
    CAPRICE94  ->SetMarkerColor(41);
    CAPRICE94  ->SetLineColor(41);
    CAPRICE98  ->SetMarkerColor(46);
    CAPRICE98  ->SetLineColor(46);
    IMAX92     ->SetMarkerColor(52);
    IMAX92     ->SetLineColor(52);
    MASS91     ->SetMarkerColor(61);
    MASS91     ->SetLineColor(61);


    AMS02->SetMarkerStyle(8);
    AMS02->SetMarkerSize(1);
    PAMELA->SetMarkerStyle(8);
    PAMELA->SetMarkerSize(1);
    BESSPolarII->SetMarkerStyle(8);
    BESSPolarII->SetMarkerSize(1);
    BESSPolarI->SetMarkerStyle(8);
    BESSPolarI->SetMarkerSize(1);
    BESSTeV->SetMarkerStyle(8);
    BESSTeV->SetMarkerSize(1);
    BESS00->SetMarkerStyle(8);
    BESS00->SetMarkerSize(1);
    BESS93->SetMarkerStyle(8);
    BESS93->SetMarkerSize(1);
    BESS95->SetMarkerStyle(8);
    BESS95->SetMarkerSize(1);
    BESS97->SetMarkerStyle(8);
    BESS97->SetMarkerSize(1);
    BESS98->SetMarkerStyle(8);
    BESS98->SetMarkerSize(1);
    BESS99->SetMarkerStyle(8);
    BESS99->SetMarkerSize(1);
    CAPRICE94->SetMarkerStyle(8);
    CAPRICE94->SetMarkerSize(1);
    CAPRICE98->SetMarkerStyle(8);
    CAPRICE98->SetMarkerSize(1);
    IMAX92->SetMarkerStyle(8);
    IMAX92->SetMarkerSize(1);
    MASS91->SetMarkerStyle(8);
    MASS91->SetMarkerSize(1);

    TAxis * xaxis = AMS02->GetXaxis();
    TAxis * yaxis = AMS02->GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#Phi_{#bar{p}} (GV m^{2} s sr)^{-1}");

    xaxis->SetTitleSize(0.045);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleSize(0.045);
    yaxis->SetTitleFont(62);

    xaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelFont(62);
    yaxis->SetLabelSize(0.05);
    //xaxis->SetMoreLogLabels();
    xaxis->SetTitleOffset(1.3);
    //yaxis->SetTitleOffset(1.);

    //gPad->SetTopMargin(0.2);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.3);

    TLegend *legend1 = new TLegend(0.72, 0.3, 1.0, 0.9); //TLegend(0.2, 0.18, 0.50, 0.38)
    TLegendEntry *legAMS02       = legend1->AddEntry(AMS02      , "AMS-02 (2011/05-2015/05)"      , "lp");
    TLegendEntry *legBESSPolarII = legend1->AddEntry(BESSPolarII, "BESS-PolarII (2007/12-2008/01)", "lp");
    TLegendEntry *legPAMELA      = legend1->AddEntry(PAMELA     , "PAMELA (2006/07-2009/12)"      , "lp");
    TLegendEntry *legBESSPolarI  = legend1->AddEntry(BESSPolarI , "BESS-PolarI (2004/12)"          , "lp");
    TLegendEntry *legBESSTeV     = legend1->AddEntry(BESSTeV    , "BESS-TeV (2002/08)"             , "lp");
    TLegendEntry *legBESS00      = legend1->AddEntry(BESS00     , "BESS00 (2000/08)"              , "lp");
    TLegendEntry *legBESS93      = legend1->AddEntry(BESS93     , "BESS93 (1993/07)"              , "lp");
    TLegendEntry *legBESS95      = legend1->AddEntry(BESS95     , "BESS95 (1995/07)"              , "lp");
    TLegendEntry *legBESS97      = legend1->AddEntry(BESS97     , "BESS97 (1997/07)"              , "lp");
    TLegendEntry *legBESS98      = legend1->AddEntry(BESS98     , "BESS98 (1998/07)"              , "lp");
    TLegendEntry *legBESS99      = legend1->AddEntry(BESS99     , "BESS99 (1999/08)"              , "lp");
    TLegendEntry *legCAPRICE94   = legend1->AddEntry(CAPRICE94  , "CAPRICE94 (1994/08)"           , "lp");
    TLegendEntry *legCAPRICE98   = legend1->AddEntry(CAPRICE98  , "CAPRICE98 (1998/05)"           , "lp");
    TLegendEntry *legIMAX92      = legend1->AddEntry(IMAX92     , "IMAX92 (1992/07)"              , "lp");
    TLegendEntry *legMASS91      = legend1->AddEntry(MASS91     , "MASS91 (1991/09)"              , "lp");

    legAMS02      ->SetTextColor(1);
    legPAMELA     ->SetTextColor(2);
    legBESSPolarII->SetTextColor(4);
    legBESSPolarI ->SetTextColor(3);
    legBESSTeV    ->SetTextColor(91);
    legBESS00     ->SetTextColor(6);
    legBESS93     ->SetTextColor(7);
    legBESS95     ->SetTextColor(8);
    legBESS97     ->SetTextColor(9);
    legBESS98     ->SetTextColor(11);
    legBESS99     ->SetTextColor(38);
    legCAPRICE94  ->SetTextColor(41);
    legCAPRICE98  ->SetTextColor(46);
    legIMAX92     ->SetTextColor(52);
    legMASS91     ->SetTextColor(61);

    legend1->SetBorderSize(0);
    legend1->SetTextFont(62);
    legend1->SetTextSize(0.02);
    legend1->Draw();


    c3->Update();

    c3->SaveAs("PbarFluxForDifferentExperiments.pdf");

}























