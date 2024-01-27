
void Plot_MeasuringTime_OnlyStormer(TH1D *GEOMETRIC_time){
    TCanvas *plot = new TCanvas();

    GEOMETRIC_time->Draw("");

    GEOMETRIC_time->SetTitle("");
    GEOMETRIC_time->SetLineColor(4);
    GEOMETRIC_time->SetLineWidth(2);

    gPad  ->SetLogx();
    gStyle->SetOptStat("00000000");

    TAxis * xaxis = GEOMETRIC_time->GetXaxis();
    TAxis * yaxis = GEOMETRIC_time->GetYaxis();
    xaxis->SetRangeUser(0, 1000);
    yaxis->SetRangeUser(0, 314150400);

    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("T / (s)");
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);

    xaxis->SetLabelOffset(0);
    yaxis->SetLabelOffset(0.02);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetTitleOffset(0);
    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.13);

    // Draw Text
    TText *t = new TText(0.28, 0.57, "Geomagnetic Cutoff");
    t->SetNDC();
    t->SetTextFont(62);
    t->SetTextSize(0.03);
    t->Draw();

    TText *t2 = new TText(0.88, 0.655, "2488 days");
    t2->SetNDC();
    t2->SetTextFont(62);
    t2->SetTextSize(0.04);
    t2->Draw();

    TText *t3 = new TText(0.88, 0.88, "3636 days");
    t3->SetNDC();
    t3->SetTextFont(62);
    t3->SetTextSize(0.04);
    t3->Draw();

    TText *t4 = new TText(0.6, 0.78, "Detector Quality Cuts and");
    t4->SetNDC();
    t4->SetTextFont(62);
    t4->SetTextSize(0.03);
    t4->Draw();

    TText *t5 = new TText(0.63, 0.74, "Live Time Fraction");
    t5->SetNDC();
    t5->SetTextFont(62);
    t5->SetTextSize(0.03);
    t5->Draw();

    // Draw Arrow
    TArrow *ar2 = new TArrow(50, 314150400, 50, 214979330, 0.02, "|>");
    ar2->Draw();

    TArrow *ar3 = new TArrow(2, 214979330, 2, 20000000, 0.02, "|>");
    ar3->Draw();

    // Draw Line
    TLine *l = new TLine(0, 214979330, 1000, 214979330);
    l->Draw();

    gPad->Update();

    plot->SaveAs("MeasuringTime.pdf");

}


void Plot_MeasuringTime(TH1D *GEOMETRIC_time, TH1D *IGRF_time){

    TCanvas *plot = new TCanvas();

    GEOMETRIC_time->Draw("");
    IGRF_time     ->Draw("same");

    GEOMETRIC_time->SetLineColor(4);
    GEOMETRIC_time->SetLineWidth(2);
    IGRF_time     ->SetLineColor(2);
    IGRF_time     ->SetLineWidth(2);

    gPad  ->SetLogx();
    gStyle->SetOptStat("00000000");

    TAxis * xaxis = GEOMETRIC_time->GetXaxis();
    TAxis * yaxis = GEOMETRIC_time->GetYaxis();
    xaxis->SetMoreLogLabels();
    xaxis->SetRangeUser(0, 100.0);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleFont(62);
    yaxis->SetTitleSize(0.045);
    xaxis->SetNoExponent();
    //xaxis->SetExponentOffset(Float_t xoff, Float_t yoff, Option_t *axis);
    //GEOMETRIC_time->SetExponentOffset(0, 1);
    //TGaxis::SetExponentOffset(-0.05, 0.9, "x");        
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("T / (s)");

    xaxis->SetLabelOffset(0);
    yaxis->SetLabelOffset(0.02);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetTitleOffset(0);

    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.16);

    GEOMETRIC_time->SetLineColor(4);
    IGRF_time     ->SetLineColor(2);
    GEOMETRIC_time->SetTitle("");

    TLegend *leg = new TLegend(0.30, 0.65, 0.5, 0.85);
    gStyle->SetLegendTextSize(0.05);
    leg->SetBorderSize(0);
    leg->SetTextFont(62);
    //leg->AddEntry(GEOMETRIC_time, "St#phimer Cutoff", "lpf"); //"St#oslashmer Cutoff"
    leg->AddEntry(GEOMETRIC_time, "St  rmer Cutoff", "lpf");
    leg->AddEntry(IGRF_time     , "IGRF Cutoff", "lpf");
    /*
    leg->AddEntry(IGRF_time       , "\u00F8"     , "lpf");
    leg->AddEntry(GEOMETRIC_time, "Størmer &oslash Cutoff","lpf");
    leg->AddEntry(GEOMETRIC_time, "#frac{Størmern}{&oslash} '\u00F8' \u0444 $\phi$ ");
    leg->AddEntry(GEOMETRIC_time, "St#ΦrmerA_{p}/A_{#bar{p}}","lpf");
    leg->AddEntry(GEOMETRIC_time, (string("St") + string("o") + string("mer")).c_str(), "lpf");
    leg->AddEntry(GEOMETRIC_time, (string("St") + #oslash + string("mer")).c_str(), "lpf");
    leg->AddEntry(GEOMETRIC_time, "St#phimer", "lpf");
    leg->AddEntry(GEOMETRIC_time, "\u00F8", "lpf");
    leg->AddEntry(GEOMETRIC_time, "\upphi", "lpf");
    yaxis_residuallow->SetTitle("#frac{This - PRL}{PRL}");
    legend_accerr->AddEntry(Effective_Acceptance_ratio, "A_{p}/A_{#bar{p}}","lpf");
    */
    leg->Draw();

    gPad->Update();

    TLatex t;
    t.SetNDC(1);
    t.SetTextAngle(0);
    t.SetTextColor(kBlack);
    t.SetTextFont(62);
    t.SetTextSize(0.05);
    //t.DrawLatex(0.3 , 0.5, "St  rmer");
    t.SetTextAngle(-10);
    t.SetTextFont(122);
    t.SetTextSize(0.055);
    //t.DrawLatex(0.327 , 0.502, "#phi");
    t.DrawLatex(0.377 , 0.782, "#phi");
    
    plot->SaveAs("IGRFvsStomer.pdf");
}

void Plot_MeasuringTimeVsLiveTime(TH1D *GEOMETRIC_fMeasuringTimeVsLiveTime){
    TCanvas *plot = new TCanvas();

    GEOMETRIC_fMeasuringTimeVsLiveTime->Draw("HIST");

    GEOMETRIC_fMeasuringTimeVsLiveTime->SetLineColor(1);
    GEOMETRIC_fMeasuringTimeVsLiveTime->SetLineWidth(2);

    GEOMETRIC_fMeasuringTimeVsLiveTime->SetTitle("");
    gStyle->SetOptStat(0);
    gPad  ->SetLogy();

    TAxis * xaxis = GEOMETRIC_fMeasuringTimeVsLiveTime->GetXaxis();
    TAxis * yaxis = GEOMETRIC_fMeasuringTimeVsLiveTime->GetYaxis();
    xaxis->SetRangeUser(0.48 , 1.0);
    yaxis->SetRangeUser(10000, 10000000);

    xaxis->SetTitle("Live time fraction");
    yaxis->SetTitle("T / (s)");
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);

    xaxis->SetLabelOffset(0);
    yaxis->SetLabelOffset(0.02);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetTitleOffset(0);

    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.16);

    gPad->Update();

    plot->SaveAs("MeasuringTimeVsLiveTime.pdf");
}

void Plot_ParticlesVsTriggers(TH2D *GEOMETRIC_fParticlesVsTriggers){
    TCanvas *plot = new TCanvas();

    GEOMETRIC_fParticlesVsTriggers->Draw("");

    GEOMETRIC_fParticlesVsTriggers->SetTitle("");
    gStyle->SetOptStat(0);

    TAxis * xaxis = GEOMETRIC_fParticlesVsTriggers->GetXaxis();
    TAxis * yaxis = GEOMETRIC_fParticlesVsTriggers->GetYaxis();
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.06);
    yaxis->SetTitleSize(0.06);

    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);

    xaxis->SetLabelOffset(0);
    yaxis->SetLabelOffset(0.02);
    xaxis->SetTitleOffset(1.2);
    yaxis->SetTitleOffset(0);

    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.16);

    gPad->Update();
    plot->SaveAs("ParticlesVsTriggers.pdf");
}


void Plot_TriggerRateVsPosition(TH3F *GEOMETRIC_fTriggerRateVsPosition){
    TCanvas *plot = new TCanvas();

    GEOMETRIC_fTriggerRateVsPosition->Draw();

    GEOMETRIC_fTriggerRateVsPosition->SetTitle("");
    gStyle->SetOptStat(0);

    TAxis * xaxis = GEOMETRIC_fTriggerRateVsPosition->GetXaxis();
    TAxis * yaxis = GEOMETRIC_fTriggerRateVsPosition->GetYaxis();
    TAxis * zaxis = GEOMETRIC_fTriggerRateVsPosition->GetZaxis();
    xaxis->SetTitleOffset(1.7);
    yaxis->SetTitleOffset(1.9);
    zaxis->SetTitleOffset(1.4);

    plot->SaveAs("TriggerRateVsPosition.pdf");
}

void Plot_CutOffRigidityVsISSPosition(TH3F *GEOMETRIC_fCutOffRigidityVsISSPosition){

    TCanvas *plot = new TCanvas();

    TProfile2D * XYProjection = GEOMETRIC_fCutOffRigidityVsISSPosition->Project3DProfile("yx");
    XYProjection->Draw("COLZ");

    XYProjection->SetTitle("");
    gStyle->SetOptStat(0);

    TAxis * xaxis = XYProjection->GetXaxis();
    TAxis * yaxis = XYProjection->GetYaxis();
    TAxis * zaxis = XYProjection->GetZaxis();
    xaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelFont(62);
    yaxis->SetLabelSize(0.05);
    zaxis->SetLabelFont(62);
    zaxis->SetLabelSize(0.05);
    xaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleFont(62);
    yaxis->SetTitleSize(0.045);
    zaxis->SetTitleFont(62);
    zaxis->SetTitleSize(0.045);

    TLatex t;
    t.SetNDC(1);
    t.SetTextAlign(22);
    t.SetTextColor(kBlack);
    t.SetTextFont(62);
    t.SetTextSize(0.05);
    t.DrawLatex(0.5 , 0.02, "Longitude / (#circ)");
    t.DrawLatex(0.25, 0.02, "West");
    t.DrawLatex(0.75, 0.02, "East");
    t.DrawLatex(0.91, 0.94, "|R| / (GV)");
    t.SetTextAngle(90);
    t.DrawLatex(0.03, 0.5 , "Latitude / (#circ)");
    t.DrawLatex(0.03, 0.75, "North");
    t.DrawLatex(0.03, 0.25, "South");      


    plot->SaveAs("CutOffRigidityVsISSPosition.pdf");
}

void Plot_LiveTimeVsISSPosition(TH3F *GEOMETRIC_fLiveTimeVsISSPosition){

    TCanvas *plot = new TCanvas();

    TProfile2D * XYProjection = GEOMETRIC_fLiveTimeVsISSPosition->Project3DProfile("yx");

    //GEOMETRIC_fLiveTimeVsISSPosition->Draw("");
    XYProjection->Draw("COLZ");

    //GEOMETRIC_fLiveTimeVsISSPosition->SetTitle("");
    //GEOMETRIC_fLiveTimeVsISSPosition->SetMaximum(1.0);
    //GEOMETRIC_fLiveTimeVsISSPosition->SetMinimum(0.0);
    XYProjection->SetTitle("");

    gStyle->SetOptStat(0);

    //TAxis * xaxis = GEOMETRIC_fLiveTimeVsISSPosition->GetXaxis();
    //TAxis * yaxis = GEOMETRIC_fLiveTimeVsISSPosition->GetYaxis();
    //TAxis * zaxis = GEOMETRIC_fLiveTimeVsISSPosition->GetZaxis();
    //xaxis->SetTitleOffset(1.7);
    //yaxis->SetTitleOffset(1.9);
    //zaxis->SetTitle("Live time fraction");

    TAxis * xaxis = XYProjection->GetXaxis();
    TAxis * yaxis = XYProjection->GetYaxis();
    TAxis * zaxis = XYProjection->GetZaxis();
    //xaxis->SetTitle("Longitude (degree)");
    //yaxis->SetTitle("Latitude (#circ)");
    xaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelFont(62);
    yaxis->SetLabelSize(0.05);
    zaxis->SetLabelFont(62);
    zaxis->SetLabelSize(0.05);
    xaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleFont(62);
    yaxis->SetTitleSize(0.045);
    zaxis->SetTitleFont(62);
    zaxis->SetTitleSize(0.045);

    TLatex t;
    t.SetNDC(1);
    t.SetTextAlign(22);
    t.SetTextColor(kBlack);
    t.SetTextFont(62);
    t.SetTextSize(0.05);
    t.DrawLatex(0.5 , 0.02, "Longitude / (#circ)");
    t.DrawLatex(0.25, 0.02, "West");   
    t.DrawLatex(0.75, 0.02, "East");
    t.DrawLatex(0.88, 0.94, "Live time fraction"); 
    t.SetTextAngle(90);
    t.DrawLatex(0.03, 0.5 , "Latitude / (#circ)");
    t.DrawLatex(0.03, 0.75, "North");
    t.DrawLatex(0.03, 0.25, "South");

    plot->SaveAs("LiveTimeFractionVsISSPosition.pdf");
}









