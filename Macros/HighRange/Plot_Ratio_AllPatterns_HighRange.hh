
void MixedPattern_Selection(
int Range1_FirstIndex, int Range1_LastIndex, int Range2_FirstIndex, int Range2_LastIndex, int Range3_FirstIndex, int Range3_LastIndex, int Range4_FirstIndex, int Range4_LastIndex, 

TH1D *h_Antiproton_number_unfolded_High_mixedPattern, TH1D *h_Proton_number_unfolded_High_mixedPattern, TH1D *h_Antiproton_number_unfolded_High_P0, TH1D *h_Antiproton_number_unfolded_High_P0VGGNN, TH1D *h_Antiproton_number_unfolded_High_P1, TH1D *h_Antiproton_number_unfolded_High_P2, TH1D *h_Antiproton_number_unfolded_High_P3, TH1D *h_Antiproton_number_unfolded_High_P4, TH1D *h_Antiproton_number_unfolded_High_P5, TH1D *h_Antiproton_number_unfolded_High_PMinus1, 

TH1D *h_Proton_number_unfolded_High_P0, TH1D *h_Proton_number_unfolded_High_P0VGGNN, TH1D *h_Proton_number_unfolded_High_P1, TH1D *h_Proton_number_unfolded_High_P2, TH1D *h_Proton_number_unfolded_High_P3, TH1D *h_Proton_number_unfolded_High_P4, TH1D *h_Proton_number_unfolded_High_P5, TH1D *h_Proton_number_unfolded_High_PMinus1, 

TGraph *g_StatisticalError_mixedPattern, TGraphErrors *g_HighResult_P0, TGraphErrors *g_HighResult_P0VGGNN, TGraphErrors *g_HighResult_P1, TGraphErrors *g_HighResult_P2, TGraphErrors *g_HighResult_P3, TGraphErrors *g_HighResult_P4, TGraphErrors *g_HighResult_P5, TGraphErrors *g_HighResult_PMinus1, 

TH1D *h_Delta_antiproton_mixedPattern, TH1D *h_Delta_antiproton_P0, TH1D *h_Delta_antiproton_P0VGGNN, TH1D *h_Delta_antiproton_P1, TH1D *h_Delta_antiproton_P2, TH1D *h_Delta_antiproton_P3, TH1D *h_Delta_antiproton_P4, TH1D *h_Delta_antiproton_P5, TH1D *h_Delta_antiproton_PMinus1, 

TH1D *h_Proton_number_Raw_High_P0, TH1D *h_Proton_number_Raw_High_P0VGGNN, TH1D *h_Proton_number_Raw_High_P1, TH1D *h_Proton_number_Raw_High_P2, TH1D *h_Proton_number_Raw_High_P3, TH1D *h_Proton_number_Raw_High_P4, TH1D *h_Proton_number_Raw_High_P5, TH1D *h_Proton_number_Raw_High_PMinus1, TH1D *h_Proton_number_Raw_High_mixedPattern,

TH1D *h_Antiproton_number_Raw_High_P0, TH1D *h_Antiproton_number_Raw_High_P0VGGNN, TH1D *h_Antiproton_number_Raw_High_P1, TH1D *h_Antiproton_number_Raw_High_P2, TH1D *h_Antiproton_number_Raw_High_P3, TH1D *h_Antiproton_number_Raw_High_P4, TH1D *h_Antiproton_number_Raw_High_P5, TH1D *h_Antiproton_number_Raw_High_PMinus1, TH1D *h_Antiproton_number_Raw_High_mixedPattern,

std::string P0Useage){

    // P0 Useage: BDT or VGGNN.
    TH1D *h_Antiproton_number_unfolded_High_P0Useage;
    TH1D *h_Proton_number_unfolded_High_P0Useage;
    TH1D *h_Delta_antiproton_P0Useage;
    TH1D *h_Proton_number_Raw_High_P0Useage;
    TH1D *h_Antiproton_number_Raw_High_P0Useage;
    if (P0Useage == "VGGNN"){
        h_Antiproton_number_unfolded_High_P0Useage = h_Antiproton_number_unfolded_High_P0VGGNN;
        h_Proton_number_unfolded_High_P0Useage     = h_Proton_number_unfolded_High_P0VGGNN;
        h_Delta_antiproton_P0Useage                = h_Delta_antiproton_P0VGGNN;
        h_Proton_number_Raw_High_P0Useage          = h_Proton_number_Raw_High_P0VGGNN;
        h_Antiproton_number_Raw_High_P0Useage      = h_Antiproton_number_Raw_High_P0VGGNN;
    }
    else if (P0Useage == "BDT"){
        h_Antiproton_number_unfolded_High_P0Useage = h_Antiproton_number_unfolded_High_P0;
        h_Proton_number_unfolded_High_P0Useage     = h_Proton_number_unfolded_High_P0;
        h_Delta_antiproton_P0Useage                = h_Delta_antiproton_P0;
        h_Proton_number_Raw_High_P0Useage          = h_Proton_number_Raw_High_P0;
        h_Antiproton_number_Raw_High_P0Useage      = h_Antiproton_number_Raw_High_P0;
    }

    // For h_Antiproton_number_unfolded_High_P0:  index1: 14.1-15.3; index13: 36.1-38.9; index27: 125-147; index28:147-175; index32:330-525; total:32points.
    // Below 38.9 GV：       neither L1 nor L9 is required.      Pattern = 0, 1, 2, 4;    3(no L2?), 5(no L2?), -1(noL2?)
    for (int q = Range1_FirstIndex; q <= Range1_LastIndex; q++){
        h_Antiproton_number_unfolded_High_mixedPattern->SetBinContent(q, h_Antiproton_number_unfolded_High_P0Useage->GetBinContent(q) + h_Antiproton_number_unfolded_High_P1->GetBinContent(q) + h_Antiproton_number_unfolded_High_P2->GetBinContent(q) + h_Antiproton_number_unfolded_High_P4->GetBinContent(q) );
        h_Proton_number_unfolded_High_mixedPattern    ->SetBinContent(q, h_Proton_number_unfolded_High_P0Useage->GetBinContent(q) + h_Proton_number_unfolded_High_P1->GetBinContent(q) + h_Proton_number_unfolded_High_P2->GetBinContent(q) + h_Proton_number_unfolded_High_P4->GetBinContent(q) );

        h_Delta_antiproton_mixedPattern          ->SetBinContent(q, sqrt( pow(h_Delta_antiproton_P0Useage->GetBinContent(q), 2) + pow(h_Delta_antiproton_P1->GetBinContent(q), 2) + pow(h_Delta_antiproton_P2->GetBinContent(q), 2) + pow(h_Delta_antiproton_P4->GetBinContent(q), 2) ) );
        h_Proton_number_Raw_High_mixedPattern    ->SetBinContent(q, h_Proton_number_Raw_High_P0Useage->GetBinContent(q) + h_Proton_number_Raw_High_P1->GetBinContent(q) + h_Proton_number_Raw_High_P2->GetBinContent(q) + h_Proton_number_Raw_High_P4->GetBinContent(q));
        h_Antiproton_number_Raw_High_mixedPattern->SetBinContent(q, h_Antiproton_number_Raw_High_P0Useage->GetBinContent(q) + h_Antiproton_number_Raw_High_P1->GetBinContent(q) + h_Antiproton_number_Raw_High_P2->GetBinContent(q) + h_Antiproton_number_Raw_High_P4->GetBinContent(q)); 

        g_StatisticalError_mixedPattern->SetPoint(q-1, g_StatisticalError_mixedPattern->GetX()[q-1], h_Delta_antiproton_mixedPattern->GetBinContent(q)/h_Proton_number_Raw_High_mixedPattern->GetBinContent(q));
    }

    // From 38.9 to 147 GV： either L1 or L9 is required.        Pattern = 0, 1, 2;       3(noL2?), 5(noL2?)
    for (int q = Range2_FirstIndex; q <= Range2_LastIndex; q++){
        h_Antiproton_number_unfolded_High_mixedPattern->SetBinContent(q, h_Antiproton_number_unfolded_High_P0Useage->GetBinContent(q) + h_Antiproton_number_unfolded_High_P1->GetBinContent(q) + h_Antiproton_number_unfolded_High_P2->GetBinContent(q) );
        h_Proton_number_unfolded_High_mixedPattern->SetBinContent(q, h_Proton_number_unfolded_High_P0Useage->GetBinContent(q) + h_Proton_number_unfolded_High_P1->GetBinContent(q) + h_Proton_number_unfolded_High_P2->GetBinContent(q) );

        h_Delta_antiproton_mixedPattern->SetBinContent(q, sqrt( pow(h_Delta_antiproton_P0Useage->GetBinContent(q), 2) + pow(h_Delta_antiproton_P1->GetBinContent(q), 2) + pow(h_Delta_antiproton_P2->GetBinContent(q), 2) ) );
        h_Proton_number_Raw_High_mixedPattern->SetBinContent(q, h_Proton_number_Raw_High_P0Useage->GetBinContent(q) + h_Proton_number_Raw_High_P1->GetBinContent(q) + h_Proton_number_Raw_High_P2->GetBinContent(q));
        h_Antiproton_number_Raw_High_mixedPattern->SetBinContent(q, h_Antiproton_number_Raw_High_P0Useage->GetBinContent(q) + h_Antiproton_number_Raw_High_P1->GetBinContent(q) + h_Antiproton_number_Raw_High_P2->GetBinContent(q));

        g_StatisticalError_mixedPattern->SetPoint(q-1, g_StatisticalError_mixedPattern->GetX()[q-1], h_Delta_antiproton_mixedPattern->GetBinContent(q)/h_Proton_number_Raw_High_mixedPattern->GetBinContent(q));
    }

    // From 147 to 175 GV：  only L9 is required.                Pattern = 0, 2;          5(noL2?)
    for (int q = Range3_FirstIndex; q <= Range3_LastIndex; q++){
        h_Antiproton_number_unfolded_High_mixedPattern->SetBinContent(q, h_Antiproton_number_unfolded_High_P0Useage->GetBinContent(q) + h_Antiproton_number_unfolded_High_P2->GetBinContent(q));
        h_Proton_number_unfolded_High_mixedPattern->SetBinContent(q, h_Proton_number_unfolded_High_P0Useage->GetBinContent(q) + h_Proton_number_unfolded_High_P2->GetBinContent(q));

        h_Delta_antiproton_mixedPattern->SetBinContent(q, sqrt( pow(h_Delta_antiproton_P0Useage->GetBinContent(q), 2) + pow(h_Delta_antiproton_P2->GetBinContent(q), 2) ) );
        h_Proton_number_Raw_High_mixedPattern->SetBinContent(q, h_Proton_number_Raw_High_P0Useage->GetBinContent(q) + h_Proton_number_Raw_High_P2->GetBinContent(q));
        h_Antiproton_number_Raw_High_mixedPattern->SetBinContent(q, h_Antiproton_number_Raw_High_P0Useage->GetBinContent(q) + h_Antiproton_number_Raw_High_P2->GetBinContent(q));

        g_StatisticalError_mixedPattern->SetPoint(q-1, g_StatisticalError_mixedPattern->GetX()[q-1], h_Delta_antiproton_mixedPattern->GetBinContent(q)/h_Proton_number_Raw_High_mixedPattern->GetBinContent(q));
    }

    // Above 175 GV：        both L1 and L9 are required.        Pattern = 0
    for (int q = Range4_FirstIndex; q <= Range4_LastIndex; q++){
        h_Antiproton_number_unfolded_High_mixedPattern->SetBinContent(q, h_Antiproton_number_unfolded_High_P0Useage->GetBinContent(q) );
        h_Proton_number_unfolded_High_mixedPattern->SetBinContent(q, h_Proton_number_unfolded_High_P0Useage->GetBinContent(q) );

        h_Delta_antiproton_mixedPattern->SetBinContent(q, sqrt( pow(h_Delta_antiproton_P0Useage->GetBinContent(q), 2)) );
        h_Proton_number_Raw_High_mixedPattern->SetBinContent(q, h_Proton_number_Raw_High_P0Useage->GetBinContent(q));
        h_Antiproton_number_Raw_High_mixedPattern->SetBinContent(q, h_Antiproton_number_Raw_High_P0Useage->GetBinContent(q));

        g_StatisticalError_mixedPattern->SetPoint(q-1, g_StatisticalError_mixedPattern->GetX()[q-1], h_Delta_antiproton_mixedPattern->GetBinContent(q)/h_Proton_number_Raw_High_mixedPattern->GetBinContent(q));
    }
}


void Reset_StatisticalError(TGraphErrors *g_HighResult_P0, TGraphErrors *g_HighResult_P0VGGNN, TGraphErrors *g_HighResult_P1, TGraphErrors *g_HighResult_P2, TGraphErrors *g_HighResult_P3, TGraphErrors *g_HighResult_P4, TGraphErrors *g_HighResult_P5, TGraphErrors *g_HighResult_PMinus1, 
TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P0, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P1, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P2, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P3, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P4, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P5, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_PMinus1, 
TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern, TGraphErrors *g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern, TGraph *g_StatisticalError_mixedPattern){

    for (int i = 0; i < g_HighResult_P0->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_P0->SetPointError(i, 0, g_HighResult_P0->GetErrorY(i) );
    }

    for (int i = 0; i < g_HighResult_P0VGGNN->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN->SetPointError(i, 0, g_HighResult_P0VGGNN->GetErrorY(i) );
    }

    for (int i = 0; i < g_HighResult_P1->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_P1->SetPointError(i, 0, g_HighResult_P1->GetErrorY(i) );
    }

    for (int i = 0; i < g_HighResult_P2->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_P2->SetPointError(i, 0, g_HighResult_P2->GetErrorY(i) );
    }

    for (int i = 0; i < g_HighResult_P3->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_P3->SetPointError(i, 0, g_HighResult_P3->GetErrorY(i) );
    }

    for (int i = 0; i < g_HighResult_P4->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_P4->SetPointError(i, 0, g_HighResult_P4->GetErrorY(i) );
    }

    for (int i = 0; i < g_HighResult_P5->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_P5->SetPointError(i, 0, g_HighResult_P5->GetErrorY(i) );
    }

    for (int i = 0; i < g_HighResult_PMinus1->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_PMinus1->SetPointError(i, 0, g_HighResult_PMinus1->GetErrorY(i) );
    }

    for (int i = 0; i < g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->GetN(); i++){
        g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->SetPointError(i, 0, g_StatisticalError_mixedPattern->GetY()[i] );
        //cout<< "Error is " << g_StatisticalError_mixedPattern->GetY()[i] <<endl;
    }

    for (int i = 0; i < g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetN(); i++){
        g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->SetPointError(i, 0, g_StatisticalError_mixedPattern->GetY()[i] );
    }

}



//// Plot
void Plot_EffectiveAcceptanceRatio_AllPatterns(TGraph *gFitFunction_P0, TGraph *gFitFunction_P1, TGraph *gFitFunction_P2, TGraph *gFitFunction_P3, TGraph *gFitFunction_P4, TGraph *gFitFunction_P5, TGraph *gFitFunction_PMinus1, TGraphErrors *Effective_Acceptance_Ratio_P0, TGraphErrors *Effective_Acceptance_Ratio_P1, TGraphErrors *Effective_Acceptance_Ratio_P2, TGraphErrors *Effective_Acceptance_Ratio_P3, TGraphErrors *Effective_Acceptance_Ratio_P4, TGraphErrors *Effective_Acceptance_Ratio_P5, TGraphErrors *Effective_Acceptance_Ratio_PMinus1){

    TCanvas c_acc("c_acc","c_acc",1000,500);
    TPad *pad_acc = new TPad("pad_acc", "The pad1 acc", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    pad_acc->Draw();
    pad_acc->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    pad_acc->SetLeftMargin(0.18);

    Effective_Acceptance_Ratio_P0->Draw("AP");
    Effective_Acceptance_Ratio_P1->Draw("same P");
    Effective_Acceptance_Ratio_P2->Draw("same P");
    //Effective_Acceptance_Ratio_P3->Draw("same P");
    Effective_Acceptance_Ratio_P4->Draw("same P");
    //Effective_Acceptance_Ratio_P5->Draw("same P");
    //Effective_Acceptance_Ratio_PMinus1->Draw("same P");

    //gFitFunction_P0->Draw("same P");
    //gFitFunction_P1->Draw("same P");
    //gFitFunction_P2->Draw("same P");
    //gFitFunction_P3->Draw("same P");
    //gFitFunction_P4->Draw("same P");
    //gFitFunction_P5->Draw("same P");
    //gFitFunction_PMinus1->Draw("same P");

    Effective_Acceptance_Ratio_P0->SetMarkerStyle(15);
    Effective_Acceptance_Ratio_P0->SetMarkerColor(2);
    Effective_Acceptance_Ratio_P1->SetMarkerStyle(15);
    Effective_Acceptance_Ratio_P1->SetMarkerColor(3);
    Effective_Acceptance_Ratio_P2->SetMarkerStyle(15);
    Effective_Acceptance_Ratio_P2->SetMarkerColor(4);
    Effective_Acceptance_Ratio_P4->SetMarkerStyle(15);
    Effective_Acceptance_Ratio_P4->SetMarkerColor(6);
    Effective_Acceptance_Ratio_P0->GetFunction("function1")->SetLineColor(2);
    Effective_Acceptance_Ratio_P1->GetFunction("function1")->SetLineColor(3);
    Effective_Acceptance_Ratio_P2->GetFunction("function1")->SetLineColor(4);
    Effective_Acceptance_Ratio_P4->GetFunction("function1")->SetLineColor(6);

    TAxis * xaxis = Effective_Acceptance_Ratio_P0->GetXaxis();
    TAxis * yaxis = Effective_Acceptance_Ratio_P0->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Effective Acceptance Ratio");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(1.0, 1.15);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    //yaxis->SetTitleOffset(0.8);

    TLegend *legend_acc = new TLegend(0.6, 0.60, 0.8, 0.89);
    legend_acc->AddEntry(Effective_Acceptance_Ratio_P0     , "Pattern 0","lpf");
    legend_acc->AddEntry(Effective_Acceptance_Ratio_P1     , "Pattern 1","lpf");
    legend_acc->AddEntry(Effective_Acceptance_Ratio_P2     , "Pattern 2","lpf");
    legend_acc->AddEntry(Effective_Acceptance_Ratio_P4     , "Pattern 4","lpf");
    legend_acc->SetTextSize(0.05);
    legend_acc->SetTextFont(62);
    legend_acc->Draw();

    // Remove StaBox out of canvas.
    TPaveStats *ps0 = (TPaveStats *)Effective_Acceptance_Ratio_P0->GetListOfFunctions()->FindObject("stats");
    ps0->SetX1NDC(4.15);
    ps0->SetX2NDC(5.55);
    ps0->SetOptStat(0);

    TPaveStats *ps1 = (TPaveStats *)Effective_Acceptance_Ratio_P1->GetListOfFunctions()->FindObject("stats");
    ps1->SetX1NDC(4.15);
    ps1->SetX2NDC(5.55);
    ps1->SetOptStat(0);

    TPaveStats *ps2 = (TPaveStats *)Effective_Acceptance_Ratio_P2->GetListOfFunctions()->FindObject("stats");
    ps2->SetX1NDC(4.15);
    ps2->SetX2NDC(5.55);
    ps2->SetOptStat(0);

    TPaveStats *ps4 = (TPaveStats *)Effective_Acceptance_Ratio_P4->GetListOfFunctions()->FindObject("stats");
    ps4->SetX1NDC(4.15); 
    ps4->SetX2NDC(5.55);
    ps4->SetOptStat(0);

    gPad->Update();

    c_acc.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("EffectiveAcceptanceRatio") + string("_LinearX.pdf")).c_str());
    gPad->SetLogx();
    c_acc.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("EffectiveAcceptanceRatio") + string("_LogX.pdf")).c_str());
}

void Plot_EffectiveAcceptanceRatio_SinglePattern(TGraphErrors *Effective_Acceptance_Ratio_PX, std::string PatternName){

    TCanvas c_acc("c_acc","c_acc",1000,500);
    TPad *pad_acc = new TPad("pad_acc", "The pad1 acc", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    pad_acc->Draw();
    pad_acc->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    pad_acc->SetLeftMargin(0.18);

    Effective_Acceptance_Ratio_PX->Draw("AP");

    TAxis * xaxis = Effective_Acceptance_Ratio_PX->GetXaxis();
    TAxis * yaxis = Effective_Acceptance_Ratio_PX->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Effective Acceptance Ratio");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(1.0, 1.15);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    //yaxis->SetTitleOffset(0.8);
    
    TLegend *legend_acc = new TLegend(0.2, 0.75, 0.4, 0.89);
    legend_acc->AddEntry(Effective_Acceptance_Ratio_PX, (PatternName).c_str(), "lpf");
    legend_acc->SetTextSize(0.05);
    legend_acc->SetTextFont(62);
    legend_acc->Draw();

    TPaveStats *psX = (TPaveStats *)Effective_Acceptance_Ratio_PX->GetListOfFunctions()->FindObject("stats");
    psX->SetX1NDC(0.55);
    psX->SetX2NDC(0.95);
    psX->SetOptStat(0);

    gPad->Update();

    c_acc.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("EffectiveAcceptanceRatio_") + PatternName + string("_LinearX.pdf")).c_str());
    gPad->SetLogx();
    c_acc.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("EffectiveAcceptanceRatio_") + PatternName + string("_LogX.pdf")).c_str());
}


void Plot_RawRatioAndUnfoldedRatio_Compare(std::string issversion, TGraphErrors *g_Raw_Unfold_Compare_P0VGGNN, TGraphErrors *g_Raw_Unfold_Compare_P1, TGraphErrors *g_Raw_Unfold_Compare_P2, TGraphErrors *g_Raw_Unfold_Compare_P4){

    TCanvas c1("c1","c1",1000,500);
    TPad *p1 = new TPad("p1", "", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1->Draw();
    p1->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);


    g_Raw_Unfold_Compare_P0VGGNN->Draw("AP");
    g_Raw_Unfold_Compare_P1->Draw("same P");
    g_Raw_Unfold_Compare_P2->Draw("same P");
    g_Raw_Unfold_Compare_P4->Draw("same P");

    g_Raw_Unfold_Compare_P0VGGNN->SetTitle("");

    g_Raw_Unfold_Compare_P0VGGNN->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P0VGGNN->SetMarkerColor(7);
    g_Raw_Unfold_Compare_P1->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P1->SetMarkerColor(41);
    g_Raw_Unfold_Compare_P2->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P2->SetMarkerColor(8);
    g_Raw_Unfold_Compare_P4->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P4->SetMarkerColor(28);


    TAxis * xaxis = g_Raw_Unfold_Compare_P0VGGNN->GetXaxis();
    TAxis * yaxis = g_Raw_Unfold_Compare_P0VGGNN->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("#frac{Raw Ratio - Unfolded Ratio}{Unfolded Ratio}");
    xaxis->SetMoreLogLabels();
    xaxis->SetNoExponent();
    yaxis->SetRangeUser(-0.5, 0.5);

    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(1.1);

    gPad->SetBottomMargin(0.18);
    gPad->SetLeftMargin(0.15);

    gPad->Update();
    c1.Update();

    TLegend *legend = new TLegend(0.70, 0.70, 0.87, 0.89);
    legend->AddEntry(g_Raw_Unfold_Compare_P0VGGNN, "FullSpan"  ,"p"); // Pattern 0
    legend->AddEntry(g_Raw_Unfold_Compare_P1     , "Inner + L1","p"); // Pattern 1
    legend->AddEntry(g_Raw_Unfold_Compare_P2     , "Inner + L9","p"); // Pattern 2
    legend->AddEntry(g_Raw_Unfold_Compare_P4     , "Inner Only","p"); // Pattern 4
    legend->SetTextSize(0.04);
    legend->SetTextFont(62);
    legend->Draw();

    c1.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("RawUnfoldedRatioCompare") + string("_") + issversion + string(".pdf")).c_str());
    


}


void Plot_Ratio_HighRange_AllPattern(TGraphErrors *gPhyrReortRatio, TGraphErrors *g_ratio_P0, TGraphErrors *g_ratio_P0VGGNN, TGraphErrors *g_ratio_P1, TGraphErrors *g_ratio_P2, TGraphErrors *g_ratio_P3, TGraphErrors *g_ratio_P4, TGraphErrors *g_ratio_P5, TGraphErrors *g_ratio_PMinus1, std::string issversion, string AcceptanceCorrectionStatue, TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern){

    TCanvas c_ratio("c_ratio","c_ratio",1000,500);
    TPad *pad_ratio = new TPad("pad_ratio", "The pad1 ratio", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    pad_ratio->Draw();
    pad_ratio->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    g_ratio_P0VGGNN->Draw("AP");
    //g_ratio_P0->Draw("same P");
    g_ratio_P1->Draw("same P");
    g_ratio_P2->Draw("same P");
    //g_ratio_P3->Draw("same P");
    g_ratio_P4->Draw("same P");
    //g_ratio_P5->Draw("same P");
    //g_ratio_PMinus1->Draw("same P");
    //g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->Draw("same P");
    gPhyrReortRatio->Draw("same P");

    g_ratio_P0VGGNN->SetTitle("");

    pad_ratio->SetLeftMargin(0.18);

    g_ratio_P0->SetMarkerStyle(15);
    g_ratio_P0->SetMarkerColor(5);
    g_ratio_P0VGGNN->SetMarkerStyle(15);
    g_ratio_P0VGGNN->SetMarkerColor(2);
    g_ratio_P1->SetMarkerStyle(15);
    g_ratio_P1->SetMarkerColor(3);
    g_ratio_P2->SetMarkerStyle(15);
    g_ratio_P2->SetMarkerColor(4);
    g_ratio_P3->SetMarkerStyle(15);
    g_ratio_P3->SetMarkerColor(9);
    g_ratio_P4->SetMarkerStyle(15);
    g_ratio_P4->SetMarkerColor(6);
    g_ratio_P5->SetMarkerStyle(15);
    g_ratio_P5->SetMarkerColor(46);
    g_ratio_PMinus1->SetMarkerStyle(15);
    g_ratio_PMinus1->SetMarkerColor(93);
    g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->SetMarkerStyle(15);
    g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->SetMarkerColor(1);
    gPhyrReortRatio->SetMarkerStyle(15);
    gPhyrReortRatio->SetMarkerColor(7);

    TAxis * xaxis = g_ratio_P0VGGNN->GetXaxis();
    TAxis * yaxis = g_ratio_P0VGGNN->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Antiproton to Proton Ratio");
    xaxis->SetLimits(10, 600);
    //yaxis->SetRangeUser(0.000000001, 0.00070);
    //yaxis->SetRangeUser(0.000000001, 0.0013);
    yaxis->SetRangeUser(0.0001, 0.0029);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(1.6);

    xaxis->SetMoreLogLabels();

    TLegend *legend_highratio = new TLegend(0.23, 0.60, 0.50, 0.89);
    //legend_highratio->AddEntry(g_ratio_P0     , "Pattern 0","p");
    legend_highratio->AddEntry(g_ratio_P0VGGNN, "Pattern 0","p");    // Here only One P0 is showed, use VGGNN.
    legend_highratio->AddEntry(g_ratio_P1     , "Pattern 1","p");
    legend_highratio->AddEntry(g_ratio_P2     , "Pattern 2","p");
    //legend_highratio->AddEntry(g_ratio_P3     , "Pattern 3","p");
    legend_highratio->AddEntry(g_ratio_P4     , "Pattern 4","p");
    //legend_highratio->AddEntry(g_ratio_P5     , "Pattern 5","p");
    //legend_highratio->AddEntry(g_ratio_PMinus1, "Pattern -1","p");
    //legend_highratio->AddEntry(g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern, "MixedPattern","p");
    legend_highratio->AddEntry(gPhyrReortRatio, "PhysicsReport2021","p");
    legend_highratio->SetTextSize(0.05);
    legend_highratio->SetTextFont(62);
    legend_highratio->Draw();

    gPad->SetLogy();

    //c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_") + AcceptanceCorrectionStatue + string("_") + issversion + string("_Linear.pdf")).c_str());
    gPad->SetLogx();
    c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_") + AcceptanceCorrectionStatue + string("_") + issversion + string("_XLog.pdf")).c_str());
}


void Plot_Ratio_HighRange_PatternSelection(TGraphErrors *gPhyrReortRatio, TGraphErrors *g_Publisehd_Ratio, int Range1_FirstIndex, int Range1_LastIndex, int Range2_FirstIndex, int Range2_LastIndex, int Range3_FirstIndex, int Range3_LastIndex, int Range4_FirstIndex, int Range4_LastIndex, 
TGraphErrors *g_ratio_P0, TGraphErrors *g_ratio_P0VGGNN, TGraphErrors *g_ratio_P1, TGraphErrors *g_ratio_P2, TGraphErrors *g_ratio_P3, TGraphErrors *g_ratio_P4, TGraphErrors *g_ratio_P5, std::string issversion, string AcceptanceCorrectionStatue, TGraphErrors *g_ratio_mixedPattern){

    for (int i=0; i < (Range4_LastIndex-Range4_FirstIndex+1); i++){
        g_ratio_P2->RemovePoint(g_ratio_P2->GetN()-1);
        g_ratio_P5->RemovePoint(g_ratio_P5->GetN()-1);
    }

    for (int i=0; i < (Range4_LastIndex-Range4_FirstIndex+1) + (Range3_LastIndex-Range3_FirstIndex+1); i++){
        g_ratio_P1->RemovePoint(g_ratio_P1->GetN()-1);
        g_ratio_P3->RemovePoint(g_ratio_P3->GetN()-1);
    }

    for (int i=0; i < (Range4_LastIndex-Range4_FirstIndex+1) + (Range3_LastIndex-Range3_FirstIndex+1) + (Range2_LastIndex-Range2_FirstIndex)+1; i++){
        g_ratio_P4->RemovePoint(g_ratio_P4->GetN()-1);
    }

    TCanvas c_ratio("c_ratio","c_ratio",1000,500);

    g_ratio_P0->Draw("AP");
    g_ratio_P0VGGNN->Draw("same P");
    g_ratio_P1->Draw("same P");
    g_ratio_P2->Draw("same P");
    //g_ratio_P3->Draw("same P");
    g_ratio_P4->Draw("same P");
    //g_ratio_P5->Draw("same P");
    //g_Publisehd_Ratio->Draw("same P");
    gPhyrReortRatio->Draw("same P");
    g_ratio_mixedPattern->Draw("same P");

    g_ratio_P0->SetTitle("");
    g_ratio_P0->SetMarkerStyle(15);
    g_ratio_P0->SetMarkerColor(2);
    g_ratio_P0VGGNN->SetMarkerStyle(15);
    g_ratio_P0VGGNN->SetMarkerColor(5);
    g_ratio_P1->SetMarkerStyle(15);
    g_ratio_P1->SetMarkerColor(3);
    g_ratio_P2->SetMarkerStyle(15);
    g_ratio_P2->SetMarkerColor(4);
    g_ratio_P3->SetMarkerStyle(15);
    g_ratio_P3->SetMarkerColor(9);
    g_ratio_P4->SetMarkerStyle(15);
    g_ratio_P4->SetMarkerColor(6);
    g_ratio_P5->SetMarkerStyle(15);
    g_ratio_P5->SetMarkerColor(46);
    g_Publisehd_Ratio->SetMarkerStyle(15);
    g_Publisehd_Ratio->SetMarkerColor(30);
    gPhyrReortRatio->SetMarkerStyle(15);
    gPhyrReortRatio->SetMarkerColor(50);
    g_ratio_mixedPattern->SetMarkerStyle(15);
    g_ratio_mixedPattern->SetMarkerColor(1);


    TAxis * xaxis = g_ratio_P0->GetXaxis();
    TAxis * yaxis = g_ratio_P0->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Antiproton to Proton Ratio");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0.00005, 0.00035);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(0.8);
    xaxis->SetMoreLogLabels();
    //pad_ratio->SetLeftMargin(0.5);
    yaxis->SetTitleOffset(0.9);

    TLegend *legend_highratio = new TLegend(0.15, 0.60, 0.44, 0.89);
    legend_highratio->AddEntry(g_ratio_P0, "Pattern 0","lpf");
    legend_highratio->AddEntry(g_ratio_P0VGGNN, "Pattern 0VGGNN","lpf");
    legend_highratio->AddEntry(g_ratio_P1, "Pattern 1","lpf");
    legend_highratio->AddEntry(g_ratio_P2, "Pattern 2","lpf");
    //legend_highratio->AddEntry(g_ratio_P3, "Pattern 3","lpf");
    legend_highratio->AddEntry(g_ratio_P4, "Pattern 4","lpf");
    //legend_highratio->AddEntry(g_ratio_P5, "Pattern 5","lpf");
    //legend_highratio->AddEntry(g_Publisehd_Ratio, "2016PRL","lpf");
    legend_highratio->AddEntry(gPhyrReortRatio, "PhysicsReport2021","lpf");
    legend_highratio->AddEntry(g_ratio_mixedPattern, "MixedPattern","lpf");
    legend_highratio->SetTextSize(0.05);
    legend_highratio->SetTextFont(62);
    legend_highratio->Draw();

    c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_PattternSelectionRange_") + AcceptanceCorrectionStatue + string("_") + issversion + string("_Linear.pdf")).c_str());
    gPad->SetLogx();
    c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_PattternSelectionRange_") + AcceptanceCorrectionStatue + string("_") + issversion + string("_XLog.pdf")).c_str());

}


void Plot_Ratio_HighRange_SinglePatternCompareReference(TGraphErrors *gPhyrReortRatio, TGraphErrors *g_ratio_PX, std::string patternname, std::string issversion, string AcceptanceCorrectionStatue) {

    TCanvas c_ratio("c_ratio","c_ratio",1000,500);

    gPhyrReortRatio->Draw("AP");
    g_ratio_PX->Draw("P same");

    gPhyrReortRatio->SetTitle("");
    gPhyrReortRatio->SetMarkerStyle(15);
    gPhyrReortRatio->SetMarkerColor(1);
    g_ratio_PX->SetMarkerStyle(15);
    g_ratio_PX->SetMarkerColor(2);

    TAxis * xaxis = gPhyrReortRatio->GetXaxis();
    TAxis * yaxis = gPhyrReortRatio->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Antiproton to Proton Ratio");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0.00005, 0.00035);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(0.8);
    xaxis->SetMoreLogLabels();
    yaxis->SetTitleOffset(0.9);

    TLegend *legend_highratio = new TLegend(0.15, 0.60, 0.44, 0.89);
    legend_highratio->AddEntry(gPhyrReortRatio, "PhysicsReport2021"  , "lpf");
    legend_highratio->AddEntry(g_ratio_PX,      (patternname).c_str(), "lpf");
    legend_highratio->SetTextSize(0.05);
    legend_highratio->SetTextFont(62);
    legend_highratio->Draw();

    c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_SinglePatttern_") + patternname + string("_") + AcceptanceCorrectionStatue + string("_") + issversion + string("_Linear.pdf")).c_str());
    gPad->SetLogx();
    c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_SinglePatttern_") + patternname + string("_") + AcceptanceCorrectionStatue + string("_") + issversion + string("_XLog.pdf")).c_str());
}



void Plot_Ratio_HighRange_mixedPattern(TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern, TGraphErrors *gPhyrReortRatio, TGraphErrors *g_Publisehd_Ratio, TGraphErrors *g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern){

    TCanvas c_ratio("c_ratio","c_ratio",1000,500);

    //g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->Draw("AP");
    g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->Draw("AP");
    g_Publisehd_Ratio->Draw("same P");

    g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->SetMarkerStyle(15);
    g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->SetMarkerColor(2);
    g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->SetMarkerSize(1.0);
    g_Publisehd_Ratio->SetMarkerStyle(15);
    g_Publisehd_Ratio->SetMarkerColor(30);
    g_Publisehd_Ratio->SetMarkerSize(1.0);
    g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->SetMarkerStyle(15);
    g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->SetMarkerColor(95);
    g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->SetMarkerSize(1.0);


    TAxis * xaxis = g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetXaxis();
    TAxis * yaxis = g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Antiproton to Proton Ratio");
    xaxis->SetLimits(10, 600);
    //yaxis->SetRangeUser(0.000000001, 0.00045);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(0.8);
    xaxis->SetMoreLogLabels();

    TLegend *legend_highratio = new TLegend(0.15, 0.65, 0.45, 0.87);
    //legend_highratio->AddEntry(g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern, "MixedPattern(NoAcceptanceCorrection)","lpf");
    legend_highratio->AddEntry(g_Publisehd_Ratio                                         , "PublisehdRatio","lpf");
    legend_highratio->AddEntry(g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern   , "MixedPattern(WithAccCorrec)");
    legend_highratio->SetTextSize(0.03);
    legend_highratio->SetTextFont(62);
    legend_highratio->Draw();

    c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_HighRange_mixedPattern") + string("_Linear.pdf")).c_str());
    gPad->SetLogx();
    c_ratio.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Ratio_HighRange_mixedPattern") + string("_XLog.pdf")).c_str());

}


void Plot_StatisticalError_ComparePatterns(TGraph *g_StatisticalError_P0, TGraph *g_StatisticalError_P0VGGNN, TGraph *g_StatisticalError_P1, TGraph *g_StatisticalError_P2, TGraph *g_StatisticalError_P3, TGraph *g_StatisticalError_P4, TGraph *g_StatisticalError_P5, TGraph *g_StatisticalError_mixedPattern, TH1D *h_Antiproton_number_unfolded_High_P0, TH1D *h_Antiproton_number_unfolded_High_P0VGGNN, TH1D *h_Antiproton_number_unfolded_High_P1, TH1D *h_Antiproton_number_unfolded_High_P2, TH1D *h_Antiproton_number_unfolded_High_P3, TH1D *h_Antiproton_number_unfolded_High_P4, TH1D *h_Antiproton_number_unfolded_High_P5){

    TCanvas c_StatisticalError("c_StatisticalError","c_StatisticalError",1000,500);

    g_StatisticalError_P0VGGNN->Draw("AP");
    //g_StatisticalError_P0->Draw("same P");
    g_StatisticalError_P1->Draw("same P");
    g_StatisticalError_P2->Draw("same P");
    //g_StatisticalError_P3->Draw("same P");
    g_StatisticalError_P4->Draw("same P");
    //g_StatisticalError_P5->Draw("same P");
    g_StatisticalError_mixedPattern->Draw("same P");

    g_StatisticalError_P0VGGNN->SetTitle("");

    g_StatisticalError_P0->SetMarkerStyle(15);
    g_StatisticalError_P0->SetMarkerColor(5);
    g_StatisticalError_P0VGGNN->SetMarkerStyle(15);
    g_StatisticalError_P0VGGNN->SetMarkerColor(2);
    g_StatisticalError_P1->SetMarkerStyle(15);
    g_StatisticalError_P1->SetMarkerColor(3);
    g_StatisticalError_P2->SetMarkerStyle(15);
    g_StatisticalError_P2->SetMarkerColor(4);
    g_StatisticalError_P3->SetMarkerStyle(15);
    g_StatisticalError_P3->SetMarkerColor(7);
    g_StatisticalError_P4->SetMarkerStyle(15);
    g_StatisticalError_P4->SetMarkerColor(6);
    g_StatisticalError_P5->SetMarkerStyle(15);
    g_StatisticalError_P5->SetMarkerColor(46);
    g_StatisticalError_mixedPattern->SetMarkerStyle(15);
    g_StatisticalError_mixedPattern->SetMarkerColor(1);

    TAxis * xaxis = g_StatisticalError_P0VGGNN->GetXaxis();
    TAxis * yaxis = g_StatisticalError_P0VGGNN->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Antiproton to Proton ratio Statistical Error");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0.000000001, 0.00020);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.047);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    //xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(1.1);
    xaxis->SetMoreLogLabels();

    TLegend *legend_StaError = new TLegend(0.55, 0.60, 0.85, 0.87);
    legend_StaError->AddEntry(g_StatisticalError_mixedPattern, "AllPattern","lpf");
    //legend_StaError->AddEntry(g_StatisticalError_P0          , "Pattern 0","lpf");
    legend_StaError->AddEntry(g_StatisticalError_P0VGGNN     , "Pattern 0","lpf");  // Here only use one P0 as representative one, namely VGGNN.
    legend_StaError->AddEntry(g_StatisticalError_P1          , "Pattern 1","lpf");
    legend_StaError->AddEntry(g_StatisticalError_P2          , "Pattern 2","lpf");
    //legend_StaError->AddEntry(g_StatisticalError_P3          , "Pattern 3","lpf");
    legend_StaError->AddEntry(g_StatisticalError_P4          , "Pattern 4","lpf");
    //legend_StaError->AddEntry(g_StatisticalError_P5          , "Pattern 5","lpf");
    legend_StaError->SetTextSize(0.05);
    legend_StaError->SetTextFont(62);
    legend_StaError->Draw();

   c_StatisticalError.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("StatisticalError_HighRange_ComparePatterns") + string(".pdf")).c_str());
} 

void Plot_StatisticalErrorCompareWithSquareRootN(TH1D *h_Antiproton_number_unfolded_High_P0, TH1D *h_Antiproton_number_unfolded_High_P0VGGNN, TH1D *h_Antiproton_number_unfolded_High_P1, TH1D *h_Antiproton_number_unfolded_High_P2, TH1D *h_Antiproton_number_unfolded_High_P3, TH1D *h_Antiproton_number_unfolded_High_P4, TH1D *h_Antiproton_number_unfolded_High_P5, TH1D *h_Delta_antiproton_P0, TH1D *h_Delta_antiproton_P0VGGNN, TH1D *h_Delta_antiproton_P1, TH1D *h_Delta_antiproton_P2, TH1D *h_Delta_antiproton_P3, TH1D *h_Delta_antiproton_P4, TH1D *h_Delta_antiproton_P5){

    TH1D *h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN = h_Antiproton_number_unfolded_High_P0VGGNN;
    TH1D *h_SquareRootCompareAntiprotonNumber_unfolded_High_P1 = h_Antiproton_number_unfolded_High_P1;
    TH1D *h_SquareRootCompareAntiprotonNumber_unfolded_High_P2 = h_Antiproton_number_unfolded_High_P2;
    TH1D *h_SquareRootCompareAntiprotonNumber_unfolded_High_P4 = h_Antiproton_number_unfolded_High_P4;

    for (int i=1; i<h_Antiproton_number_unfolded_High_P0VGGNN->GetNbinsX(); i++){
        h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN->SetBinContent(i, h_Delta_antiproton_P0VGGNN->GetBinContent(i) / sqrt(h_Antiproton_number_unfolded_High_P0VGGNN->GetBinContent(i)));
    }
    for (int i=1; i<h_Antiproton_number_unfolded_High_P1->GetNbinsX(); i++){
        h_SquareRootCompareAntiprotonNumber_unfolded_High_P1->SetBinContent(i, h_Delta_antiproton_P1->GetBinContent(i) / sqrt(h_Antiproton_number_unfolded_High_P1->GetBinContent(i)));
    }
    for (int i=1; i<h_Antiproton_number_unfolded_High_P2->GetNbinsX(); i++){
        h_SquareRootCompareAntiprotonNumber_unfolded_High_P2->SetBinContent(i, h_Delta_antiproton_P2->GetBinContent(i) / sqrt(h_Antiproton_number_unfolded_High_P2->GetBinContent(i)));
    }
    for (int i=1; i<h_Antiproton_number_unfolded_High_P4->GetNbinsX(); i++){
        h_SquareRootCompareAntiprotonNumber_unfolded_High_P4->SetBinContent(i, h_Delta_antiproton_P4->GetBinContent(i) / sqrt(h_Antiproton_number_unfolded_High_P4->GetBinContent(i)));
    }

    h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN->SetStats(0);
    h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN->SetTitle("");

    TCanvas c_StaCom("c_StaCom","c_StaCom",1000,500);

    h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN->Draw("HIST");
    h_SquareRootCompareAntiprotonNumber_unfolded_High_P1->Draw("same HIST");
    h_SquareRootCompareAntiprotonNumber_unfolded_High_P2->Draw("same HIST");
    h_SquareRootCompareAntiprotonNumber_unfolded_High_P4->Draw("same HIST");

    h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN->SetLineColor(2);
    h_SquareRootCompareAntiprotonNumber_unfolded_High_P1->SetLineColor(4);
    h_SquareRootCompareAntiprotonNumber_unfolded_High_P2->SetLineColor(6);
    h_SquareRootCompareAntiprotonNumber_unfolded_High_P4->SetLineColor(9);

    TAxis * xaxis = h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN->GetXaxis();
    TAxis * yaxis = h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Fit / SquareRoot");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0, 6);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(0.8);
    xaxis->SetMoreLogLabels();

    TLegend *legend_StaCompaSquareRoot = new TLegend(0.15, 0.60, 0.4, 0.87);
    legend_StaCompaSquareRoot->AddEntry(h_SquareRootCompareAntiprotonNumber_unfolded_High_P0VGGNN, "Pattern 0","lpf");
    legend_StaCompaSquareRoot->AddEntry(h_SquareRootCompareAntiprotonNumber_unfolded_High_P1, "Pattern 1","lpf");
    legend_StaCompaSquareRoot->AddEntry(h_SquareRootCompareAntiprotonNumber_unfolded_High_P2, "Pattern 2","lpf");
    legend_StaCompaSquareRoot->AddEntry(h_SquareRootCompareAntiprotonNumber_unfolded_High_P4, "Pattern 4","lpf");
    legend_StaCompaSquareRoot->SetTextSize(0.05);
    legend_StaCompaSquareRoot->SetTextFont(62);
    legend_StaCompaSquareRoot->Draw();

    TLine line1(38.9, 0, 38.9, 6);
    line1.Draw("same");
    TLine line2(147, 0, 147, 6);
    line2.Draw("same");
    TLine line3(175, 0, 175, 6);
    line3.Draw("same");
    line1.SetLineStyle(2);
    line2.SetLineStyle(2);
    line3.SetLineStyle(2);

    gPad->SetLogx();

    c_StaCom.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("StatisticalErrorComareWithSuareRoot_HighRange_ComparePatterns") + string(".pdf")).c_str());



}


void Plot_Chi2(int Range1_FirstIndex, int Range1_LastIndex, int Range2_FirstIndex, int Range2_LastIndex, int Range3_FirstIndex, int Range3_LastIndex, int Range4_FirstIndex, int Range4_LastIndex, TGraph *g_FitChi2_P0, TGraph *g_FitChi2_P1, TGraph *g_FitChi2_P2, TGraph *g_FitChi2_P3, TGraph *g_FitChi2_P4, TGraph *g_FitChi2_P5){

    /*
    // Only show used range.
    for (int i=0; i < (Range4_LastIndex-Range4_FirstIndex+1); i++){
        g_FitChi2_P2->RemovePoint(g_FitChi2_P2->GetN()-1);
        g_FitChi2_P5->RemovePoint(g_FitChi2_P5->GetN()-1);
    }

    for (int i=0; i < (Range4_LastIndex-Range4_FirstIndex+1) + (Range3_LastIndex-Range3_FirstIndex+1); i++){
        g_FitChi2_P1->RemovePoint(g_FitChi2_P1->GetN()-1);
        g_FitChi2_P3->RemovePoint(g_FitChi2_P3->GetN()-1);
    }

    for (int i=0; i < (Range4_LastIndex-Range4_FirstIndex+1) + (Range3_LastIndex-Range3_FirstIndex+1) + (Range2_LastIndex-Range2_FirstIndex)+1; i++){
        g_FitChi2_P4->RemovePoint(g_FitChi2_P4->GetN()-1);
    }
    */

    // Plot
    TCanvas c_chi2("c_chii2","c_chi2",1000,500);

    g_FitChi2_P0->Draw("AP");
    g_FitChi2_P1->Draw("same P");
    g_FitChi2_P2->Draw("same P");
    //g_FitChi2_P3->Draw("same P");
    g_FitChi2_P4->Draw("same P");
    //g_FitChi2_P5->Draw("same P");

    g_FitChi2_P0->SetTitle("");

    g_FitChi2_P0->SetMarkerStyle(15);
    g_FitChi2_P0->SetMarkerColor(2);
    g_FitChi2_P1->SetMarkerStyle(15);
    g_FitChi2_P1->SetMarkerColor(3);
    g_FitChi2_P2->SetMarkerStyle(15);
    g_FitChi2_P2->SetMarkerColor(4);
    g_FitChi2_P3->SetMarkerStyle(15);
    g_FitChi2_P3->SetMarkerColor(7);
    g_FitChi2_P4->SetMarkerStyle(15);
    g_FitChi2_P4->SetMarkerColor(6);
    g_FitChi2_P5->SetMarkerStyle(15);
    g_FitChi2_P5->SetMarkerColor(46);


    TAxis * xaxis = g_FitChi2_P0->GetXaxis();
    TAxis * yaxis = g_FitChi2_P0->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("chi2/dof");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0, 20);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.05);
    yaxis->SetTitleSize(0.05);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(0.8);
    xaxis->SetMoreLogLabels();


    TLegend *legend_CHI2 = new TLegend(0.55, 0.60, 0.85, 0.87);
    legend_CHI2->AddEntry(g_FitChi2_P0          , "Pattern 0","lpf");
    legend_CHI2->AddEntry(g_FitChi2_P1          , "Pattern 1","lpf");
    legend_CHI2->AddEntry(g_FitChi2_P2          , "Pattern 2","lpf");
    //legend_CHI2->AddEntry(g_FitChi2_P3          , "Pattern 3","lpf");
    legend_CHI2->AddEntry(g_FitChi2_P4          , "Pattern 4","lpf");
    //legend_CHI2->AddEntry(g_FitChi2_P5          , "Pattern 5","lpf");
    legend_CHI2->SetTextSize(0.05);
    legend_CHI2->SetTextFont(62);
    legend_CHI2->Draw();

    gPad->SetLogx();

   c_chi2.SaveAs( ( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/plot_Allpatterns/") + string("Chi2dof_HighRange_ComparePatterns") + string(".pdf")).c_str());

}






