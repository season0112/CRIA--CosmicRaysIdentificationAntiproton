
std::tuple<TGraphErrors, TH1D, TH1D, TH1D, TH1D, TH1D, TH1D, TH1D> Load2016PRLResult(std::vector<double> subrange_450){
    // PRL2016 Result
    std::vector<double> RigidityBinPoint_Published( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.end());
    TGraphErrors gPublishedRatioError = TGraphErrors(57, RigidityBinPoint_Published.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.data(), 0, AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.data());

    // PRL2016 Result Error
    std::vector<double> Publishedratio_Error_SHORT (AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin()+26, AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin()+57); // 33 points, 14.1-last point.
    std::vector<double> Publishedratio_Statistic_error_SHORT (AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticErrorPRL.begin()+26, AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticErrorPRL.begin()+57);
    std::vector<double> Publishedratio_Systematic_error_SHORT (AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicErrorPRL.begin()+26, AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicErrorPRL.begin()+57);
    std::vector<double> Publishedratio_SHORT (AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin()+26, AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin()+57);

    TH1D hPublished_Statistic_error            = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Statistic_error_relative   = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Systematic_error           = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Systematic_error_relative  = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_total_error                = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_total_error_relative       = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Statistic_error_Proportion = TH1D("", "", 31, subrange_450.data());

    for (int n=0; n<31; n++){
        hPublished_total_error.SetBinContent(n+1, Publishedratio_Error_SHORT.at(n));
        hPublished_Statistic_error.SetBinContent(n+1, Publishedratio_Statistic_error_SHORT.at(n));
        hPublished_Statistic_error_relative.SetBinContent(n+1, Publishedratio_Statistic_error_SHORT.at(n)/Publishedratio_SHORT.at(n)*100);
        hPublished_Systematic_error.SetBinContent(n+1, Publishedratio_Systematic_error_SHORT.at(n));
        hPublished_Systematic_error_relative.SetBinContent(n+1, Publishedratio_Systematic_error_SHORT.at(n)/Publishedratio_SHORT.at(n)*100);
        hPublished_Statistic_error_Proportion.SetBinContent(n+1, Publishedratio_Statistic_error_SHORT.at(n)/Publishedratio_Error_SHORT.at(n));
        hPublished_total_error_relative.SetBinContent(n+1, Publishedratio_Error_SHORT.at(n)/Publishedratio_SHORT.at(n)*100 );
    }

    return {gPublishedRatioError, hPublished_Statistic_error, hPublished_Statistic_error_relative, hPublished_Systematic_error, hPublished_Systematic_error_relative, hPublished_total_error, hPublished_total_error_relative, hPublished_Statistic_error_Proportion};
}


std::tuple<TGraphAsymmErrors *, TGraphAsymmErrors *, TH1D, TH1D> LoadEffectiveAcceptance(std::string binningversion, std::string pattern, std::vector<double> subrange_450, std::vector<double> subrange_525){
    TFile *f5 = new TFile();
    TFile *f6 = new TFile();
    if ( binningversion == "450version" ){
        f5 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_Pattern") + pattern + std::string("_450version.root")).c_str() );
        f6 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_Pattern") + pattern + std::string("_450version.root")).c_str() ) ;
    }
    else if( binningversion == "525version" ){
        f5 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_Pattern") + pattern + std::string("_525version.root")).c_str() );
        f6 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_Pattern") + pattern + std::string("_525version.root")).c_str() );
    }
    TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton = (TGraphAsymmErrors*)f5->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Proton     = (TGraphAsymmErrors*)f6->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");

    TH1D h_EffectiveAcceptance_Antiproton;
    TH1D h_EffectiveAcceptance_Proton;
    if ( binningversion == "450version" ){
        h_EffectiveAcceptance_Antiproton = TH1D("", "", 31, subrange_450.data());  //31 points. 14.1-525
        h_EffectiveAcceptance_Proton     = TH1D("", "", 31, subrange_450.data());  //31 points. 14.1-525
    }
    else if( binningversion == "525version" ){
        h_EffectiveAcceptance_Antiproton = TH1D("", "", 32, subrange_525.data());  //32 points. 14.1-525
        h_EffectiveAcceptance_Proton     = TH1D("", "", 32, subrange_525.data());  //32 points. 14.1-525
    }
    Utilities::ConvertToHistogram ( g_EffectiveAcceptance_Antiproton, h_EffectiveAcceptance_Antiproton);
    Utilities::ConvertToHistogram ( g_EffectiveAcceptance_Proton    , h_EffectiveAcceptance_Proton);

    return {g_EffectiveAcceptance_Antiproton, g_EffectiveAcceptance_Proton, h_EffectiveAcceptance_Antiproton, h_EffectiveAcceptance_Proton};
}


std::tuple<TGraphErrors *> CalculateAcceptanceRatio(int bins, TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton, TGraphAsymmErrors *g_EffectiveAcceptance_Proton){
 
   Double_t xa[bins], Aa[bins], ea[bins], xp[bins], Ap[bins], ep[bins];

    TGraphErrors *g_EffectiveAcceptanceRatio_withManualError = new TGraphErrors(bins);   //bins=31,32 for 450,525 versions.
    for (int p = 27; p < 27+bins; ++p) {   // 27:14.7(14.1-15.3),  58:427.5(330-525)
        g_EffectiveAcceptance_Antiproton->GetPoint(p,xa[p-27],Aa[p-27]);
        g_EffectiveAcceptance_Proton->GetPoint(p,xp[p-27],Ap[p-27]);
        ea[p-27] = g_EffectiveAcceptance_Antiproton->GetErrorY(p);
        ep[p-27] = g_EffectiveAcceptance_Proton->GetErrorY(p);

        std::cout<< "Aa:" << Aa[p-27] << "  ea:" << ea[p-27] << std::endl;
        std::cout<< "Ap:" << Ap[p-27] << "  ep:" << ep[p-27] << std::endl;
        std::cout<< "\n" << std::endl;

        g_EffectiveAcceptanceRatio_withManualError->SetPoint(p-27,xa[p-27],Ap[p-27]/Aa[p-27]);
        g_EffectiveAcceptanceRatio_withManualError->SetPointError(p-27, 0.0, sqrt(pow(ep[p-27],2)/pow(Aa[p-27],2)+pow(Ap[p-27],2)/pow(Aa[p-27],4)*pow(ea[p-27],2)));
    }

    return {g_EffectiveAcceptanceRatio_withManualError};
}



void Plot_EffectiveAcceptance_Parametrilised(TGraphErrors *eff_A, std::string NNsuffix){
    TCanvas cc1("cc1","cc1",1000,500);
    eff_A->Draw("AP *");
    eff_A->SetTitle("");
    gStyle->SetOptFit(1);
    gPad->SetLogx();
    TAxis * cc1xaxis = eff_A->GetXaxis();
    TAxis * cc1yaxis = eff_A->GetYaxis();
    cc1xaxis->SetTitle("Rigidity (GV)");
    cc1yaxis->SetTitle("A_{p}/A_{#bar{p}}");
    cc1xaxis->SetMoreLogLabels();
    cc1.Update();
    cc1.SaveAs("plot_ratio/effective_acceptance_parametrilised.pdf");
}

void Plot_UnfoledRatio(TGraphErrors gRatio_unfolded, TGraphErrors g_ratio_raw_with_effective_correction, TGraph gRatio_published, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c_ratio("c_ratio","c_ratio",1000,500);

    gPad->SetLogx();

    gRatio_published.Draw("AP *");
    gRatio_unfolded.Draw("same P");
    g_ratio_raw_with_effective_correction.Draw("same P");

    gRatio_published.SetMarkerStyle(15);
    gRatio_published.SetMarkerColor(2);
    gRatio_unfolded.SetMarkerStyle(15);
    gRatio_unfolded.SetMarkerColor(4);
    g_ratio_raw_with_effective_correction.SetMarkerStyle(15);
    g_ratio_raw_with_effective_correction.SetMarkerColor(6);

    TAxis * xaxis = gRatio_published.GetXaxis();
    TAxis * yaxis = gRatio_published.GetYaxis();
    xaxis->SetMoreLogLabels();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("Antiproton to Proton ratio");
    xaxis->SetLimits(1.0,600);
    yaxis->SetRangeUser(0.000000001,0.00035);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.035);
    yaxis->SetTitleSize(0.04);
    xaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelFont(62);
    yaxis->SetLabelSize(0.05);

    TLegend *legend1 = new TLegend(0.45,0.15,0.78,0.35);
    legend1->AddEntry(&gRatio_published, "PRL paper 2016","lpf");
    legend1->AddEntry(&gRatio_unfolded, "Unfolded","lpf");
    legend1->AddEntry(&g_ratio_raw_with_effective_correction, "Raw","lpf");
    legend1->SetTextSize(0.05);
    legend1->SetTextFont(62);
    legend1->Draw();

    c_ratio.SaveAs( (std::string("plot_ratio/Ratio_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str() );
}


void Plot_Raw_Unfolded_Comparison(TGraphErrors graw_unfold_com, std::string issversion, std::string pattern, std::string NNsuffix){

    TCanvas c32("c32","c32",1000,500);

    gPad->SetLogx();

    graw_unfold_com.Draw("AP *");

    graw_unfold_com.SetMarkerStyle(15);
    graw_unfold_com.SetMarkerColor(2);
    graw_unfold_com.SetTitle("");

    TAxis * xaxis = graw_unfold_com.GetXaxis();
    TAxis * yaxis = graw_unfold_com.GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("#frac{Raw Ratio - Unfolded Ratio}{Unfolded Ratio}");
    //yaxis->SetTitle("#splitline{Raw Ratio-Unfolded Ratio}{Unfolded Ratio}");
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
        
    //yaxis->SetRangeUser(-1, 1.5);
    gPad->Update();
    c32.Update();

    //c32.SaveAs( (std::string("plot_ratio/RawUnfoldedComparisonPattern") + pattern + NNsuffix + std::string(".pdf")).c_str() );
    c32.SaveAs( (std::string("plot_ratio/raw_unfolded_comparison_Pattern_") + pattern + std::string("_") + issversion + std::string("_") + NNsuffix + std::string(".pdf")).c_str() );
}

void Plot_FitChi2(TGraph *g_fitchi2, std::string issversion, std::string pattern, std::string NNsuffix){

    TCanvas c_chi2("c_chi2","c_chi2",1000,500);

    g_fitchi2->Draw("AP");

    g_fitchi2->SetTitle("");
    g_fitchi2->GetXaxis()->SetTitle("Rigidity (GV)");
    g_fitchi2->GetYaxis()->SetTitle("Chi2/dof");

    g_fitchi2->GetXaxis()->SetTitleFont(62);
    g_fitchi2->GetXaxis()->SetTitleSize(0.045);
    g_fitchi2->GetYaxis()->SetTitleFont(62);
    g_fitchi2->GetYaxis()->SetTitleSize(0.045);

    g_fitchi2->GetXaxis()->SetLabelFont(62);
    g_fitchi2->GetXaxis()->SetLabelSize(0.05);
    g_fitchi2->GetYaxis()->SetLabelFont(62);
    g_fitchi2->GetYaxis()->SetLabelSize(0.05);
    g_fitchi2->SetMarkerStyle(15);
    g_fitchi2->SetMarkerColor(1);

    g_fitchi2->GetYaxis()->SetRangeUser(0,2);
    c_chi2.Update();

    c_chi2.SaveAs( (std::string("plot_ratio/Chi2dof_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str() );

}


void Plot_MM(TH2D *hMigrationMatrix, std::string pattern, std::string NNsuffix){

    /*
    // Normalize
    for (int q = 1; q <= hMigrationMatrix->GetNcells(); q++){    
        hMigrationMatrix->SetBinContent(q, hMigrationMatrix->GetBinContent(q) / hMigrationMatrix->GetMaximum());
    }
    */

    // Normalize2 (the sum of probabilities in the projection along Y axis is 1)
    for (int q=1; q<=hMigrationMatrix->GetNbinsX(); q++ ){

        double sum0 = 0;
        for (int p=1; p<=hMigrationMatrix->GetNbinsY(); p++){
            sum0 = sum0 + hMigrationMatrix->GetBinContent(q,p);
        }

        for (int p=1; p<=hMigrationMatrix->GetNbinsY(); p++){
            hMigrationMatrix->SetBinContent(q, p, hMigrationMatrix->GetBinContent(q,p)/sum0);
        }
    }
    /*
    //check:
    std::cout<< "check2:" <<std::endl;
    for (int q=1; q<=hMigrationMatrix->GetNbinsX(); q++ ){
        double sum2=0;
        for (int p=1; p<=hMigrationMatrix->GetNbinsY(); p++){
            sum2 = sum2+hMigrationMatrix->GetBinContent(q, p);
        }
        std::cout<< "sum2:" << sum2 <<std::endl;
    }
    */

    // Plot
    TCanvas c4("c4", "c4", 1000, 800);

    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();
    gStyle->SetOptStat(0);

    hMigrationMatrix->Draw("COL");

    double min = hMigrationMatrix->GetYaxis()->GetBinLowEdge(1);
    double max = hMigrationMatrix->GetYaxis()->GetBinUpEdge(hMigrationMatrix->GetXaxis()->GetNbins());
    TPaletteAxis *p = new TPaletteAxis(max+50, min, max+250, max, hMigrationMatrix);
    p->Draw();

    TAxis* xaxis = hMigrationMatrix->GetXaxis();
    TAxis* yaxis = hMigrationMatrix->GetYaxis();
    TAxis* zaxis = hMigrationMatrix->GetZaxis();
    xaxis->SetTitle("|R|_{Reconstructed} / (GV)");
    yaxis->SetTitle("|R|_{True} / (GV)");
    zaxis->SetTitle("Probability");
    xaxis->SetMoreLogLabels();
    yaxis->SetMoreLogLabels();
    xaxis->SetNoExponent();
    yaxis->SetNoExponent();

    zaxis->SetRangeUser(0,1);

    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    zaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.037);
    yaxis->SetTitleSize(0.037);
    zaxis->SetTitleSize(0.037);
    xaxis->SetLabelSize(0.035);
    yaxis->SetLabelSize(0.035);
    zaxis->SetLabelSize(0.035);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    zaxis->SetLabelFont(62);

    xaxis->SetLabelOffset(0);
    yaxis->SetLabelOffset(0);
    xaxis->SetTitleOffset(1.1);
    yaxis->SetTitleOffset(1.4);
    zaxis->SetTitleOffset(1.3);

    gPad->SetBottomMargin(0.13);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.2);

    gPad->Update();
    c4.Update();

    c4.SaveAs( (std::string("plot_ratio/MM_Pattern_") + pattern + std::string(".pdf")).c_str() );
}

void Plot_Statistic_Error(TH1D hStatistic_error_relative, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c5("c5","c5",1000,500);
    hStatistic_error_relative.SetFillColor(0);
    hStatistic_error_relative.SetLineColor(1);
    hStatistic_error_relative.Draw("");
    TAxis* xAxis5 = hStatistic_error_relative.GetXaxis();
    TAxis* yAxis5 = hStatistic_error_relative.GetYaxis();
    xAxis5->SetTitle("Rigidity (GV)");
    yAxis5->SetTitle("Relative Error (%)");
    xAxis5->SetMoreLogLabels();
    yAxis5->SetMoreLogLabels();
    TText *t5 = new TText(530,27,"Statistic error");
    t5->SetTextAlign(12);
    t5->SetTextColor(1);
    t5->SetTextFont(43);
    t5->SetTextSize(15);
    t5->SetTextAngle(0);
    t5->Draw();
    c5.Update();
    c5.SaveAs( (std::string("plot_ratio/Statistic_error_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}

void Plot_Statistic_Error_Absolute(TH1D hStatistic_error, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c6("c6","c6",1000,500);
    hStatistic_error.SetFillColor(0);
    hStatistic_error.SetLineColor(1);
    hStatistic_error.Draw("");
    TAxis* xAxis6 = hStatistic_error.GetXaxis();
    TAxis* yAxis6 = hStatistic_error.GetYaxis();
    xAxis6->SetTitle("Rigidity (GV)");
    yAxis6->SetTitle("Statistic Error (%)");
    xAxis6->SetMoreLogLabels();
    yAxis6->SetMoreLogLabels();
    TText *t6 = new TText(530,27,"Statistic error");
    t6->SetTextAlign(12);
    t6->SetTextColor(1);
    t6->SetTextFont(43);
    t6->SetTextSize(15);
    t6->SetTextAngle(0);
    t6->Draw();
    c6.Update();
    c6.SaveAs( (std::string("plot_ratio/Statistic_error_absolute_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}

void Plot_Statistic_Error_Compare(TH1D hStatistic_error, TH1D hPublished_Statistic_error, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c7("c7","c7",1000,500);
    hStatistic_error.SetFillColor(0);
    hStatistic_error.SetLineColor(1);
    hStatistic_error.Draw("");
    hPublished_Statistic_error.SetLineColor(4);
    hPublished_Statistic_error.Draw("same");
    TAxis* xAxis7 = hStatistic_error.GetXaxis();
    TAxis* yAxis7 = hStatistic_error.GetYaxis();
    xAxis7->SetTitle("Rigidity (GV)");
    yAxis7->SetTitle("Statistic Error");
    xAxis7->SetMoreLogLabels();
    yAxis7->SetMoreLogLabels();
    TLegend * leg7 = new TLegend(0.17,0.75,0.47,0.9);
    leg7->SetFillColor(0);
    leg7->AddEntry(&hStatistic_error,"Statistic_error (This analysis Full Span)","lp");
    leg7->AddEntry(&hPublished_Statistic_error,"Statistic_error (2016PRL)","lp");
    leg7->Draw();
    TText *t7 = new TText(530,27,"Statistic error");
    t7->SetTextAlign(12);
    t7->SetTextColor(1);
    t7->SetTextFont(43);
    t7->SetTextSize(15);
    t7->SetTextAngle(0);
    t7->Draw();
    c7.Update();
    c7.SaveAs( (std::string("plot_ratio/Statistic_error_compare_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}

void Plot_Error_Compare(TH1D hStatistic_error, TH1D hPublished_Statistic_error, TH1D hPublished_Systematic_error, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c8("c8","c8",1000,500);
    hStatistic_error.SetFillColor(0);
    hStatistic_error.SetLineColor(1);
    hStatistic_error.Draw("");
    hPublished_Statistic_error.SetLineColor(4);
    hPublished_Statistic_error.Draw("same");
    hPublished_Systematic_error.SetLineColor(6);
    hPublished_Systematic_error.Draw("same");
    TAxis* xAxis8 = hStatistic_error.GetXaxis();
    TAxis* yAxis8 = hStatistic_error.GetYaxis();
    xAxis8->SetTitle("Rigidity (GV)");
    yAxis8->SetTitle("Error");
    xAxis8->SetMoreLogLabels();
    yAxis8->SetMoreLogLabels();
    TLegend * leg8 = new TLegend(0.17,0.75,0.47,0.9);
    leg8->SetFillColor(0);
    leg8->AddEntry(&hStatistic_error,"Statistic_error(This analysis Full Span)","lp");
    leg8->AddEntry(&hPublished_Statistic_error,"Statistic_error (2016PRL)","lp");
    leg8->AddEntry(&hPublished_Systematic_error,"Systematic_error (2016PRL)","lp");
    leg8->Draw();
    TText *t8 = new TText(530,27,"Statistic error");
    t8->SetTextAlign(12);
    t8->SetTextColor(1);
    t8->SetTextFont(43);
    t8->SetTextSize(15);
    t8->SetTextAngle(0);
    t8->Draw();
    c8.Update();
    c8.SaveAs( (std::string("plot_ratio/Error_compare_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}

void Plot_Error_compare_with_total(TH1D hStatistic_error, TH1D hPublished_total_error, TH1D hPublished_Statistic_error, TH1D hPublished_Systematic_error, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c82("c82","c82",1000,500);
    hStatistic_error.SetFillColor(0);
    hStatistic_error.SetLineColor(1);
    hStatistic_error.Draw("");
    hPublished_total_error.SetLineColor(2);
    hPublished_total_error.Draw("same");
    hPublished_Statistic_error.SetLineColor(4);
    hPublished_Statistic_error.Draw("same");
    hPublished_Systematic_error.SetLineColor(6);
    hPublished_Systematic_error.Draw("same");
    TAxis* xAxis82 = hStatistic_error.GetXaxis();
    TAxis* yAxis82 = hStatistic_error.GetYaxis();
    xAxis82->SetTitle("Rigidity (GV)");
    yAxis82->SetTitle("Error");
    xAxis82->SetMoreLogLabels();
    yAxis82->SetMoreLogLabels();
    yAxis82->SetRangeUser(0,0.00006);
    yAxis82->SetLimits(0,0.00006);
    TLegend * leg82 = new TLegend(0.13,0.75,0.39,0.9);
    leg82->SetFillColor(0);
    leg82->AddEntry(&hStatistic_error,"Statistic_error(This analysis Full Span)","lp");
    leg82->AddEntry(&hPublished_Statistic_error,"Statistic_error (2016PRL)","lp");
    leg82->AddEntry(&hPublished_Systematic_error,"Systematic_error (2016PRL)","lp");
    leg82->AddEntry(&hPublished_total_error,"Total_error (2016PRL)","lp");
    leg82->Draw();
    c82.Update();
    c82.SaveAs( (std::string("plot_ratio/Error_compare_with_total_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}

void Plot_Statistic_Error_Compare_Relative(TH1D hStatistic_error_relative, TH1D hPublished_Statistic_error_relative, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c9("c9","c9",1000,500);
    hStatistic_error_relative.SetFillColor(0);
    hStatistic_error_relative.SetLineColor(2);
    hStatistic_error_relative.Draw("");
    hPublished_Statistic_error_relative.SetLineColor(1);
    hPublished_Statistic_error_relative.Draw("same");
    TAxis* xAxis9 = hStatistic_error_relative.GetXaxis();
    TAxis* yAxis9 = hStatistic_error_relative.GetYaxis();
    xAxis9->SetTitle("Rigidity (GV)");
    yAxis9->SetTitle("Statistic Error (%)");
    yAxis9->SetRangeUser(0,28);
    xAxis9->SetMoreLogLabels();
    yAxis9->SetMoreLogLabels();
    TLegend * leg9 = new TLegend(0.17,0.75,0.43,0.9);
    leg9->SetFillColor(0);
    leg9->AddEntry(&hStatistic_error_relative,"This analysis: Full Span","lp");
    leg9->AddEntry(&hPublished_Statistic_error_relative,"AMS PRL paper in 2016","lp");
    leg9->Draw();
    //TText *t9 = new TText(530,27,"Statistic error");
    //t9->SetTextAlign(12);
    //t9->SetTextColor(1);
    //t9->SetTextFont(43);
    //t9->SetTextSize(15);
    //t9->SetTextAngle(0);
    //t9->Draw();
    c9.Update();
    c9.SaveAs( (std::string("plot_ratio/Statistic_error_compare_relative_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}


void Plot_Error_Compare_Relative(TH1D hStatistic_error_relative, TH1D hPublished_Statistic_error_relative, TH1D hPublished_Systematic_error_relative, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c10("c10","c10",1000,500);
    hStatistic_error_relative.SetFillColor(0);
    hStatistic_error_relative.SetLineColor(1);
    hStatistic_error_relative.Draw("");
    hPublished_Statistic_error_relative.SetLineColor(4);
    hPublished_Statistic_error_relative.Draw("same");
    hPublished_Systematic_error_relative.SetLineColor(6);
    hPublished_Systematic_error_relative.Draw("same");
    TAxis* xAxis10 = hStatistic_error_relative.GetXaxis();
    TAxis* yAxis10 = hStatistic_error_relative.GetYaxis();
    xAxis10->SetTitle("Rigidity (GV)");
    yAxis10->SetTitle("Error (%)");
    xAxis10->SetMoreLogLabels();
    yAxis10->SetMoreLogLabels();
    TLegend * leg10 = new TLegend(0.17,0.75,0.47,0.9);
    leg10->SetFillColor(0);
    leg10->AddEntry(&hStatistic_error_relative,"Statistic_error(This analysis Full Span)","lp");
    leg10->AddEntry(&hPublished_Statistic_error_relative,"Statistic_error (2016PRL)","lp");
    leg10->AddEntry(&hPublished_Systematic_error_relative,"Statistic_error (2016PRL)","lp");
    leg10->Draw();
    TText *t10 = new TText(530,27,"Statistic error");
    t10->SetTextAlign(12);
    t10->SetTextColor(1);
    t10->SetTextFont(43);
    t10->SetTextSize(15);
    t10->SetTextAngle(0);
    t10->Draw();
    c10.Update();
    c10.SaveAs( (std::string("plot_ratio/Error_compare_relative_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}


void Plot_Error_Compare_Relative_With_Total(TH1D hStatistic_error_relative, TH1D hPublished_Statistic_error_relative, TH1D hPublished_Systematic_error_relative, TH1D hPublished_total_error_relative, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c101("c101","c101",1000,500);
    hStatistic_error_relative.SetFillColor(0);
    hStatistic_error_relative.SetLineColor(1);
    hStatistic_error_relative.Draw("");
    hPublished_Statistic_error_relative.SetLineColor(4);
    hPublished_Statistic_error_relative.Draw("same");
    hPublished_Systematic_error_relative.SetLineColor(6);
    hPublished_Systematic_error_relative.Draw("same");
    hPublished_total_error_relative.SetLineColor(2);
    hPublished_total_error_relative.Draw("same");
    TAxis* xAxis101 = hStatistic_error_relative.GetXaxis();
    TAxis* yAxis101 = hStatistic_error_relative.GetYaxis();
    xAxis101->SetTitle("Rigidity (GV)");
    yAxis101->SetTitle("Error (%)");
    xAxis101->SetMoreLogLabels();
    yAxis101->SetMoreLogLabels();
    yAxis101->SetLimits(0,42);
    yAxis101->SetRangeUser(0,42);
    TLegend * leg101 = new TLegend(0.17,0.75,0.47,0.9);
    leg101->SetFillColor(0);
    leg101->AddEntry(&hStatistic_error_relative,"Statistic_error(This analysis Full Span)","lp");
    leg101->AddEntry(&hPublished_Statistic_error_relative,"Statistic_error (2016PRL)","lp");
    leg101->AddEntry(&hPublished_Systematic_error_relative,"Statistic_error (2016PRL)","lp");
    leg101->Draw();
    c101.Update();
    c101.SaveAs( (std::string("plot_ratio/Error_compare_relative_with_total_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}


void Plot_Statistic_Error_Proportion_Published(TH1D hPublished_Statistic_error_Proportion, std::string issversion, std::string pattern, std::string NNsuffix){
    TCanvas c11("c11","c11",1000,500);
    hPublished_Statistic_error_Proportion.SetFillColor(0);
    hPublished_Statistic_error_Proportion.SetLineColor(1);
    hPublished_Statistic_error_Proportion.Draw("");
    TAxis* xAxis11 = hPublished_Statistic_error_Proportion.GetXaxis();
    TAxis* yAxis11 = hPublished_Statistic_error_Proportion.GetYaxis();
    xAxis11->SetTitle("Rigidity (GV)");
    yAxis11->SetTitle("Statistic Error Proportion");
    xAxis11->SetMoreLogLabels();
    yAxis11->SetMoreLogLabels();
    c11.Update();
    c11.SaveAs( (std::string("plot_ratio/Statistic_error_Proportion_published_Pattern_") + pattern + std::string("_") + issversion + NNsuffix + std::string(".pdf")).c_str());
}

    









    
