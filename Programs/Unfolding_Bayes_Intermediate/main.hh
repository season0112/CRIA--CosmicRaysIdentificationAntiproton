template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

void Plot_Unfolded_Ratio(TGraphErrors *gRatio_unfolded, TGraphErrors gPublishedRatio, std::string issversion){
    TCanvas c1("c1","c1",1000,500);

    gRatio_unfolded->Draw("AP *");
    gPublishedRatio.Draw("same P");

    gRatio_unfolded->SetMarkerStyle(15);
    gRatio_unfolded->SetMarkerColor(2);

    gRatio_unfolded->GetXaxis()->SetTitle("Rigidity (GV)");
    gRatio_unfolded->GetYaxis()->SetTitle("Ratio");

    c1.SaveAs( (std::string("Time_Averaged_ratio/unfolding/unfolded_ratio") + issversion + std::string(".pdf")).c_str() );
}


void Plot_EffectiveAcceptance(TGraphErrors *Effective_Acceptance){
    TCanvas c2("c2","c2",1000,500);
    Effective_Acceptance->Draw("");
    gPad->SetLogx();
    c2.SaveAs( (std::string("Time_Averaged_ratio/unfolding/Effective_Acceptance") + std::string(".pdf")).c_str() );

}


void Plot_Error(TGraph *g_error, TGraphErrors *g_SystematicError_TRD, std::string issversion){
    TCanvas c2("c2","c2",1000,500);

    g_error->Draw("AP *");
    g_SystematicError_TRD->Draw("same P");

    g_error->GetYaxis()->SetLimits(0.00000001, 0.000003);
    g_error->GetYaxis()->SetRangeUser(0.00000001, 0.000003);

    g_error->SetMarkerStyle(15);
    g_error->SetMarkerColor(2);

    g_SystematicError_TRD->SetMarkerStyle(15);
    g_SystematicError_TRD->SetMarkerColor(3);

    g_error->GetXaxis()->SetTitle("Rigidity (GV)");
    g_error->GetYaxis()->SetTitle("Error");

    c2.SaveAs( (std::string("Time_Averaged_ratio/unfolding/Error_") + issversion + std::string(".pdf")).c_str() );
}




