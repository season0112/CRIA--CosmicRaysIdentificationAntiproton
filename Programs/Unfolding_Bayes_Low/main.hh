template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


void Plot_unfolded_ratio(TGraphErrors *gRatio_unfolded, TGraphErrors gPublishedRatio, std::string issversion){

    TCanvas c1("c1","c1",1000,500);
    gRatio_unfolded->SetMarkerStyle(15);
    gRatio_unfolded->SetMarkerColor(2);

    gPublishedRatio.GetXaxis()->SetTitle("Rigidity (GV)");
    gPublishedRatio.GetYaxis()->SetTitle("Ratio");
    gRatio_unfolded->Draw("AP *");
    gRatio_unfolded->GetXaxis()->SetLimits(0.5,10);
    gRatio_unfolded->GetYaxis()->SetRangeUser(0,0.0002);
    gPublishedRatio.Draw("same P");
    c1.SaveAs( (std::string("Time_Averaged_ratio_Low/plots/unfolded_ratio") + issversion + std::string(".pdf")).c_str() );

}


void Plot_Effective_Acceptance(TGraphErrors *Effective_Acceptance){

    TCanvas c2("c2","c2",1000,500);
    Effective_Acceptance->Draw("");
    c2.SaveAs( (std::string("Time_Averaged_ratio_Low/plots/Effective_Acceptance") + std::string(".pdf")).c_str() );

}


std::tuple<TH1D *, TGraphErrors *, TGraphErrors *, TF1 *, TGraph *> CalculateEffectiveAcceptanceRatio_Parametrization(TH1D *hAcceptance_proton, TH1D *hAcceptance_antiproton, TGraphAsymmErrors *Acceptance_antiproton_all, TGraphAsymmErrors *Acceptance_proton_all, std::vector<double> subrangepointused){

    //// Calculate Effective_Acceptance ratio
    TH1D *hEffective_Acceptance = new TH1D (*hAcceptance_proton);
    hEffective_Acceptance->Divide(hAcceptance_antiproton);
    TGraphErrors *Effective_Acceptance = new TGraphErrors(hEffective_Acceptance);

    TGraphErrors *g_Effective_Acceptance_all = new TGraphErrors(Acceptance_antiproton_all->GetN()); //up to 18 GV
    double ep, Ap, ea, Aa;
    for (int acc_index = 0; acc_index < Acceptance_antiproton_all->GetN(); acc_index++){
        Aa = Acceptance_antiproton_all->GetY()[acc_index];
        Ap = Acceptance_proton_all    ->GetY()[acc_index];
        ea = Acceptance_antiproton_all->GetErrorY(acc_index);
        ep = Acceptance_proton_all->GetErrorY(acc_index);
        g_Effective_Acceptance_all->SetPoint(acc_index, Acceptance_antiproton_all->GetX()[acc_index], Ap/Aa);
        g_Effective_Acceptance_all->SetPointError(acc_index, 0, sqrt(pow(ep,2)/pow(Aa,2)+pow(Ap,2)/pow(Aa,4)*pow(ea,2)) );
    }

    //// Reset Effective_Acceptance error
    for (long unsigned int k = 1; k <= subrangepointused.size(); ++k){
        Aa = hAcceptance_antiproton->GetBinContent(k);
        Ap = hAcceptance_proton->GetBinContent(k);
        ea = hAcceptance_antiproton->GetBinError(k);
        ep = hAcceptance_proton->GetBinError(k);
        //cout<< "Aa:" << Aa << "  ea:" << ea <<endl;
        //cout<< "Ap:" << Ap << "  ep:" << ep <<endl;
        //cout<< "\n" <<endl;
        Effective_Acceptance->SetPointError(k-1, 0.0, sqrt(pow(ep,2)/pow(Aa,2)+pow(Ap,2)/pow(Aa,4)*pow(ea,2)));
    }    

    // Parametrization of ratio of effective acceptance
    Effective_Acceptance->RemovePoint(0); // remove 0.8-1.0GV point
    TF1  *function1 = new TF1("function1","[0]*log(log(x))+[1]",0,10); //to be fixed
    //TF1  *function1 = new TF1("function1","[0]*log(x)+[1]"); // to be fixed
    Effective_Acceptance->Fit(function1);
    TF1 *fittedfuction1 = Effective_Acceptance->GetFunction("function1");
    TGraph *gFitFunction = new TGraph(fittedfuction1);    

    return {hEffective_Acceptance, Effective_Acceptance, g_Effective_Acceptance_all, fittedfuction1, gFitFunction};
}



std::tuple<TGraphErrors *, TGraphErrors *, TGraphErrors *, TH1D *> CalculateRatio_And_ResetError(std::vector<double> subrangepointused, std::vector<double> subrange_low, TH1D *hantiproton, TH1D *hproton, TH1D *hAntiproton_unfolded, TH1D *hProton_unfolded, TGraph *gFitFunction, TGraphErrors *g_error){

    TH1D *ratio_raw = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    ratio_raw->Divide(hantiproton, hproton, 1, 1);
    TH1D *ratio_unfolded = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    ratio_unfolded->Divide(hAntiproton_unfolded, hProton_unfolded, 1, 1);

    TH1D *ratio_unfolded_with_effective_correction = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    TH1D *ratio_raw_with_effective_correction      = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    for (long unsigned int q = 1; q <= subrangepointused.size(); q++){
        ratio_unfolded_with_effective_correction->SetBinContent(q, ratio_unfolded->GetBinContent(q) *  gFitFunction->Eval( subrangepointused.at(q-1) ));
        ratio_raw_with_effective_correction->SetBinContent(q, ratio_raw->GetBinContent(q) *  gFitFunction->Eval( subrangepointused.at(q-1) ));
    }

    TGraphErrors *gRatio_unfolded = new TGraphErrors (ratio_unfolded_with_effective_correction);
    TGraphErrors *gRatio_Raw      = new TGraphErrors (ratio_raw_with_effective_correction);

    //// Set points error and fill StatisticalRelError
    TGraphErrors *g_StatisticalRelError = new TGraphErrors(ratio_unfolded_with_effective_correction);
    TH1D *h_StatisticalRelError         = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    for (long unsigned int q = 1; q <= subrangepointused.size(); q++){
        gRatio_unfolded->SetPointError( q-1, 0, g_error->GetY()[q-1]/100000 );
        gRatio_Raw->SetPointError( q-1, 0, g_error->GetY()[q-1]/100000 );
        g_StatisticalRelError->SetPoint (q-1, g_error->GetX()[q-1], g_error->GetY()[q-1]/100000/gRatio_unfolded->GetY()[q-1]*100); //in %.
        h_StatisticalRelError->SetBinContent(q, g_error->GetY()[q-1]/100000/gRatio_unfolded->GetY()[q-1]*100); //in %.
    }    

    return {gRatio_unfolded, gRatio_Raw, g_StatisticalRelError, h_StatisticalRelError};
}


