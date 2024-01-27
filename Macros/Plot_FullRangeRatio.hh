#include <tuple>

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

std::tuple<TGraph *> GraphDivision(TGraph *GraphA, TGraph *GraphB){
    TGraph *GraphRatio = new TGraph();
    for (int i=0; i<GraphA->GetN(); i++){
        GraphRatio->SetPoint(i, GraphA->GetX()[i], GraphA->GetY()[i]/GraphB->GetY()[i]);
    }
    return {GraphRatio};
}

//// Load

// Load Reference pbar Number
std::tuple<TGraph *, TGraph *> Load_Reference_pbar_Number(){
    std::vector<double> PublishedRatioPRL_pbarNumber (AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL_pbarNumber.begin(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL_pbarNumber.end());
    std::vector<double> bincenter (AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins());
    bincenter.erase(bincenter.begin()); // remove first bin (0.8-1.0GV)
    TGraph *g_Published_pbarNumber = new TGraph(PublishedRatioPRL_pbarNumber.size(), bincenter.data(), PublishedRatioPRL_pbarNumber.data());

    std::vector<double> PhysicsReport_pbarNumber (AntiprotonNewBinning::AntiprotonResults::PhysicsReport_pbarNumber.begin(), AntiprotonNewBinning::AntiprotonResults::PhysicsReport_pbarNumber.end());
    std::vector<double> bincenter_525version ( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_zhili525.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_zhili525.end()); // since first point (0.8-1.0) don't have result, therefore it should be removed.
    TGraph *g_PhysicsReport_pbarNumber = new TGraph(PhysicsReport_pbarNumber.size(), bincenter_525version.data(), PhysicsReport_pbarNumber.data());
    return {g_Published_pbarNumber, g_PhysicsReport_pbarNumber};
}


// Load EffectiveAcceptance correction
std::tuple<TGraphErrors *, TGraphErrors *, TGraphErrors *, TGraph *, TH1D> Load_EffectiveAcceptance_Correction(std::vector<double> PhysicsReportBinEdge){
    TFile *f_eff_B1220_antiproton         = new TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_antipr.pl1ph.021000.qgsp_bic_ams_7.8_all.root");
    TFile *f_eff_B1220_antiproton_minus10 = new TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_antipr.pl1ph.021000.qgsp_bic_ams.minus10_7.8_all.root");
    TFile *f_eff_B1220_antiproton_plus10  = new TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_antipr.pl1ph.021000.qgsp_bic_ams.plus10_7.8_all.root");
    TFile *f_eff_B1220_proton             = new TFile("/hpcwork/jara0052/sichen/Acceptance/EffectiveAcceptance_B1220_pr.pl1ph.021000_7.8_all.root");

    TGraphAsymmErrors *Acceptance_B1220_antiproton         = (TGraphAsymmErrors*)f_eff_B1220_antiproton        ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_B1220_antiproton_minus10 = (TGraphAsymmErrors*)f_eff_B1220_antiproton_minus10->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_B1220_antiproton_plus10  = (TGraphAsymmErrors*)f_eff_B1220_antiproton_plus10 ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_B1220_proton             = (TGraphAsymmErrors*)f_eff_B1220_proton            ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");

    // Remove last point to avoid NaN
    Acceptance_B1220_antiproton        ->RemovePoint(61);
    Acceptance_B1220_antiproton_minus10->RemovePoint(61);
    Acceptance_B1220_antiproton_plus10 ->RemovePoint(61);
    Acceptance_B1220_proton            ->RemovePoint(61);

    // Fill EffectiveAcceptance ratios in TGraphErrors and calculate the Systematic Uncertainties.
    int number = Acceptance_B1220_antiproton->GetN();

    TGraphErrors *Effective_Acceptance_ratio         = new TGraphErrors(number-1);
    TGraphErrors *Effective_Acceptance_ratio_minus10 = new TGraphErrors(number-1);
    TGraphErrors *Effective_Acceptance_ratio_plus10  = new TGraphErrors(number-1);
    TGraph *g_SysUncertaintyRel_ACCratio             = new TGraph(number-1);

    double ep, Ap, ea, Aa, x, Aa_minus10, Aa_plus10, ea_minus10, ea_plus10;
    for (int p = 0; p < number-1; p++) {
        //(1). Fill EffectiveAcceptance ratios
        // value taken start from 1, because 0.8-1.0GV is not used in this analysis.
        Aa         = Acceptance_B1220_antiproton        ->GetY()[p+1];
        Aa_minus10 = Acceptance_B1220_antiproton_minus10->GetY()[p+1];
        Aa_plus10  = Acceptance_B1220_antiproton_plus10 ->GetY()[p+1];
        Ap         = Acceptance_B1220_proton            ->GetY()[p+1];
        ea         = Acceptance_B1220_antiproton        ->GetErrorY(p+1);
        ea_minus10 = Acceptance_B1220_antiproton_minus10->GetErrorY(p+1);
        ea_plus10  = Acceptance_B1220_antiproton_plus10 ->GetErrorY(p+1);
        ep         = Acceptance_B1220_proton            ->GetErrorY(p+1);
         x         = Acceptance_B1220_antiproton        ->GetX()[p+1];

        Effective_Acceptance_ratio        ->SetPoint     (p, x, Ap/Aa);
        Effective_Acceptance_ratio        ->SetPointError(p, 0, sqrt(pow(ep,2)/pow(Aa,2)+pow(Ap,2)/pow(Aa,4)*pow(ea,2)) );
        Effective_Acceptance_ratio_minus10->SetPoint     (p, x, Ap/Aa_minus10);
        Effective_Acceptance_ratio_minus10->SetPointError(p, 0, sqrt(pow(ep,2)/pow(Aa_minus10,2)+pow(Ap,2)/pow(Aa_minus10,4)*pow(ea_minus10,2)) );
        Effective_Acceptance_ratio_plus10 ->SetPoint     (p, x, Ap/Aa_plus10);
        Effective_Acceptance_ratio_plus10 ->SetPointError(p, 0, sqrt(pow(ep,2)/pow(Aa_plus10,2)+pow(Ap,2)/pow(Aa_plus10,4)*pow(ea_plus10,2)) );

        //(2). calculate the Systematic Uncertainties due to ACC.
        //g_SysUncertaintyRel_ACCratio->SetPoint(p, x, ( (Ap/Aa_plus10-Ap/Aa) + (Ap/Aa-Ap/Aa_minus10) )/2 *100);   // Wrong but very close approximation. in %.
        g_SysUncertaintyRel_ACCratio->SetPoint(p, x, Ap / pow(Aa,2) * ((Aa-Aa_plus10) + (Aa_minus10-Aa))/2 *100);  // Correct Error propagation. in %.

    }

    TH1D h_SysUncertaintyRel_ACCratio = TH1D("", "", g_SysUncertaintyRel_ACCratio->GetN()-2, PhysicsReportBinEdge.data()); //Up to 525: GetN()-2
    Utilities::ConvertToHistogram ( g_SysUncertaintyRel_ACCratio, h_SysUncertaintyRel_ACCratio);
    return {Effective_Acceptance_ratio, Effective_Acceptance_ratio_minus10, Effective_Acceptance_ratio_plus10, g_SysUncertaintyRel_ACCratio, h_SysUncertaintyRel_ACCratio};
}
    
    
//// Apply Relative Uncertainty From B1220 version To B1042 version
std::tuple<TGraphErrors *, TGraphErrors *, TGraphAsymmErrors *> ApplyRelativeUncertaintyToB1042(TGraphErrors *Effective_Acceptance_ratio_B1220, TGraphErrors *Effective_Acceptance_ratio_B1220_minus10, TGraphErrors *Effective_Acceptance_ratio_B1220_plus10, TGraphErrors *Effective_Acceptance_High_B1042, TGraphErrors *Effective_Acceptance_Low_B1042, TGraphErrors *Effective_Acceptance_Intermediate_B1042, TGraphErrors *Effective_Acceptance_Intermediate_inHigh_B1042){

    int TotalPointsNumber = Effective_Acceptance_Low_B1042->GetN() + Effective_Acceptance_Intermediate_B1042->GetN() + Effective_Acceptance_High_B1042->GetN();

    TGraphErrors *Effective_Acceptance_ratio_B1042_minus10_AllRange = new TGraphErrors(TotalPointsNumber);
    TGraphErrors *Effective_Acceptance_ratio_B1042_plus10_AllRange  = new TGraphErrors(TotalPointsNumber);
    TGraphAsymmErrors *Effective_Acceptance_ratio_B1042_AllRange    = new TGraphAsymmErrors(TotalPointsNumber);

    // Low Range
    for (int i=0; i<Effective_Acceptance_Low_B1042->GetN(); i++ ){
        double factor_minus10 = Effective_Acceptance_ratio_B1220_minus10->GetY()[i+1] / Effective_Acceptance_ratio_B1220->GetY()[i+1];
        Effective_Acceptance_ratio_B1042_minus10_AllRange->SetPoint     (i, Effective_Acceptance_Low_B1042->GetX()[i], Effective_Acceptance_Low_B1042->GetY()[i] * factor_minus10);
        Effective_Acceptance_ratio_B1042_minus10_AllRange->SetPointError(i, 0, 0);

        double factor_plus10 = Effective_Acceptance_ratio_B1220_plus10->GetY()[i+1] / Effective_Acceptance_ratio_B1220->GetY()[i+1];
        Effective_Acceptance_ratio_B1042_plus10_AllRange->SetPoint     (i, Effective_Acceptance_Low_B1042->GetX()[i], Effective_Acceptance_Low_B1042->GetY()[i] * factor_plus10); 
        Effective_Acceptance_ratio_B1042_plus10_AllRange->SetPointError(i, 0, 0);

        Effective_Acceptance_ratio_B1042_AllRange->SetPoint     (i, Effective_Acceptance_Low_B1042->GetX()[i], Effective_Acceptance_Low_B1042->GetY()[i]);
        Effective_Acceptance_ratio_B1042_AllRange->SetPointError(i, 0, 0,
                                                                 Effective_Acceptance_ratio_B1042_AllRange->GetY()[i]-Effective_Acceptance_ratio_B1042_minus10_AllRange->GetY()[i],
                                                                 Effective_Acceptance_ratio_B1042_plus10_AllRange->GetY()[i]-Effective_Acceptance_ratio_B1042_AllRange->GetY()[i]); // (exl, exh, eyl. eyh)
        //cout<< Effective_Acceptance_ratio_B1042_plus10_AllRange->GetY()[i]-Effective_Acceptance_Low_B1042->GetY()[i] << endl;

        cout<< Effective_Acceptance_ratio_B1220_plus10->GetX()[i+1] << " and " << Effective_Acceptance_ratio_B1220_minus10->GetX()[i+1] << " and " << Effective_Acceptance_Low_B1042->GetX()[i] << endl;

    }

    // Intermediate Range
    for (int i=0+Effective_Acceptance_Low_B1042->GetN(); i<Effective_Acceptance_Low_B1042->GetN()+Effective_Acceptance_Intermediate_B1042->GetN(); i++ ){

        double factor_minus10 = Effective_Acceptance_ratio_B1220_minus10->GetY()[i+1] / Effective_Acceptance_ratio_B1220->GetY()[i+1];
        Effective_Acceptance_ratio_B1042_minus10_AllRange->SetPoint     (i, Effective_Acceptance_Intermediate_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()], Effective_Acceptance_Intermediate_B1042->GetY()[i-Effective_Acceptance_Low_B1042->GetN()] * factor_minus10 );
        Effective_Acceptance_ratio_B1042_minus10_AllRange->SetPointError(i, 0, 0);

        double factor_plus10 = Effective_Acceptance_ratio_B1220_plus10->GetY()[i+1] / Effective_Acceptance_ratio_B1220->GetY()[i+1];
        Effective_Acceptance_ratio_B1042_plus10_AllRange->SetPoint      (i, Effective_Acceptance_Intermediate_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()], Effective_Acceptance_Intermediate_B1042->GetY()[i-Effective_Acceptance_Low_B1042->GetN()] * factor_plus10 );
        Effective_Acceptance_ratio_B1042_plus10_AllRange->SetPointError (i, 0, 0);

        Effective_Acceptance_ratio_B1042_AllRange->SetPoint     (i, Effective_Acceptance_Intermediate_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()], Effective_Acceptance_Intermediate_B1042->GetY()[i-Effective_Acceptance_Low_B1042->GetN()]);
        Effective_Acceptance_ratio_B1042_AllRange->SetPointError(i, 0, 0,
                                                                 Effective_Acceptance_ratio_B1042_AllRange->GetY()[i]-Effective_Acceptance_ratio_B1042_minus10_AllRange->GetY()[i],
                                                                 Effective_Acceptance_ratio_B1042_plus10_AllRange->GetY()[i]-Effective_Acceptance_ratio_B1042_AllRange->GetY()[i]); // (exl, exh, eyl. eyh)

        cout<< Effective_Acceptance_ratio_B1220_plus10->GetX()[i+1] << " and " << Effective_Acceptance_ratio_B1220_minus10->GetX()[i+1] << " and " << Effective_Acceptance_Intermediate_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()] << endl;

    }

    // High Range
    for (int i=0+Effective_Acceptance_Low_B1042->GetN()+Effective_Acceptance_Intermediate_B1042->GetN(); i<Effective_Acceptance_Low_B1042->GetN()+Effective_Acceptance_Intermediate_B1042->GetN()+Effective_Acceptance_High_B1042->GetN(); i++ ){

        double factor_minus10 = Effective_Acceptance_ratio_B1220_minus10->GetY()[i+1] / Effective_Acceptance_ratio_B1220->GetY()[i+1];
        Effective_Acceptance_ratio_B1042_minus10_AllRange->SetPoint     (i, Effective_Acceptance_High_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()-Effective_Acceptance_Intermediate_B1042->GetN()], Effective_Acceptance_High_B1042->GetY()[i-Effective_Acceptance_Low_B1042->GetN()-Effective_Acceptance_Intermediate_B1042->GetN()] * factor_minus10 );
        Effective_Acceptance_ratio_B1042_minus10_AllRange->SetPointError(i, 0, 0);

        double factor_plus10 = Effective_Acceptance_ratio_B1220_plus10->GetY()[i+1] / Effective_Acceptance_ratio_B1220->GetY()[i+1];
        Effective_Acceptance_ratio_B1042_plus10_AllRange->SetPoint      (i, Effective_Acceptance_High_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()-Effective_Acceptance_Intermediate_B1042->GetN()], Effective_Acceptance_High_B1042->GetY()[i-Effective_Acceptance_Low_B1042->GetN()-Effective_Acceptance_Intermediate_B1042->GetN()] * factor_plus10);
        Effective_Acceptance_ratio_B1042_plus10_AllRange->SetPointError (i, 0, 0);

        Effective_Acceptance_ratio_B1042_AllRange->SetPoint     (i, Effective_Acceptance_High_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()-Effective_Acceptance_Intermediate_B1042->GetN()], Effective_Acceptance_High_B1042->GetY()[i-Effective_Acceptance_Low_B1042->GetN()-Effective_Acceptance_Intermediate_B1042->GetN()]);
        Effective_Acceptance_ratio_B1042_AllRange->SetPointError(i, 0, 0,
                                                                 Effective_Acceptance_ratio_B1042_AllRange->GetY()[i]-Effective_Acceptance_ratio_B1042_minus10_AllRange->GetY()[i],
                                                                 Effective_Acceptance_ratio_B1042_plus10_AllRange->GetY()[i]-Effective_Acceptance_ratio_B1042_AllRange->GetY()[i]); // (exl, exh, eyl. eyh)

        cout<< Effective_Acceptance_ratio_B1220_plus10->GetX()[i+1] << " and " << Effective_Acceptance_ratio_B1220_minus10->GetX()[i+1] << " and " << Effective_Acceptance_High_B1042->GetX()[i-Effective_Acceptance_Low_B1042->GetN()-Effective_Acceptance_Intermediate_B1042->GetN()] << endl;

    }

    return {Effective_Acceptance_ratio_B1042_minus10_AllRange, Effective_Acceptance_ratio_B1042_plus10_AllRange, Effective_Acceptance_ratio_B1042_AllRange};
}


//// Calculate "Total error" and "Total Systematic error" in This analysis
void Calculate_TotalError_TotalSystematicError_LowRange(TGraphErrors *g_StatisticalError_Low, TGraph *g_SysUncertaintyRel_ACCratio, TGraph *g_System_ACC_Low, TGraphErrors *g_SystematicError_Shape_Low, TGraphErrors *g_LowResult, TGraph *g_TotalSysError_Low, TGraph *g_TotalSysRelError_Low, TH1D *h_TotalSysRelError_Low, TGraph *g_TotalError_Low, TGraph *g_TotalRelError_Low, TH1D *h_TotalRelError_Low, TH1D *h_SystematicError_Shape_Low){
    // Low
    for (int i=0; i<g_StatisticalError_Low->GetN(); i++){
        // Check
        if (to_string_with_precision(g_SysUncertaintyRel_ACCratio->GetX()[i]) != to_string_with_precision(g_StatisticalError_Low->GetX()[i])){
            cout<< "!!! Attention:" << g_SysUncertaintyRel_ACCratio->GetX()[i] << "!=" << g_StatisticalError_Low->GetX()[i] <<endl;}
       
        // SystematicRelativeError_fromShape
        h_SystematicError_Shape_Low->SetBinContent(i+1, g_SystematicError_Shape_Low->GetY()[i] / g_LowResult->GetY()[i] * 100);
        cout<< "checking:::::::::" << g_SystematicError_Shape_Low->GetY()[i] / g_LowResult->GetY()[i] * 100 <<endl; 

        // Total systematic error            
        g_System_ACC_Low   ->SetPoint(i, g_StatisticalError_Low->GetX()[i], pow( (pow(g_SysUncertaintyRel_ACCratio->GetY()[i]/100*g_LowResult->GetY()[i], 2)), 0.5));
        g_TotalSysError_Low->SetPoint(i, g_StatisticalError_Low->GetX()[i], pow( (pow(g_System_ACC_Low->GetY()[i], 2) + pow(g_SystematicError_Shape_Low->GetY()[i] ,2)), 0.5) );

        g_TotalSysRelError_Low->SetPoint(i, g_StatisticalError_Low->GetX()[i], g_TotalSysError_Low->GetY()[i]/g_LowResult->GetY()[i] * 100);
        h_TotalSysRelError_Low->SetBinContent(i+1, g_TotalSysRelError_Low->GetY()[i]);

        // Total error
        g_TotalError_Low->SetPoint(i, g_StatisticalError_Low->GetX()[i],    pow( (pow(g_TotalSysError_Low->GetY()[i], 2) + pow(g_StatisticalError_Low->GetY()[i],2)), 0.5));
        g_TotalRelError_Low->SetPoint(i, g_StatisticalError_Low->GetX()[i], pow( (pow(g_TotalSysError_Low->GetY()[i], 2) + pow(g_StatisticalError_Low->GetY()[i],2)), 0.5)/g_LowResult->GetY()[i] * 100 );
        h_TotalRelError_Low->SetBinContent(i+1, g_TotalRelError_Low->GetY()[i]);
    }
}


void Calculate_TotalError_TotalSystematicError_IntermediateRange(TGraph *g_StatisticError_Intermediate, TGraph *g_SysUncertaintyRel_ACCratio, TGraphErrors *g_StatisticalError_Low, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, TGraph *g_System_ACC_Intermediate, TGraphErrors *g_IntermediateResult, TGraph *g_TotalSysError_Intermediate, TGraph *g_TotalSysRelError_Intermediate, TH1D *h_TotalSysRelError_Intermediate, TGraph *g_TotalError_Intermediate, TGraph *g_TotalRelError_Intermediate, TH1D *h_TotalRelError_Intermediate, TGraphErrors *g_SystematicError_TRD_Intermediate, TH1D *h_SystematicRelativeError_TRD){
    // Intermediate
    for (int i=0; i<g_StatisticError_Intermediate->GetN(); i++){
        // Check
        if (to_string_with_precision(g_SysUncertaintyRel_ACCratio->GetX()[i+g_StatisticalError_Low->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin]) != to_string_with_precision(g_StatisticError_Intermediate->GetX()[i])){
            cout<< "!!! Attention:" << g_SysUncertaintyRel_ACCratio->GetX()[i+g_StatisticalError_Low->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin] << "=?" << g_StatisticError_Intermediate->GetX()[i] <<endl;}

        // Fill h_SystematicRelativeError_TRD
        h_SystematicRelativeError_TRD ->SetBinContent(i+1, g_SystematicError_TRD_Intermediate->GetY()[i] / g_IntermediateResult->GetY()[i] * 100 ); 

        // Total Systematic Error
        g_System_ACC_Intermediate->SetPoint(i,        g_StatisticError_Intermediate->GetX()[i], pow( (pow(g_SysUncertaintyRel_ACCratio->GetY()[i+g_StatisticalError_Low->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin]/100*g_IntermediateResult->GetY()[i], 2) ), 0.5));
        g_TotalSysError_Intermediate->SetPoint(i,     g_StatisticError_Intermediate->GetX()[i], pow( (pow(g_System_ACC_Intermediate->GetY()[i], 2) + pow(g_SystematicError_TRD_Intermediate->GetY()[i] ,2)), 0.5) );

        g_TotalSysRelError_Intermediate->SetPoint(i,  g_StatisticError_Intermediate->GetX()[i], g_TotalSysError_Intermediate->GetY()[i]/g_IntermediateResult->GetY()[i] * 100 );
        h_TotalSysRelError_Intermediate->SetBinContent(i+1, g_TotalSysRelError_Intermediate->GetY()[i]);         

        // Total Error
        g_TotalError_Intermediate->SetPoint(i,    g_StatisticError_Intermediate->GetX()[i], pow(  (pow(g_TotalSysError_Intermediate->GetY()[i], 2) + pow(g_StatisticError_Intermediate->GetY()[i],2)), 0.5));
        g_TotalRelError_Intermediate->SetPoint(i, g_StatisticError_Intermediate->GetX()[i], pow(  (pow(g_TotalSysError_Intermediate->GetY()[i], 2) + pow(g_StatisticError_Intermediate->GetY()[i],2)), 0.5)/g_IntermediateResult->GetY()[i] * 100 );
        h_TotalRelError_Intermediate->SetBinContent(i+1, g_TotalRelError_Intermediate->GetY()[i]);
    }
}


void Calculate_TotalError_TotalSystematicError_HighRange(TGraph *g_StatisticalError_High, TGraph *g_SysUncertaintyRel_ACCratio_B1220, TGraphErrors *g_StatisticalError_Low, TGraph *g_StatisticError_Intermediate, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin, TGraph *g_System_ACC_High, TGraph *g_TotalSysError_High, TGraphErrors *g_HighResult, TGraph *g_System_CC_High, TGraph *g_TotalSysRelError_High, TGraph *g_System_CC_relative_High, TH1D *h_TotalSysRelError_High, TH1D *h_System_CC_relative_High, TGraph *g_TotalError_High, TGraph *g_TotalRelError_High, TH1D *h_TotalRelError_High){
    // High
    for (int i=0; i<g_StatisticalError_High->GetN(); i++){
        if (to_string_with_precision(g_SysUncertaintyRel_ACCratio_B1220->GetX()[i+g_StatisticalError_Low->GetN()+g_StatisticError_Intermediate->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin-IntermediateRangeRemoveAtEnd-HighRangeRemovedAtBegin]) != to_string_with_precision(g_StatisticalError_High->GetX()[i])){
            cout<< "!!! Attention:" << g_SysUncertaintyRel_ACCratio_B1220->GetX()[i+g_StatisticalError_Low->GetN()+g_StatisticError_Intermediate->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin-IntermediateRangeRemoveAtEnd-HighRangeRemovedAtBegin] << "=?" << g_StatisticalError_High->GetX()[i] <<endl;}
        // total systematic error
        g_System_ACC_High->SetPoint(i,       g_StatisticalError_High->GetX()[i],  g_SysUncertaintyRel_ACCratio_B1220->GetY()[i+g_StatisticalError_Low->GetN()+g_StatisticError_Intermediate->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin-IntermediateRangeRemoveAtEnd-HighRangeRemovedAtBegin]/100*g_HighResult->GetY()[i] );
        g_System_CC_relative_High->SetPoint(i, g_StatisticalError_High->GetX()[i], g_System_CC_High->GetY()[i]/g_HighResult->GetY()[i] * 100);
        g_TotalSysError_High->SetPoint(i,    g_StatisticalError_High->GetX()[i],  pow( (pow(g_SysUncertaintyRel_ACCratio_B1220->GetY()[i+g_StatisticalError_Low->GetN()+g_StatisticError_Intermediate->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin-IntermediateRangeRemoveAtEnd-HighRangeRemovedAtBegin]/100*g_HighResult->GetY()[i], 2) + pow(g_System_CC_High->GetY()[i],2)) , 0.5));
        g_TotalSysRelError_High->SetPoint(i, g_StatisticalError_High->GetX()[i],  pow( (pow(g_SysUncertaintyRel_ACCratio_B1220->GetY()[i+g_StatisticalError_Low->GetN()+g_StatisticError_Intermediate->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin-IntermediateRangeRemoveAtEnd-HighRangeRemovedAtBegin]/100*g_HighResult->GetY()[i], 2) + pow(g_System_CC_High->GetY()[i],2)) , 0.5)/g_HighResult->GetY()[i] * 100);
        h_TotalSysRelError_High->SetBinContent(i+1, g_TotalSysRelError_High->GetY()[i]);
        h_System_CC_relative_High->SetBinContent(i+1, g_System_CC_relative_High->GetY()[i]);
        // total error
        g_TotalError_High->SetPoint(i,    g_StatisticalError_High->GetX()[i],   pow( (pow(g_TotalSysError_High->GetY()[i], 2) + pow(g_StatisticalError_High->GetY()[i],2)), 0.5));
        g_TotalRelError_High->SetPoint(i, g_StatisticalError_High->GetX()[i],   pow( (pow(g_TotalSysError_High->GetY()[i], 2) + pow(g_StatisticalError_High->GetY()[i],2)), 0.5)/g_HighResult->GetY()[i] * 100);
        h_TotalRelError_High->SetBinContent(i+1, g_TotalRelError_High->GetY()[i]);
    }
}


//// Binning shift.
void BinningShift(double binshift, TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, TGraphErrors *g_HighResult){
    for (int i=0; i<g_LowResult->GetN(); i++){
        g_LowResult->SetPoint(i, g_LowResult->GetX()[i] * binshift, g_LowResult->GetY()[i]);
    }
    for (int i=0; i<g_IntermediateResult->GetN(); i++){
        g_IntermediateResult->SetPoint(i, g_IntermediateResult->GetX()[i] * binshift, g_IntermediateResult->GetY()[i]);
    }
    for (int i=0; i<g_HighResult->GetN(); i++){
        g_HighResult->SetPoint(i, g_HighResult->GetX()[i] * binshift, g_HighResult->GetY()[i]);
    }
}


//// Reset low,intermedaite,high ratio graph error to total errot from statistical error.
void Reset_Error_From_Statistical_To_Total(TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, TGraphErrors *g_HighResult, TGraph *g_TotalError_Low, TGraph *g_TotalError_Intermediate, TGraph *g_TotalError_High){
    for (int i=0; i<g_LowResult->GetN(); i++){
        g_LowResult->SetPointError(i, g_LowResult->GetErrorX(i), g_TotalError_Low->GetY()[i]);
    }
    for (int i=0; i<g_IntermediateResult->GetN(); i++){
        g_IntermediateResult->SetPointError(i, g_IntermediateResult->GetErrorX(i), g_TotalError_Intermediate->GetY()[i]);
    }
    for (int i=0; i<g_HighResult->GetN(); i++){
        g_HighResult->SetPointError(i, g_HighResult->GetErrorX(i), g_TotalError_High->GetY()[i]);
    }
}


//// Make High Result With Scaled AcceptanceError
void MakeResultWithScaledAcceptanceError(TGraphErrors *g_LowResult_ScaledAcceptanceError, TGraphErrors *g_IntermediateResult_ScaledAcceptanceError, TGraphErrors *g_HighResult_ScaledAcceptanceError, TGraph *g_StatisticalError_Low, TGraph *g_StatisticError_Intermediate, TGraph *g_StatisticalError_High, TGraph *g_SysUncertaintyRel_ACCratio, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin, TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, TGraphErrors *g_HighResult, TGraph *g_System_CC_High, TGraph *g_System_ACC_Low, TGraph *g_System_ACC_Intermediate, TGraph *g_System_ACC_High, TGraphErrors *g_SystematicError_Shape_Low, TGraphErrors *g_SystematicError_TRD_Intermediate, double Scaler){

    // Low
    for (int i=0; i<g_LowResult_ScaledAcceptanceError->GetN(); i++){
        g_LowResult_ScaledAcceptanceError->SetPoint(i, g_LowResult->GetX()[i], g_LowResult->GetY()[i]);
        double StatisticError              = g_StatisticalError_Low->GetY()[i];
        double TotalSystematicError_scaled = pow( (pow(g_System_ACC_Low->GetY()[i] * Scaler, 2) + pow(g_SystematicError_Shape_Low->GetY()[i] ,2)), 0.5);
        double TotalError_scaled           = pow( (pow(TotalSystematicError_scaled, 2) + pow(StatisticError,2)), 0.5);
        g_LowResult_ScaledAcceptanceError->SetPointError(i, g_LowResult_ScaledAcceptanceError->GetErrorX(i), TotalError_scaled); 
    }

    // Intermediate
    for (int i=0; i<g_StatisticError_Intermediate->GetN(); i++){
        g_IntermediateResult_ScaledAcceptanceError->SetPoint(i, g_IntermediateResult->GetX()[i], g_IntermediateResult->GetY()[i]); 
        double StatisticError              = g_StatisticError_Intermediate->GetY()[i];         
        double TotalSystematicError_scaled = pow( (pow(g_System_ACC_Intermediate->GetY()[i] * Scaler, 2) + pow(g_SystematicError_TRD_Intermediate->GetY()[i] ,2)), 0.5); 
        double TotalError_scaled           = pow( (pow(TotalSystematicError_scaled, 2) + pow(StatisticError,2)), 0.5);
        g_IntermediateResult_ScaledAcceptanceError->SetPointError(i, g_IntermediateResult_ScaledAcceptanceError->GetErrorX(i), TotalError_scaled); 
    }

    // High
    for (int i=0; i<g_HighResult_ScaledAcceptanceError->GetN(); i++){
        g_HighResult_ScaledAcceptanceError->SetPoint(i, g_HighResult->GetX()[i], g_HighResult->GetY()[i]); 
        double StatisticError              = g_StatisticalError_High->GetY()[i];
        double TotalSystematicError_scaled = pow( (pow(g_SysUncertaintyRel_ACCratio->GetY()[i+g_StatisticalError_Low->GetN()+g_StatisticError_Intermediate->GetN()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin-IntermediateRangeRemoveAtEnd-HighRangeRemovedAtBegin]/100*g_HighResult->GetY()[i] * Scaler, 2) + pow(g_System_CC_High->GetY()[i],2)) , 0.5);
        //double TotalSystematicError_scaled = pow( (pow(g_System_ACC_High->GetY()[i] * Scaler, 2) + pow(g_System_CC_High->GetY()[i],2)) , 0.5);
        double TotalError_scaled           = pow( (pow(TotalSystematicError_scaled, 2) + pow(StatisticError,2)), 0.5);
        g_HighResult_ScaledAcceptanceError->SetPointError(i, g_HighResult_ScaledAcceptanceError->GetErrorX(i), TotalError_scaled); 
    }
}


//// Overlap Range
// Deal with overlap range for Relative Error plot
void Reset_Error_For_OverlapedRange(TH1D *h_StatisticalRelError_Low, TH1D *h_StatisticalRelError_Intermediate, TH1D *h_TotalRelError_Low, TH1D *h_TotalRelError_Intermediate, TH1D *h_TotalSysRelError_Low, TH1D *h_TotalSysRelError_Intermediate, TH1D *h_Statistic_Error_Relative_High, TH1D *h_TotalRelError_High, TH1D *h_TotalSysRelError_High, TH1D *h_System_CC_relative_High, TH1D *h_SystematicError_Shape_Low, TH1D *h_SystematicRelativeError_TRD, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin){

    for (int i=0; i<LowRangeRemoveAtEnd; i++){
        h_StatisticalRelError_Low->SetBinContent( h_StatisticalRelError_Low->GetNbinsX()-LowRangeRemoveAtEnd+i, h_StatisticalRelError_Intermediate->GetBinContent(IntermediateRangeRemoveAtBegin+i) );
        h_TotalRelError_Low      ->SetBinContent( h_TotalRelError_Low->GetNbinsX()-LowRangeRemoveAtEnd+i, h_TotalRelError_Intermediate->GetBinContent(IntermediateRangeRemoveAtBegin+i) );
        h_TotalSysRelError_Low   ->SetBinContent( h_TotalSysRelError_Low->GetNbinsX()-LowRangeRemoveAtEnd+i, h_TotalSysRelError_Intermediate->GetBinContent(IntermediateRangeRemoveAtBegin+i) );
        h_SystematicError_Shape_Low->SetBinContent( h_SystematicError_Shape_Low->GetNbinsX()-LowRangeRemoveAtEnd+i, h_SystematicRelativeError_TRD->GetBinContent(IntermediateRangeRemoveAtBegin+i) );
    }
    for (int i=0; i<IntermediateRangeRemoveAtBegin; i++){
        h_StatisticalRelError_Intermediate->SetBinContent( 1+i, h_StatisticalRelError_Low->GetBinContent( h_StatisticalRelError_Low->GetNbinsX()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin+1+i ) );
        h_TotalRelError_Intermediate      ->SetBinContent( 1+i, h_TotalRelError_Low->GetBinContent( h_TotalRelError_Low->GetNbinsX()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin+1+i ) );
        h_TotalSysRelError_Intermediate   ->SetBinContent( 1+i, h_TotalSysRelError_Low->GetBinContent( h_TotalSysRelError_Low->GetNbinsX()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin+1+i ) );
        h_SystematicRelativeError_TRD     ->SetBinContent( 1+i, h_SystematicError_Shape_Low->GetBinContent( h_SystematicError_Shape_Low->GetNbinsX()-LowRangeRemoveAtEnd-IntermediateRangeRemoveAtBegin+1+i ) );
    }
    for (int i=0; i<IntermediateRangeRemoveAtEnd; i++){
        h_StatisticalRelError_Intermediate->SetBinContent( h_StatisticalRelError_Intermediate->GetNbinsX()-i, h_Statistic_Error_Relative_High->GetBinContent(IntermediateRangeRemoveAtEnd+HighRangeRemovedAtBegin-i));
        h_TotalRelError_Intermediate      ->SetBinContent( h_TotalRelError_Intermediate->GetNbinsX()-i, h_TotalRelError_High->GetBinContent( IntermediateRangeRemoveAtEnd+HighRangeRemovedAtBegin-i));
        h_TotalSysRelError_Intermediate   ->SetBinContent( h_TotalSysRelError_Intermediate->GetNbinsX()-i, h_TotalSysRelError_High->GetBinContent( IntermediateRangeRemoveAtEnd+HighRangeRemovedAtBegin-i));
    }
    for (int i=0; i<HighRangeRemovedAtBegin; i++){
        h_Statistic_Error_Relative_High->SetBinContent( 1+i, h_StatisticalRelError_Intermediate->GetBinContent( h_StatisticalRelError_Intermediate->GetNbinsX()-2-i));
        h_TotalRelError_High           ->SetBinContent( 1+i, h_TotalRelError_Intermediate->GetBinContent( h_TotalRelError_Intermediate->GetNbinsX()-2-i));
        h_System_CC_relative_High->SetBinContent( 1+i, 0 ); // assume first a few points in high range has 0 CC uncertainty.
        h_TotalSysRelError_High  ->SetBinContent( 1+i, h_TotalSysRelError_Intermediate->GetBinContent( h_TotalSysRelError_Intermediate->GetNbinsX()-2-i));
    }
}


// Remove overlaped range for Antiproton and Proton Numbers and Acceptance.
void Remove_OverlapedRange_Antiproton_Proton_Numbers_and_Acceptance(
TGraph *g_Antiproton_number_unfolded_Low         , TGraph *g_Proton_number_unfolded_Low         , TGraphErrors *Effective_Acceptance_Low, 
TGraph *g_Antiproton_number_unfolded_Intermediate, TGraph *g_Proton_number_unfolded_Intermediate, TGraphErrors *Effective_Acceptance_Intermediate, 
TGraph *g_Antiproton_number_unfolded_High        , TGraph *g_Proton_number_unfolded_High        , TGraphErrors *Effective_Acceptance_Intermediate_inHigh,
TGraph *g_Antiproton_number_unfolded_High_P0VGGNN, TGraph *g_Antiproton_number_unfolded_High_P0 , TGraph *g_Antiproton_number_unfolded_High_P1, TGraph *g_Antiproton_number_unfolded_High_P2, TGraph *g_Antiproton_number_unfolded_High_P4,
TGraph *g_Proton_number_unfolded_High_P0VGGNN    , TGraph *g_Proton_number_unfolded_High_P0     , TGraph *g_Proton_number_unfolded_High_P1    , TGraph *g_Proton_number_unfolded_High_P2    , TGraph *g_Proton_number_unfolded_High_P4,
TGraph *g_Antiproton_number_raw_High_P0VGGNN     , TGraph *g_Antiproton_number_raw_High_P0      , TGraph *g_Antiproton_number_raw_High_P1     , TGraph *g_Antiproton_number_raw_High_P2     , TGraph *g_Antiproton_number_raw_High_P4,
TGraph *g_Proton_number_raw_High_P0VGGNN         , TGraph *g_Proton_number_raw_High_P0          , TGraph *g_Proton_number_raw_High_P1         , TGraph *g_Proton_number_raw_High_P2         , TGraph *g_Proton_number_raw_High_P4,
int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin, TGraph *Effective_Acceptance_High, 
TGraph *g_antiproton_number_raw_Low, TGraph *g_Antiproton_number_raw_Intermediate, TGraph *g_Antiproton_number_raw_High,
TGraph *g_Proton_number_raw_Low, TGraph *g_Proton_number_raw_Intermediate, TGraph *g_Proton_number_raw_High){

    for (int i=0; i<LowRangeRemoveAtEnd; i++){
        g_antiproton_number_raw_Low     ->RemovePoint(g_antiproton_number_raw_Low->GetN()-1);
        g_Antiproton_number_unfolded_Low->RemovePoint(g_Antiproton_number_unfolded_Low->GetN()-1);
        g_Proton_number_unfolded_Low    ->RemovePoint(g_Proton_number_unfolded_Low->GetN()-1);
        g_Proton_number_raw_Low         ->RemovePoint(g_Proton_number_raw_Low->GetN()-1);
        Effective_Acceptance_Low        ->RemovePoint(Effective_Acceptance_Low->GetN()-1);
    }
    for (int i=0; i<IntermediateRangeRemoveAtBegin; i++){
        g_Antiproton_number_raw_Intermediate     ->RemovePoint(0);
        g_Antiproton_number_unfolded_Intermediate->RemovePoint(0);
        g_Proton_number_unfolded_Intermediate    ->RemovePoint(0);
        g_Proton_number_raw_Intermediate         ->RemovePoint(0);
        Effective_Acceptance_Intermediate        ->RemovePoint(0);
    }
    for (int i=0; i<IntermediateRangeRemoveAtEnd; i++){
        g_Antiproton_number_raw_Intermediate     ->RemovePoint(g_Antiproton_number_raw_Intermediate->GetN()-1);
        g_Antiproton_number_unfolded_Intermediate->RemovePoint(g_Antiproton_number_unfolded_Intermediate->GetN()-1);
        g_Proton_number_unfolded_Intermediate    ->RemovePoint(g_Proton_number_unfolded_Intermediate->GetN()-1);
        g_Proton_number_raw_Intermediate         ->RemovePoint(g_Proton_number_raw_Intermediate->GetN()-1); 
        Effective_Acceptance_Intermediate        ->RemovePoint(Effective_Acceptance_Intermediate->GetN()-1);
    }
    for (int i=0; i<HighRangeRemovedAtBegin; i++){
        g_Antiproton_number_raw_High             ->RemovePoint(0);
        g_Antiproton_number_unfolded_High        ->RemovePoint(0);

        g_Antiproton_number_unfolded_High_P0VGGNN->RemovePoint(0);
        g_Antiproton_number_unfolded_High_P0     ->RemovePoint(0);
        g_Antiproton_number_unfolded_High_P1     ->RemovePoint(0);
        g_Antiproton_number_unfolded_High_P2     ->RemovePoint(0);
        g_Antiproton_number_unfolded_High_P4     ->RemovePoint(0);

        g_Antiproton_number_raw_High_P0VGGNN     ->RemovePoint(0);
        g_Antiproton_number_raw_High_P0          ->RemovePoint(0);
        g_Antiproton_number_raw_High_P1          ->RemovePoint(0);
        g_Antiproton_number_raw_High_P2          ->RemovePoint(0);
        g_Antiproton_number_raw_High_P4          ->RemovePoint(0);

        g_Proton_number_unfolded_High            ->RemovePoint(0);
        g_Proton_number_raw_High                 ->RemovePoint(0);

        g_Proton_number_unfolded_High_P0VGGNN    ->RemovePoint(0);
        g_Proton_number_unfolded_High_P0         ->RemovePoint(0);
        g_Proton_number_unfolded_High_P1         ->RemovePoint(0);
        g_Proton_number_unfolded_High_P2         ->RemovePoint(0);
        g_Proton_number_unfolded_High_P4         ->RemovePoint(0);

        g_Proton_number_raw_High_P0VGGNN         ->RemovePoint(0);   
        g_Proton_number_raw_High_P0              ->RemovePoint(0);
        g_Proton_number_raw_High_P1              ->RemovePoint(0);
        g_Proton_number_raw_High_P2              ->RemovePoint(0);
        g_Proton_number_raw_High_P4              ->RemovePoint(0);

        Effective_Acceptance_Intermediate_inHigh ->RemovePoint(0);
        Effective_Acceptance_High                ->RemovePoint(0);
    }
}


// Remove overlaped range for Ratio.
void Remove_OverlapedRange_Ratio(TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, TGraphErrors *g_HighResult, TGraphErrors *g_LowResult_ScaledAcceptanceError, TGraphErrors *g_IntermediateResult_ScaledAcceptanceError, TGraphErrors *g_HighResult_ScaledAcceptanceError, TGraph *g_System_ACC_Low, TGraph *g_System_ACC_Intermediate, TGraph *g_System_ACC_High, TGraph *g_TotalSysError_Low, TGraph *g_TotalSysError_Intermediate, TGraph *g_TotalSysError_High, TGraphErrors *g_StatisticalError_Low, TGraph *g_StatisticError_Intermediate, TGraph *g_StatisticalError_High, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin){

    cout << "test:" <<endl;
    cout<< g_System_ACC_Low->GetN() << endl;
    cout<< g_LowResult->GetN() << endl;
    cout<< g_TotalSysError_Low->GetN() << endl;

    cout<< g_System_ACC_Intermediate->GetN() << endl;
    cout<< g_IntermediateResult->GetN() << endl;
    cout<< g_TotalSysError_Intermediate->GetN() << endl;

    cout<< g_System_ACC_High->GetN() << endl;
    cout<< g_HighResult->GetN() << endl;
    cout<< g_TotalSysError_High->GetN() << endl;

    // For Ratio: first two points in the intermediate range, which have been removed in "Plot_Ratio_ThisAnalysisOnly_Overlapped.
    for (int i=0; i<LowRangeRemoveAtEnd; i++){
        g_LowResult->RemovePoint(g_LowResult->GetN()-1);
    }
    for (int i=0; i<IntermediateRangeRemoveAtBegin-2; i++){
        g_IntermediateResult->RemovePoint(0);
    }
    for (int i=0; i<IntermediateRangeRemoveAtEnd; i++){
        g_IntermediateResult->RemovePoint(g_IntermediateResult->GetN()-1);
    }
    for (int i=0; i<HighRangeRemovedAtBegin; i++){
        g_HighResult->RemovePoint(0);
    }

    // For Error:
    for (int i=0; i<LowRangeRemoveAtEnd; i++){
        g_System_ACC_Low      ->RemovePoint(g_System_ACC_Low->GetN()-1);    
        g_TotalSysError_Low   ->RemovePoint(g_TotalSysError_Low->GetN()-1);
        g_StatisticalError_Low->RemovePoint(g_StatisticalError_Low->GetN()-1);
        g_LowResult_ScaledAcceptanceError->RemovePoint(g_LowResult_ScaledAcceptanceError->GetN()-1);
    }
    for (int i=0; i<IntermediateRangeRemoveAtBegin; i++){
        g_System_ACC_Intermediate    ->RemovePoint(0);
        g_TotalSysError_Intermediate ->RemovePoint(0);
        g_StatisticError_Intermediate->RemovePoint(0);
        g_IntermediateResult_ScaledAcceptanceError->RemovePoint(0);
    }
    for (int i=0; i<IntermediateRangeRemoveAtEnd; i++){
        g_System_ACC_Intermediate    ->RemovePoint(g_System_ACC_Intermediate->GetN()-1);
        g_TotalSysError_Intermediate ->RemovePoint(g_TotalSysError_Intermediate->GetN()-1);
        g_StatisticError_Intermediate->RemovePoint(g_StatisticError_Intermediate->GetN()-1);
        g_IntermediateResult_ScaledAcceptanceError->RemovePoint(g_IntermediateResult_ScaledAcceptanceError->GetN()-1);
    }
    for (int i=0; i<HighRangeRemovedAtBegin; i++){
        g_System_ACC_High      ->RemovePoint(0);
        g_TotalSysError_High   ->RemovePoint(0);
        g_StatisticalError_High->RemovePoint(0);
        g_HighResult_ScaledAcceptanceError->RemovePoint(0);
    }

}


//// Plot


// Calculate Total antiproton and proton Numbers in Aachen result.
void CalculateAntiprotonAndProtonNumbers(TGraph *g_Antiproton_number_unfolded_Low, TGraph *g_Proton_number_unfolded_Low, TGraph *g_Antiproton_number_unfolded_Intermediate, TGraph *g_Proton_number_unfolded_Intermediate, TGraph *g_Antiproton_number_unfolded_High, TGraph *g_Proton_number_unfolded_High, TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber){
    cout<< "\n" <<endl;
    cout<< "Total Number Printing:" <<endl;
    double totalnumber_antiproton = 0;
    double totalnumber_proton = 0;
    double totalnumber_antiproton_High = 0;
    double totalnumber_proton_High = 0;
    double totalnumber_antiproton_Intermediate = 0;
    double totalnumber_proton_Intermediate = 0;
    double totalnumber_antiproton_Low = 0;
    double totalnumber_proton_Low = 0;

    for (int i=0; i<g_Antiproton_number_unfolded_Low->GetN(); i++){
        //cout<< "Low range number X:" << g_Antiproton_number_unfolded_Low->GetX()[i] <<endl;
        //cout<< "Low range number:" << g_Antiproton_number_unfolded_Low->GetY()[i] <<endl;
        totalnumber_antiproton      = totalnumber_antiproton + g_Antiproton_number_unfolded_Low->GetY()[i];
        totalnumber_proton          = totalnumber_proton + g_Proton_number_unfolded_Low->GetY()[i];
        totalnumber_antiproton_Low  = totalnumber_antiproton_Low + g_Antiproton_number_unfolded_Low->GetY()[i];
        totalnumber_proton_Low      = totalnumber_proton_Low + g_Proton_number_unfolded_Low->GetY()[i];
    }
    for (int i=0; i<g_Antiproton_number_unfolded_Intermediate->GetN(); i++){
        //cout<< "Intermediate range number X:" << g_Antiproton_number_unfolded_Intermediate->GetX()[i] <<endl;
        //cout<< "Intermediate range number:" << g_Antiproton_number_unfolded_Intermediate->GetY()[i] <<endl;
        totalnumber_antiproton              = totalnumber_antiproton + g_Antiproton_number_unfolded_Intermediate->GetY()[i];
        totalnumber_proton                  = totalnumber_proton + g_Proton_number_unfolded_Intermediate->GetY()[i];
        totalnumber_antiproton_Intermediate = totalnumber_antiproton_Intermediate + g_Antiproton_number_unfolded_Intermediate->GetY()[i];
        totalnumber_proton_Intermediate     = totalnumber_proton_Intermediate + g_Proton_number_unfolded_Intermediate->GetY()[i];
    }
    for (int i=0; i<g_Antiproton_number_unfolded_High->GetN(); i++){
        //cout<< "High range number X:" << g_Antiproton_number_unfolded_High->GetX()[i] <<endl;
        //cout<< "High range number:" << g_Antiproton_number_unfolded_High->GetY()[i] <<endl;
        totalnumber_antiproton      = totalnumber_antiproton + g_Antiproton_number_unfolded_High->GetY()[i];
        totalnumber_proton          = totalnumber_proton + g_Proton_number_unfolded_High->GetY()[i];
        totalnumber_antiproton_High = totalnumber_antiproton_High + g_Antiproton_number_unfolded_High->GetY()[i];
        totalnumber_proton_High     = totalnumber_proton_High + g_Proton_number_unfolded_High->GetY()[i];
    }

    // g_PhysicsReport_pbarNumber->GetN()=58
    // g_Published_pbarNumber->GetN()=57
    // g_PhysicsReport_pbarNumber->GetX()[27]=15.95
    // g_Published_pbarNumber->GetX()[27]=15.95
    double totalnumber_PhysicsReport = 0;
    double totalnumber_PublishedPRL = 0;
    double totalnumber_PhysicsReport_High = 0;
    double totalnumber_PublishedPRL_High = 0;

    for (int i=0; i<g_PhysicsReport_pbarNumber->GetN(); i++){
        totalnumber_PhysicsReport = totalnumber_PhysicsReport + g_PhysicsReport_pbarNumber->GetY()[i]; 
    }

    for (int i=27; i<g_PhysicsReport_pbarNumber->GetN(); i++){
        totalnumber_PhysicsReport_High = totalnumber_PhysicsReport_High + g_PhysicsReport_pbarNumber->GetY()[i];
    }

    for (int i=0; i<g_Published_pbarNumber->GetN(); i++){
        totalnumber_PublishedPRL = totalnumber_PublishedPRL + g_Published_pbarNumber->GetY()[i];
    }

    for (int i=27; i<g_Published_pbarNumber->GetN(); i++){
        totalnumber_PublishedPRL_High = totalnumber_PublishedPRL_High + g_Published_pbarNumber->GetY()[i];
    }

    cout << "Low Range: "               << g_Antiproton_number_unfolded_Low->GetX()[0] << " to " << g_Antiproton_number_unfolded_Low->GetX()[g_Antiproton_number_unfolded_Low->GetN()] << " GV" <<endl; 
    cout << "IntermediateRange Range: " << g_Antiproton_number_unfolded_Intermediate->GetX()[0] << " to " << g_Antiproton_number_unfolded_Intermediate->GetX()[g_Antiproton_number_unfolded_Intermediate->GetN()] << " GV" <<endl;
    cout << "High Range: "              << g_Antiproton_number_unfolded_High->GetX()[0] << " to " << g_Antiproton_number_unfolded_High->GetX()[g_Antiproton_number_unfolded_High->GetN()]  << " GV" <<endl;

    cout<< "Total Number antiproton (LowRange): "          << totalnumber_antiproton_Low          << endl;
    cout<< "Total Number proton (LowRange): "              << totalnumber_proton_Low              << endl;
    cout<< "Total Number antiproton (IntermediateRange): " << totalnumber_antiproton_Intermediate << endl;
    cout<< "Total Number proton (IntermediateRange): "     << totalnumber_proton_Intermediate     << endl;

    cout<< "Total Number antiproton (PhysicsReport): "     << totalnumber_PhysicsReport           << endl;
    cout<< "Total Number antiproton (PublishedPRL): "      << totalnumber_PublishedPRL            << endl;
    cout<< "Total Number antiproton: "                     << totalnumber_antiproton              << endl;
    cout<< "Total Number proton: "                         << totalnumber_proton                  << endl;

    cout<< "Total Number antiproton (PhysicsReport) (HighRange): " << totalnumber_PhysicsReport_High << endl;
    cout<< "Total Number antiproton (PublishedPRL) (HighRange): "  << totalnumber_PublishedPRL_High  << endl;
    cout<< "Total Number antiproton (This analysis) (HighRange): " << totalnumber_antiproton_High         << endl;
    cout<< "Total Number antiproton (This analysis/PhysicsReport) (HighRange): " << totalnumber_antiproton_High/totalnumber_PhysicsReport_High         << endl;
    cout<< "Total Number antiproton (This analysis/PublishedPRL) (HighRange): "  << totalnumber_antiproton_High/totalnumber_PublishedPRL_High           << endl;

    cout<< "\n" <<endl;
}


// Tool: Deviding by bin width
void DividingBinWidth(TGraph *InputGraph, TH1D *InputHisto){
    for (int i=0; i<InputGraph->GetN(); i++){
        InputGraph->SetPoint(i, InputGraph->GetX()[i], InputGraph->GetY()[i]/InputHisto->GetBinWidth(i+1));
    }
}


// Plot Antiproton and Proton Number.
void Plot_Antiproton_and_Proton_Numbers_OverRigidityBinWidth(TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber, 
TGraph *g_Antiproton_number_unfolded_Low         , TGraph *g_Proton_number_unfolded_Low, 
TGraph *g_Antiproton_number_unfolded_Intermediate, TGraph *g_Proton_number_unfolded_Intermediate, 
TGraph *g_Antiproton_number_unfolded_High        , TGraph *g_Antiproton_number_unfolded_High_P0, TGraph *g_Antiproton_number_unfolded_High_P0VGGNN, TGraph *g_Antiproton_number_unfolded_High_P1, TGraph *g_Antiproton_number_unfolded_High_P2, TGraph *g_Antiproton_number_unfolded_High_P3, TGraph *g_Antiproton_number_unfolded_High_P4, TGraph *g_Antiproton_number_unfolded_High_P5, 
TGraph *g_Proton_number_unfolded_High            , TGraph *g_Proton_number_unfolded_High_P0         , TGraph *g_Proton_number_unfolded_High_P0VGGNN, TGraph *g_Proton_number_unfolded_High_P1   , TGraph *g_Proton_number_unfolded_High_P2    , TGraph *g_Proton_number_unfolded_High_P3    , TGraph *g_Proton_number_unfolded_High_P4    , TGraph *g_Proton_number_unfolded_High_P5, 
TGraph *g_antiproton_number_raw_Low         , TGraph *g_Proton_number_raw_Low, 
TGraph *g_Antiproton_number_raw_Intermediate, TGraph *g_Proton_number_raw_Intermediate, 
TGraph *g_Antiproton_number_raw_High        , TGraph *g_Antiproton_number_raw_High_P0, TGraph *g_Antiproton_number_raw_High_P0VGGNN, TGraph *g_Antiproton_number_raw_High_P1, TGraph *g_Antiproton_number_raw_High_P2, TGraph *g_Antiproton_number_raw_High_P3, TGraph *g_Antiproton_number_raw_High_P4, TGraph *g_Antiproton_number_raw_High_P5, 
TGraph *g_Proton_number_raw_High, TGraph *g_Proton_number_raw_High_P0, TGraph *g_Proton_number_raw_High_P0VGGNN, TGraph *g_Proton_number_raw_High_P1, TGraph *g_Proton_number_raw_High_P2, TGraph *g_Proton_number_raw_High_P3, TGraph *g_Proton_number_raw_High_P4, TGraph *g_Proton_number_raw_High_P5,
TH1D *h_Antiproton_number_unfolded_Low, TH1D *h_Antiproton_number_unfolded_Intermediate, TH1D *h_Antiproton_number_unfolded_High_P0VGGNN, std::string issversion, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin){   // Here P0 useage depends on Input Parameter. No worry about BDT or VGGNN here. 


    //// Dividing by Bin Width
    TGraph *g_Antiproton_number_unfolded_Low_Over_RigidityBin          = new TGraph(*g_Antiproton_number_unfolded_Low);
    TGraph *g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin = new TGraph(*g_Antiproton_number_unfolded_Intermediate);
    TGraph *g_Antiproton_number_unfolded_High_Over_RigidityBin         = new TGraph(*g_Antiproton_number_unfolded_High);
    TGraph *g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin = new TGraph(*g_Antiproton_number_unfolded_High_P0VGGNN);
    TGraph *g_Antiproton_number_unfolded_High_P1_Over_RigidityBin      = new TGraph(*g_Antiproton_number_unfolded_High_P1);
    TGraph *g_Antiproton_number_unfolded_High_P2_Over_RigidityBin      = new TGraph(*g_Antiproton_number_unfolded_High_P2);
    TGraph *g_Antiproton_number_unfolded_High_P4_Over_RigidityBin      = new TGraph(*g_Antiproton_number_unfolded_High_P4);

    TGraph *g_antiproton_number_raw_Low_Over_RigidityBin               = new TGraph(*g_antiproton_number_raw_Low);
    TGraph *g_Antiproton_number_raw_Intermediate_Over_RigidityBin      = new TGraph(*g_Antiproton_number_raw_Intermediate);
    TGraph *g_Antiproton_number_raw_High_Over_RigidityBin              = new TGraph(*g_Antiproton_number_raw_High);
    TGraph *g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin      = new TGraph(*g_Antiproton_number_raw_High_P0VGGNN);
    TGraph *g_Antiproton_number_raw_High_P1_Over_RigidityBin           = new TGraph(*g_Antiproton_number_raw_High_P1);
    TGraph *g_Antiproton_number_raw_High_P2_Over_RigidityBin           = new TGraph(*g_Antiproton_number_raw_High_P2);
    TGraph *g_Antiproton_number_raw_High_P4_Over_RigidityBin           = new TGraph(*g_Antiproton_number_raw_High_P4);

    TGraph *g_Proton_number_unfolded_Low_Over_RigidityBin              = new TGraph(*g_Proton_number_unfolded_Low);
    TGraph *g_Proton_number_unfolded_Intermediate_Over_RigidityBin     = new TGraph(*g_Proton_number_unfolded_Intermediate);

    TGraph *g_Proton_number_raw_Low_Over_RigidityBin                   = new TGraph(*g_Proton_number_raw_Low);
    TGraph *g_Proton_number_raw_Intermediate_Over_RigidityBin          = new TGraph(*g_Proton_number_raw_Intermediate);

    for (int i=0; i<g_Antiproton_number_unfolded_Low->GetN(); i++){
        g_Antiproton_number_unfolded_Low_Over_RigidityBin->SetPoint(i, g_Antiproton_number_unfolded_Low_Over_RigidityBin->GetX()[i], g_Antiproton_number_unfolded_Low_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Low->GetBinWidth(i+1));
        g_antiproton_number_raw_Low_Over_RigidityBin     ->SetPoint(i, g_antiproton_number_raw_Low_Over_RigidityBin->GetX()[i]     , g_antiproton_number_raw_Low_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Low->GetBinWidth(i+1));
        g_Proton_number_unfolded_Low_Over_RigidityBin    ->SetPoint(i, g_Proton_number_unfolded_Low_Over_RigidityBin->GetX()[i]    , g_Proton_number_unfolded_Low_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Low->GetBinWidth(i+1));
        g_Proton_number_raw_Low_Over_RigidityBin         ->SetPoint(i, g_Proton_number_raw_Low_Over_RigidityBin->GetX()[i]         , g_Proton_number_raw_Low_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Low->GetBinWidth(i+1));
    }

    for (int i=0; i<g_Antiproton_number_unfolded_Intermediate->GetN(); i++){
        g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->SetPoint(i, g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->GetX()[i], g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Intermediate->GetBinWidth(i+1+IntermediateRangeRemoveAtEnd));
        g_Antiproton_number_raw_Intermediate_Over_RigidityBin     ->SetPoint(i, g_Antiproton_number_raw_Intermediate_Over_RigidityBin     ->GetX()[i], g_Antiproton_number_raw_Intermediate_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Intermediate->GetBinWidth(i+1+IntermediateRangeRemoveAtEnd));
        g_Proton_number_unfolded_Intermediate_Over_RigidityBin    ->SetPoint(i, g_Proton_number_unfolded_Intermediate_Over_RigidityBin    ->GetX()[i], g_Proton_number_unfolded_Intermediate_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Intermediate->GetBinWidth(i+1+IntermediateRangeRemoveAtEnd));
        g_Proton_number_raw_Intermediate_Over_RigidityBin         ->SetPoint(i, g_Proton_number_raw_Intermediate_Over_RigidityBin         ->GetX()[i], g_Proton_number_raw_Intermediate_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_Intermediate->GetBinWidth(i+1+IntermediateRangeRemoveAtEnd));
    }

    for (int i=0; i<g_Antiproton_number_unfolded_High->GetN(); i++){
        g_Antiproton_number_unfolded_High_Over_RigidityBin->SetPoint(i, g_Antiproton_number_unfolded_High_Over_RigidityBin->GetX()[i], g_Antiproton_number_unfolded_High_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Antiproton_number_raw_High_Over_RigidityBin->SetPoint(i,g_Antiproton_number_raw_High_Over_RigidityBin->GetX()[i], g_Antiproton_number_raw_High_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));

        g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin->SetPoint(i, g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin->GetX()[i], g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Antiproton_number_raw_High_P1_Over_RigidityBin     ->SetPoint(i, g_Antiproton_number_raw_High_P1_Over_RigidityBin->GetX()[i], g_Antiproton_number_raw_High_P1_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Antiproton_number_raw_High_P2_Over_RigidityBin     ->SetPoint(i, g_Antiproton_number_raw_High_P2_Over_RigidityBin->GetX()[i], g_Antiproton_number_raw_High_P2_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Antiproton_number_raw_High_P4_Over_RigidityBin     ->SetPoint(i, g_Antiproton_number_raw_High_P4_Over_RigidityBin->GetX()[i], g_Antiproton_number_raw_High_P4_Over_RigidityBin->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));

    }


    //// index for used range
    int index_38  = 13;
    int index_147 = 27;
    int index_175 = 28;

    //// Calculating Antiproton_number/RigidityBinWidth, and plot
    TGraph *g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used = new TGraph(g_Antiproton_number_unfolded_High_P0VGGNN->GetN());
    TGraph *g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used      = new TGraph(index_147);
    TGraph *g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used      = new TGraph(index_175);
    TGraph *g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used      = new TGraph(index_38);

    TGraph *g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used      = new TGraph(g_Antiproton_number_unfolded_High_P0VGGNN->GetN());
    TGraph *g_Antiproton_number_raw_High_P1_Over_RigidityBin_used           = new TGraph(index_147);
    TGraph *g_Antiproton_number_raw_High_P2_Over_RigidityBin_used           = new TGraph(index_175);
    TGraph *g_Antiproton_number_raw_High_P4_Over_RigidityBin_used           = new TGraph(index_38);

    TGraph *g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used     = new TGraph(g_Proton_number_unfolded_High_P0VGGNN->GetN()); 
    TGraph *g_Proton_number_unfolded_High_P1_Over_RigidityBin_used          = new TGraph(g_Proton_number_unfolded_High_P1->GetN());
    TGraph *g_Proton_number_unfolded_High_P2_Over_RigidityBin_used          = new TGraph(g_Proton_number_unfolded_High_P2->GetN());
    TGraph *g_Proton_number_unfolded_High_P4_Over_RigidityBin_used          = new TGraph(g_Proton_number_unfolded_High_P4->GetN());

    TGraph *g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used          = new TGraph(g_Proton_number_raw_High_P0VGGNN->GetN());
    TGraph *g_Proton_number_raw_High_P1_Over_RigidityBin_used               = new TGraph(g_Proton_number_raw_High_P1->GetN());
    TGraph *g_Proton_number_raw_High_P2_Over_RigidityBin_used               = new TGraph(g_Proton_number_raw_High_P2->GetN());
    TGraph *g_Proton_number_raw_High_P4_Over_RigidityBin_used               = new TGraph(g_Proton_number_raw_High_P4->GetN());

    //Used Area
    for (int i=0; i<g_Antiproton_number_unfolded_High_P0->GetN(); i++){
        g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_unfolded_High_P0VGGNN->GetX()[i], g_Antiproton_number_unfolded_High_P0VGGNN->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used     ->SetPoint(i, g_Antiproton_number_raw_High_P0VGGNN->GetX()[i] , g_Antiproton_number_raw_High_P0VGGNN->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));  
        g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used    ->SetPoint(i, g_Proton_number_unfolded_High_P0VGGNN->GetX()[i], g_Proton_number_unfolded_High_P0VGGNN->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used         ->SetPoint(i, g_Proton_number_raw_High_P0VGGNN->GetX()[i]     , g_Proton_number_raw_High_P0VGGNN->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));   
    } 
    for (int i=0; i<index_147; i++){
        g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_unfolded_High_P1->GetX()[i], g_Antiproton_number_unfolded_High_P1->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin)); 
        g_Antiproton_number_raw_High_P1_Over_RigidityBin_used     ->SetPoint(i, g_Antiproton_number_raw_High_P1->GetX()[i]     , g_Antiproton_number_raw_High_P1->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));  
        g_Proton_number_unfolded_High_P1_Over_RigidityBin_used    ->SetPoint(i, g_Proton_number_unfolded_High_P1->GetX()[i]    , g_Proton_number_unfolded_High_P1->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin)); 
        g_Proton_number_raw_High_P1_Over_RigidityBin_used         ->SetPoint(i, g_Proton_number_raw_High_P1->GetX()[i]         , g_Proton_number_raw_High_P1->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin)); 
    }
    for (int i=0; i<index_175; i++){
        g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_unfolded_High_P2->GetX()[i], g_Antiproton_number_unfolded_High_P2->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Antiproton_number_raw_High_P2_Over_RigidityBin_used     ->SetPoint(i, g_Antiproton_number_raw_High_P2->GetX()[i]     , g_Antiproton_number_raw_High_P2->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_unfolded_High_P2_Over_RigidityBin_used    ->SetPoint(i, g_Proton_number_unfolded_High_P2->GetX()[i]    , g_Proton_number_unfolded_High_P2->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_raw_High_P2_Over_RigidityBin_used         ->SetPoint(i, g_Proton_number_raw_High_P2->GetX()[i]         , g_Proton_number_raw_High_P2->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin)); 
    }
    for (int i=0; i<index_38; i++){
        g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_unfolded_High_P4->GetX()[i], g_Antiproton_number_unfolded_High_P4->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Antiproton_number_raw_High_P4_Over_RigidityBin_used     ->SetPoint(i, g_Antiproton_number_raw_High_P4->GetX()[i]     , g_Antiproton_number_raw_High_P4->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_unfolded_High_P4_Over_RigidityBin_used    ->SetPoint(i, g_Proton_number_unfolded_High_P4->GetX()[i]    , g_Proton_number_unfolded_High_P4->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin)); 
        g_Proton_number_raw_High_P4_Over_RigidityBin_used         ->SetPoint(i, g_Proton_number_raw_High_P4->GetX()[i]         , g_Proton_number_raw_High_P4->GetY()[i]/h_Antiproton_number_unfolded_High_P0VGGNN->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }


    //// Plot: (1). Antiproton numbers
    TCanvas c_number("c_number","c_number",1000, 500);

    TPad *pad_number   = new TPad("", "", 0.0, 0.28, 1.0, 1.0 , 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    TPad *pad_residual = new TPad("", "", 0.0, 0.0 , 1.0, 0.355, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)

    pad_number->Draw();
    pad_residual->Draw();

    // upper plot: number plot
    pad_number->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    g_Antiproton_number_unfolded_Low_Over_RigidityBin         ->Draw("AP");
    g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->Draw("same P");
    //g_Antiproton_number_unfolded_High_Over_RigidityBin        ->Draw("same P");

    //g_PhysicsReport_pbarNumber->Draw("AP");
    //g_Published_pbarNumber->Draw("same P");

    g_antiproton_number_raw_Low_Over_RigidityBin         ->Draw("same P");               
    g_Antiproton_number_raw_Intermediate_Over_RigidityBin->Draw("same P");
    //g_Antiproton_number_raw_High_Over_RigidityBin        ->Draw("same P");

    /*
    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin ->Draw("same P");
    g_Antiproton_number_raw_High_P1_Over_RigidityBin      ->Draw("same P");
    g_Antiproton_number_raw_High_P2_Over_RigidityBin      ->Draw("same P");
    g_Antiproton_number_raw_High_P4_Over_RigidityBin      ->Draw("same P");
    */

    g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->Draw("same P");
    g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used      ->Draw("same P");
    g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used      ->Draw("same P");
    g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used      ->Draw("same P");

    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->Draw("same P");
    g_Antiproton_number_raw_High_P1_Over_RigidityBin_used      ->Draw("same P");
    g_Antiproton_number_raw_High_P2_Over_RigidityBin_used      ->Draw("same P");
    g_Antiproton_number_raw_High_P4_Over_RigidityBin_used      ->Draw("same P");

    TAxis * xaxis_number = g_Antiproton_number_unfolded_Low_Over_RigidityBin->GetXaxis();
    TAxis * yaxis_number = g_Antiproton_number_unfolded_Low_Over_RigidityBin->GetYaxis();
    xaxis_number->SetLimits(1, 600); //1, 600
    yaxis_number->SetRangeUser(0.1, 90000); //1000, 90000 
    //yaxis_number->SetMoreLogLabels();
    xaxis_number->SetTitle("|R| / (GV)");
    yaxis_number->SetTitle("N_{P} / #Delta R");

    TLine *line_Pbar = new TLine(0.02, 0.76, 0.02, 0.78);
    line_Pbar->SetLineColor(kBlack);
    line_Pbar->SetLineWidth(1);
    line_Pbar->SetNDC(1);
    line_Pbar->Draw();

    xaxis_number->SetTitleFont(62);
    yaxis_number->SetTitleFont(62);
    xaxis_number->SetTitleSize(0.05);
    yaxis_number->SetTitleSize(0.05);
    xaxis_number->SetLabelFont(62);
    xaxis_number->SetLabelSize(0.05);
    yaxis_number->SetLabelFont(62);
    yaxis_number->SetLabelSize(0.05);

    xaxis_number->SetTitleOffset(1.3);
    //yaxis_number->SetTitleOffset(1.5);
    gPad->SetBottomMargin(0.13);
    //gPad->SetLeftMargin(0.16);

    xaxis_number->SetLabelOffset(999); // remove x axis
    xaxis_number->SetLabelSize(0);     // remove x label

    g_Published_pbarNumber->SetMarkerStyle(15);
    g_Published_pbarNumber->SetMarkerColor(0); //1
    g_PhysicsReport_pbarNumber->SetMarkerStyle(15);
    g_PhysicsReport_pbarNumber->SetMarkerColor(0);//3
    g_PhysicsReport_pbarNumber->SetMarkerSize(0);

    g_antiproton_number_raw_Low_Over_RigidityBin         ->SetMarkerStyle(15);
    g_antiproton_number_raw_Low_Over_RigidityBin         ->SetMarkerColor(6);
    g_antiproton_number_raw_Low_Over_RigidityBin         ->SetMarkerSize(0.9);
    g_Antiproton_number_raw_Intermediate_Over_RigidityBin->SetMarkerStyle(15);
    g_Antiproton_number_raw_Intermediate_Over_RigidityBin->SetMarkerColor(6);
    g_Antiproton_number_raw_Intermediate_Over_RigidityBin->SetMarkerSize(0.9);
    g_Antiproton_number_raw_High_Over_RigidityBin        ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_Over_RigidityBin        ->SetMarkerColor(6); 
    g_Antiproton_number_raw_High_Over_RigidityBin        ->SetMarkerSize(0.9);

    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin ->SetMarkerColor(46);
    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin ->SetMarkerSize(0.9);
    g_Antiproton_number_raw_High_P1_Over_RigidityBin      ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P1_Over_RigidityBin      ->SetMarkerColor(41);
    g_Antiproton_number_raw_High_P1_Over_RigidityBin      ->SetMarkerSize(0.9);
    g_Antiproton_number_raw_High_P2_Over_RigidityBin      ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P2_Over_RigidityBin      ->SetMarkerColor(38);
    g_Antiproton_number_raw_High_P2_Over_RigidityBin      ->SetMarkerSize(0.9);
    g_Antiproton_number_raw_High_P4_Over_RigidityBin      ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P4_Over_RigidityBin      ->SetMarkerColor(49);
    g_Antiproton_number_raw_High_P4_Over_RigidityBin      ->SetMarkerSize(0.9);

    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerColor(46);
    g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerSize(0.9);
    g_Antiproton_number_raw_High_P1_Over_RigidityBin_used      ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P1_Over_RigidityBin_used      ->SetMarkerColor(41);
    g_Antiproton_number_raw_High_P1_Over_RigidityBin_used      ->SetMarkerSize(0.9);
    g_Antiproton_number_raw_High_P2_Over_RigidityBin_used      ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P2_Over_RigidityBin_used      ->SetMarkerColor(38);
    g_Antiproton_number_raw_High_P2_Over_RigidityBin_used      ->SetMarkerSize(0.9);
    g_Antiproton_number_raw_High_P4_Over_RigidityBin_used      ->SetMarkerStyle(15);
    g_Antiproton_number_raw_High_P4_Over_RigidityBin_used      ->SetMarkerColor(49);
    g_Antiproton_number_raw_High_P4_Over_RigidityBin_used      ->SetMarkerSize(0.9);

    g_Antiproton_number_unfolded_Low_Over_RigidityBin->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_Low_Over_RigidityBin->SetMarkerColor(4);
    g_Antiproton_number_unfolded_Low_Over_RigidityBin->SetMarkerSize(0.9);
    g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->SetMarkerColor(4);
    g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->SetMarkerSize(0.9);
    g_Antiproton_number_unfolded_High_Over_RigidityBin->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_Over_RigidityBin->SetMarkerColor(4);
    g_Antiproton_number_unfolded_High_Over_RigidityBin->SetMarkerSize(0.9);

    g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerColor(8);
    g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerSize(0.9);
    g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used ->SetMarkerColor(9);
    g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used ->SetMarkerSize(0.9);
    g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used ->SetMarkerColor(12);
    g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used ->SetMarkerSize(0.9);
    g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used ->SetMarkerColor(30);
    g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used ->SetMarkerSize(0.9);

    g_Antiproton_number_unfolded_Low_Over_RigidityBin->SetTitle("");
    gPad->SetLogy();
    gPad->SetLogx();

    TLegend *legend_number = new TLegend(0.18, 0.15, 0.68, 0.42);
    legend_number->AddEntry(g_antiproton_number_raw_Low_Over_RigidityBin              , "InnerCentral (Raw)"      , "p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_Low_Over_RigidityBin         , "InnerCentral (Unfolded)" , "p");
    legend_number->AddEntry(g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used, "FullSpan (Raw)"          , "p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used, "FullSpan (Unfolded)", "p");
    legend_number->AddEntry(g_Antiproton_number_raw_High_P1_Over_RigidityBin_used     , "Inner + L1 (Raw)"        , "p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used, "Inner + L1 (Unfolded)"   , "p");
    legend_number->AddEntry(g_Antiproton_number_raw_High_P2_Over_RigidityBin_used     , "Inner + L9 (Raw)"        , "p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used, "Inner + L9 (Unfolded)"   , "p");
    legend_number->AddEntry(g_Antiproton_number_raw_High_P4_Over_RigidityBin_used     , "Inner Only (Raw)"        , "p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used, "Inner Only (Unfolded)"   , "p");

    //legend_number->AddEntry(g_Published_pbarNumber                   , "PRL paper 2016","p");
    //legend_number->AddEntry(g_PhysicsReport_pbarNumber               , "PhysicsReport 2021", "p");
    //legend_number->AddEntry(g_Antiproton_number_unfolded_Low         , "This analysis (LowRange)","p");
    //legend_number->AddEntry(g_Antiproton_number_unfolded_Intermediate, "This analysis (IntermediateRange)","p");
    //legend_number->AddEntry(g_Antiproton_number_unfolded_High        , "This analysis","p");
    legend_number->SetTextSize(0.04);
    legend_number->SetTextFont(62);
    legend_number->SetBorderSize(0);
    legend_number->SetNColumns(2);
    legend_number->Draw();

    CalculateAntiprotonAndProtonNumbers(g_Antiproton_number_unfolded_Low, g_Proton_number_unfolded_Low, g_Antiproton_number_unfolded_Intermediate, g_Proton_number_unfolded_Intermediate, g_Antiproton_number_unfolded_High, g_Proton_number_unfolded_High, g_PhysicsReport_pbarNumber, g_Published_pbarNumber);

    // lower plot: Number Ratio Plot
    pad_residual->cd();

    TGraph *g_PbarRawUnfoldedRatio_Low           = new TGraph();
    TGraph *g_PbarRawUnfoldedRatio_Intermediate  = new TGraph();
    TGraph *g_PbarRawUnfoldedRatio_P0VGGNN       = new TGraph();
    TGraph *g_PbarRawUnfoldedRatio_P1            = new TGraph();
    TGraph *g_PbarRawUnfoldedRatio_P2            = new TGraph();
    TGraph *g_PbarRawUnfoldedRatio_P4            = new TGraph();

    tie(g_PbarRawUnfoldedRatio_Low)          = GraphDivision(g_Antiproton_number_unfolded_Low_Over_RigidityBin         , g_antiproton_number_raw_Low_Over_RigidityBin);
    tie(g_PbarRawUnfoldedRatio_Intermediate) = GraphDivision(g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin, g_Antiproton_number_raw_Intermediate_Over_RigidityBin);
    tie(g_PbarRawUnfoldedRatio_P0VGGNN)      = GraphDivision(g_Antiproton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used, g_Antiproton_number_raw_High_P0VGGNN_Over_RigidityBin_used);
    tie(g_PbarRawUnfoldedRatio_P1)           = GraphDivision(g_Antiproton_number_unfolded_High_P1_Over_RigidityBin_used, g_Antiproton_number_raw_High_P1_Over_RigidityBin_used);
    tie(g_PbarRawUnfoldedRatio_P2)           = GraphDivision(g_Antiproton_number_unfolded_High_P2_Over_RigidityBin_used, g_Antiproton_number_raw_High_P2_Over_RigidityBin_used);
    tie(g_PbarRawUnfoldedRatio_P4)           = GraphDivision(g_Antiproton_number_unfolded_High_P4_Over_RigidityBin_used, g_Antiproton_number_raw_High_P4_Over_RigidityBin_used);
     
    g_PbarRawUnfoldedRatio_Low         ->Draw("AP");
    g_PbarRawUnfoldedRatio_Intermediate->Draw("same P");
    g_PbarRawUnfoldedRatio_P0VGGNN     ->Draw("same P");
    g_PbarRawUnfoldedRatio_P1          ->Draw("same P");
    g_PbarRawUnfoldedRatio_P2          ->Draw("same P");
    g_PbarRawUnfoldedRatio_P4          ->Draw("same P");

    g_PbarRawUnfoldedRatio_Low->SetMarkerStyle(15);
    g_PbarRawUnfoldedRatio_Low->SetMarkerColor(4);
    g_PbarRawUnfoldedRatio_Low->SetMarkerSize(0.9);

    g_PbarRawUnfoldedRatio_Intermediate->SetMarkerStyle(15);
    g_PbarRawUnfoldedRatio_Intermediate->SetMarkerColor(4);
    g_PbarRawUnfoldedRatio_Intermediate->SetMarkerSize(0.9);

    g_PbarRawUnfoldedRatio_P0VGGNN->SetMarkerStyle(15);
    g_PbarRawUnfoldedRatio_P0VGGNN->SetMarkerColor(8);
    g_PbarRawUnfoldedRatio_P0VGGNN->SetMarkerSize(0.9);
    g_PbarRawUnfoldedRatio_P1->SetMarkerStyle(15);
    g_PbarRawUnfoldedRatio_P1->SetMarkerColor(9);
    g_PbarRawUnfoldedRatio_P1->SetMarkerSize(0.9);
    g_PbarRawUnfoldedRatio_P2->SetMarkerStyle(15);
    g_PbarRawUnfoldedRatio_P2->SetMarkerColor(12);
    g_PbarRawUnfoldedRatio_P2->SetMarkerSize(0.9);
    g_PbarRawUnfoldedRatio_P4->SetMarkerStyle(15);
    g_PbarRawUnfoldedRatio_P4->SetMarkerColor(30);
    g_PbarRawUnfoldedRatio_P4->SetMarkerSize(0.9);


    TAxis * xaxis_RawUnfoldedPbarCountRatio = g_PbarRawUnfoldedRatio_Low->GetXaxis();
    TAxis * yaxis_RawUnfoldedPbarCountRatio = g_PbarRawUnfoldedRatio_Low->GetYaxis();
    xaxis_RawUnfoldedPbarCountRatio->SetLimits(1, 600); //1, 600
    yaxis_RawUnfoldedPbarCountRatio->SetRangeUser(0, 2.5);
    xaxis_RawUnfoldedPbarCountRatio->SetTitle("|R| / (GV)");
    yaxis_RawUnfoldedPbarCountRatio->SetTitle("Unfolded / Raw");
    xaxis_RawUnfoldedPbarCountRatio->SetTitleFont(62);
    yaxis_RawUnfoldedPbarCountRatio->SetTitleFont(62);
    xaxis_RawUnfoldedPbarCountRatio->SetTitleSize(0.1);
    yaxis_RawUnfoldedPbarCountRatio->SetTitleSize(0.09);
    xaxis_RawUnfoldedPbarCountRatio->SetLabelFont(62);
    xaxis_RawUnfoldedPbarCountRatio->SetLabelSize(0.1);
    yaxis_RawUnfoldedPbarCountRatio->SetLabelFont(62);
    yaxis_RawUnfoldedPbarCountRatio->SetLabelSize(0.1);
    xaxis_RawUnfoldedPbarCountRatio->SetMoreLogLabels();

    xaxis_RawUnfoldedPbarCountRatio->SetTitleOffset(1.2);
    yaxis_RawUnfoldedPbarCountRatio->SetTitleOffset(0.45);
    pad_residual->SetBottomMargin(0.3);
    //gPad->SetLeftMargin(0.16);

    gPad->SetLogx();

    TLegend *legend_RawUnfoldedPbarCountRatio = new TLegend(0.18, 0.70, 0.68, 0.85);
    legend_RawUnfoldedPbarCountRatio->AddEntry(g_PbarRawUnfoldedRatio_Low         , "InnerCentral"  , "p");
    //legend_RawUnfoldedPbarCountRatio->AddEntry(g_PbarRawUnfoldedRatio_Intermediate, "InnerCentral", "p");
    legend_RawUnfoldedPbarCountRatio->AddEntry(g_PbarRawUnfoldedRatio_P0VGGNN     , "FullSpan"      , "p");
    legend_RawUnfoldedPbarCountRatio->AddEntry(g_PbarRawUnfoldedRatio_P1          , "Inner + L1"    , "p");
    legend_RawUnfoldedPbarCountRatio->AddEntry(g_PbarRawUnfoldedRatio_P2          , "Inner + L9"    , "p");
    legend_RawUnfoldedPbarCountRatio->AddEntry(g_PbarRawUnfoldedRatio_P4          , "Inner Only"    , "p");

    legend_RawUnfoldedPbarCountRatio->SetTextSize(0.09);
    legend_RawUnfoldedPbarCountRatio->SetTextFont(62);
    legend_RawUnfoldedPbarCountRatio->SetBorderSize(0);
    legend_RawUnfoldedPbarCountRatio->SetNColumns(3);
    legend_RawUnfoldedPbarCountRatio->Draw();

    c_number.SaveAs( (string("NumberPlot_") + issversion + string(".pdf")).c_str());

    //// Plot: (2). Proton numbers
    TCanvas c_p_number("","",1000,500);

    TPad *pad_p_number   = new TPad("", "", 0.0, 0.28, 1.0, 1.0 , 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    TPad *pad_p_residual = new TPad("", "", 0.0, 0.0 , 1.0, 0.355, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)

    pad_p_number->Draw();
    pad_p_residual->Draw();

    // upper plot: number plot
    pad_p_number->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    g_Proton_number_unfolded_Low_Over_RigidityBin->Draw("AP");
    g_Proton_number_unfolded_Intermediate_Over_RigidityBin->Draw("same P");

    g_Proton_number_raw_Low_Over_RigidityBin->Draw("same P"); 
    g_Proton_number_raw_Intermediate_Over_RigidityBin->Draw("same P");

    g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->Draw("same P");
    g_Proton_number_unfolded_High_P1_Over_RigidityBin_used      ->Draw("same P");
    g_Proton_number_unfolded_High_P2_Over_RigidityBin_used      ->Draw("same P");
    g_Proton_number_unfolded_High_P4_Over_RigidityBin_used      ->Draw("same P");

    g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->Draw("same P");
    g_Proton_number_raw_High_P1_Over_RigidityBin_used      ->Draw("same P");
    g_Proton_number_raw_High_P2_Over_RigidityBin_used      ->Draw("same P");
    g_Proton_number_raw_High_P4_Over_RigidityBin_used      ->Draw("same P");
    
    g_Proton_number_raw_Low_Over_RigidityBin         ->SetMarkerStyle(15);
    g_Proton_number_raw_Low_Over_RigidityBin         ->SetMarkerColor(6);
    g_Proton_number_raw_Low_Over_RigidityBin         ->SetMarkerSize(0.9);
    g_Proton_number_raw_Intermediate_Over_RigidityBin->SetMarkerStyle(15);
    g_Proton_number_raw_Intermediate_Over_RigidityBin->SetMarkerColor(6);
    g_Proton_number_raw_Intermediate_Over_RigidityBin->SetMarkerSize(0.9);

    g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerColor(46);
    g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerSize(0.9);
    g_Proton_number_raw_High_P1_Over_RigidityBin_used      ->SetMarkerStyle(15);
    g_Proton_number_raw_High_P1_Over_RigidityBin_used      ->SetMarkerColor(41);
    g_Proton_number_raw_High_P1_Over_RigidityBin_used      ->SetMarkerSize(0.9);
    g_Proton_number_raw_High_P2_Over_RigidityBin_used      ->SetMarkerStyle(15);
    g_Proton_number_raw_High_P2_Over_RigidityBin_used      ->SetMarkerColor(38);
    g_Proton_number_raw_High_P2_Over_RigidityBin_used      ->SetMarkerSize(0.9);
    g_Proton_number_raw_High_P4_Over_RigidityBin_used      ->SetMarkerStyle(15);
    g_Proton_number_raw_High_P4_Over_RigidityBin_used      ->SetMarkerColor(49);
    g_Proton_number_raw_High_P4_Over_RigidityBin_used      ->SetMarkerSize(0.9);

    g_Proton_number_unfolded_Low_Over_RigidityBin->SetMarkerStyle(15);
    g_Proton_number_unfolded_Low_Over_RigidityBin->SetMarkerColor(4);
    g_Proton_number_unfolded_Low_Over_RigidityBin->SetMarkerSize(0.9);
    g_Proton_number_unfolded_Intermediate_Over_RigidityBin->SetMarkerStyle(15);
    g_Proton_number_unfolded_Intermediate_Over_RigidityBin->SetMarkerColor(4);
    g_Proton_number_unfolded_Intermediate_Over_RigidityBin->SetMarkerSize(0.9);

    g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerColor(8);
    g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used ->SetMarkerSize(0.9);
    g_Proton_number_unfolded_High_P1_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Proton_number_unfolded_High_P1_Over_RigidityBin_used ->SetMarkerColor(9);
    g_Proton_number_unfolded_High_P1_Over_RigidityBin_used ->SetMarkerSize(0.9);
    g_Proton_number_unfolded_High_P2_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Proton_number_unfolded_High_P2_Over_RigidityBin_used ->SetMarkerColor(12);
    g_Proton_number_unfolded_High_P2_Over_RigidityBin_used ->SetMarkerSize(0.9); 
    g_Proton_number_unfolded_High_P4_Over_RigidityBin_used ->SetMarkerStyle(15);
    g_Proton_number_unfolded_High_P4_Over_RigidityBin_used ->SetMarkerColor(30);
    g_Proton_number_unfolded_High_P4_Over_RigidityBin_used ->SetMarkerSize(0.9);

    TAxis * xaxis_p = g_Proton_number_unfolded_Low_Over_RigidityBin->GetXaxis();
    TAxis * yaxis_p = g_Proton_number_unfolded_Low_Over_RigidityBin->GetYaxis();
    xaxis_p->SetLimits(1, 600); //1, 600
    yaxis_p->SetRangeUser(100, 10000000000); //1000, 90000
    xaxis_p->SetTitle("|R| / (GV)");
    yaxis_p->SetTitle("N_{P} / #Delta R");
    xaxis_p->SetTitleFont(62);
    yaxis_p->SetTitleFont(62);
    xaxis_p->SetTitleSize(0.05);
    yaxis_p->SetTitleSize(0.05);
    xaxis_p->SetLabelFont(62);
    xaxis_p->SetLabelSize(0.05);
    yaxis_p->SetLabelFont(62);
    yaxis_p->SetLabelSize(0.05);

    xaxis_p->SetTitleOffset(1.3);
    //yaxis_p->SetTitleOffset(1.5);
    gPad->SetBottomMargin(0.13);
    //gPad->SetLeftMargin(0.16);

    xaxis_p->SetLabelOffset(999); // remove x axis
    xaxis_p->SetLabelSize(0);     // remove x label

    g_Proton_number_unfolded_Low_Over_RigidityBin->SetTitle("");
    gPad->SetLogy();
    gPad->SetLogx();
 
    TLegend *legend_number_p = new TLegend(0.18, 0.15, 0.68, 0.42);
    legend_number_p->AddEntry(g_Proton_number_raw_Low_Over_RigidityBin                   , "InnerCentral (Raw)"     , "p");
    legend_number_p->AddEntry(g_Proton_number_unfolded_Low_Over_RigidityBin              , "InnerCentral (Unfolded)", "p");
    legend_number_p->AddEntry(g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used     , "FullSpan (Raw)"         , "p");
    legend_number_p->AddEntry(g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used, "FullSpan (Unfolded)"    , "p");
    legend_number_p->AddEntry(g_Proton_number_raw_High_P1_Over_RigidityBin_used          , "Inner + L1 (Raw)"       , "p");
    legend_number_p->AddEntry(g_Proton_number_unfolded_High_P1_Over_RigidityBin_used     , "Inner + L1 (Unfolded)"  , "p");
    legend_number_p->AddEntry(g_Proton_number_raw_High_P2_Over_RigidityBin_used          , "Inner + L9 (Raw)"       , "p");
    legend_number_p->AddEntry(g_Proton_number_unfolded_High_P2_Over_RigidityBin_used     , "Inner + L9 (Unfolded)"  , "p");
    legend_number_p->AddEntry(g_Proton_number_raw_High_P4_Over_RigidityBin_used          , "Inner Only (Raw)"       , "p");
    legend_number_p->AddEntry(g_Proton_number_unfolded_High_P4_Over_RigidityBin_used     , "Inner Only (Unfolded)"  , "p");

    legend_number_p->SetTextSize(0.04);
    legend_number_p->SetTextFont(62);
    legend_number_p->SetBorderSize(0);
    legend_number_p->SetNColumns(2);
    legend_number_p->Draw();


    // lower plot: Number Ratio Plot 
    pad_p_residual->cd();

    TGraph *g_ProtonRawUnfoldedRatio_Low           = new TGraph();
    TGraph *g_ProtonRawUnfoldedRatio_Intermediate  = new TGraph();
    TGraph *g_ProtonRawUnfoldedRatio_P0VGGNN       = new TGraph();
    TGraph *g_ProtonRawUnfoldedRatio_P1            = new TGraph();
    TGraph *g_ProtonRawUnfoldedRatio_P2            = new TGraph();
    TGraph *g_ProtonRawUnfoldedRatio_P4            = new TGraph();

    tie(g_ProtonRawUnfoldedRatio_Low)          = GraphDivision(g_Proton_number_unfolded_Low_Over_RigidityBin         , g_Proton_number_raw_Low_Over_RigidityBin);
    tie(g_ProtonRawUnfoldedRatio_Intermediate) = GraphDivision(g_Proton_number_unfolded_Intermediate_Over_RigidityBin, g_Proton_number_raw_Intermediate_Over_RigidityBin);
    tie(g_ProtonRawUnfoldedRatio_P0VGGNN)      = GraphDivision(g_Proton_number_unfolded_High_P0VGGNN_Over_RigidityBin_used, g_Proton_number_raw_High_P0VGGNN_Over_RigidityBin_used);
    tie(g_ProtonRawUnfoldedRatio_P1)           = GraphDivision(g_Proton_number_unfolded_High_P1_Over_RigidityBin_used, g_Proton_number_raw_High_P1_Over_RigidityBin_used);
    tie(g_ProtonRawUnfoldedRatio_P2)           = GraphDivision(g_Proton_number_unfolded_High_P2_Over_RigidityBin_used, g_Proton_number_raw_High_P2_Over_RigidityBin_used);
    tie(g_ProtonRawUnfoldedRatio_P4)           = GraphDivision(g_Proton_number_unfolded_High_P4_Over_RigidityBin_used, g_Proton_number_raw_High_P4_Over_RigidityBin_used);

    g_ProtonRawUnfoldedRatio_Low         ->Draw("AP");
    g_ProtonRawUnfoldedRatio_Intermediate->Draw("same P");
    g_ProtonRawUnfoldedRatio_P0VGGNN     ->Draw("same P");
    g_ProtonRawUnfoldedRatio_P1          ->Draw("same P");
    g_ProtonRawUnfoldedRatio_P2          ->Draw("same P");
    g_ProtonRawUnfoldedRatio_P4          ->Draw("same P");

    g_ProtonRawUnfoldedRatio_Low->SetMarkerStyle(15);
    g_ProtonRawUnfoldedRatio_Low->SetMarkerColor(4);
    g_ProtonRawUnfoldedRatio_Low->SetMarkerSize(0.9);

    g_ProtonRawUnfoldedRatio_Intermediate->SetMarkerStyle(15);
    g_ProtonRawUnfoldedRatio_Intermediate->SetMarkerColor(4);
    g_ProtonRawUnfoldedRatio_Intermediate->SetMarkerSize(0.9);

    g_ProtonRawUnfoldedRatio_P0VGGNN->SetMarkerStyle(15);
    g_ProtonRawUnfoldedRatio_P0VGGNN->SetMarkerColor(8);
    g_ProtonRawUnfoldedRatio_P0VGGNN->SetMarkerSize(0.9);
    g_ProtonRawUnfoldedRatio_P1->SetMarkerStyle(15);
    g_ProtonRawUnfoldedRatio_P1->SetMarkerColor(9);
    g_ProtonRawUnfoldedRatio_P1->SetMarkerSize(0.9);
    g_ProtonRawUnfoldedRatio_P2->SetMarkerStyle(15);
    g_ProtonRawUnfoldedRatio_P2->SetMarkerColor(12);
    g_ProtonRawUnfoldedRatio_P2->SetMarkerSize(0.9);
    g_ProtonRawUnfoldedRatio_P4->SetMarkerStyle(15);
    g_ProtonRawUnfoldedRatio_P4->SetMarkerColor(30);
    g_ProtonRawUnfoldedRatio_P4->SetMarkerSize(0.9);

    TAxis * xaxis_RawUnfoldedProtonCountRatio = g_ProtonRawUnfoldedRatio_Low->GetXaxis();
    TAxis * yaxis_RawUnfoldedProtonCountRatio = g_ProtonRawUnfoldedRatio_Low->GetYaxis();
    xaxis_RawUnfoldedProtonCountRatio->SetLimits(1, 600); //1, 600
    yaxis_RawUnfoldedProtonCountRatio->SetRangeUser(0, 2.5);
    xaxis_RawUnfoldedProtonCountRatio->SetTitle("|R| / (GV)");
    yaxis_RawUnfoldedProtonCountRatio->SetTitle("Unfolded / Raw");
    xaxis_RawUnfoldedProtonCountRatio->SetTitleFont(62);
    yaxis_RawUnfoldedProtonCountRatio->SetTitleFont(62);
    xaxis_RawUnfoldedProtonCountRatio->SetTitleSize(0.1);
    yaxis_RawUnfoldedProtonCountRatio->SetTitleSize(0.09);
    xaxis_RawUnfoldedProtonCountRatio->SetLabelFont(62);
    xaxis_RawUnfoldedProtonCountRatio->SetLabelSize(0.1);
    yaxis_RawUnfoldedProtonCountRatio->SetLabelFont(62);
    yaxis_RawUnfoldedProtonCountRatio->SetLabelSize(0.1);
    xaxis_RawUnfoldedProtonCountRatio->SetMoreLogLabels();
    
    xaxis_RawUnfoldedProtonCountRatio->SetTitleOffset(1.2);
    yaxis_RawUnfoldedProtonCountRatio->SetTitleOffset(0.45);
    gPad->SetBottomMargin(0.3);
    //gPad->SetLeftMargin(0.16);

    gPad->SetLogx();

    TLegend *legend_RawUnfoldedProtonCountRatio = new TLegend(0.18, 0.70, 0.68, 0.85);
    legend_RawUnfoldedProtonCountRatio->AddEntry(g_ProtonRawUnfoldedRatio_Low         , "InnerCentral"  , "p");
    //legend_RawUnfoldedProtonCountRatio->AddEntry(g_ProtonRawUnfoldedRatio_Intermediate, "InnerCentral", "p");
    legend_RawUnfoldedProtonCountRatio->AddEntry(g_ProtonRawUnfoldedRatio_P0VGGNN     , "FullSpan"      , "p");
    legend_RawUnfoldedProtonCountRatio->AddEntry(g_ProtonRawUnfoldedRatio_P1          , "Inner + L1"    , "p");
    legend_RawUnfoldedProtonCountRatio->AddEntry(g_ProtonRawUnfoldedRatio_P2          , "Inner + L9"    , "p");
    legend_RawUnfoldedProtonCountRatio->AddEntry(g_ProtonRawUnfoldedRatio_P4          , "Inner Only"    , "p");

    legend_RawUnfoldedProtonCountRatio->SetTextSize(0.09);
    legend_RawUnfoldedProtonCountRatio->SetTextFont(62);
    legend_RawUnfoldedProtonCountRatio->SetBorderSize(0);
    legend_RawUnfoldedProtonCountRatio->SetNColumns(3);
    legend_RawUnfoldedProtonCountRatio->Draw();    

    c_p_number.SaveAs( (string("ProtonNumberPlot_") + issversion + string(".pdf")).c_str());


    //// Plot: (3). Raw/Unfold Number Ratio 
    TCanvas c_number_ratio("c_number_ratio","c_number_ratio",1000,500);    

    TGraph *g_LowAntiprotonNumberRawUnfoldedRatio          = new TGraph(g_Antiproton_number_unfolded_Low_Over_RigidityBin->GetN());
    TGraph *g_IntermediateAntiprotonNumberRawUnfoldedRatio = new TGraph(g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->GetN());
    TGraph *g_HighAntiprotonNumberRawUnfoldedRatio         = new TGraph(g_Antiproton_number_unfolded_High_Over_RigidityBin->GetN());

    TGraph *g_HighAntiprotonNumberRawUnfoldedRatio_P0VGGNN = new TGraph(g_Antiproton_number_raw_High_P0VGGNN->GetN());
    TGraph *g_HighAntiprotonNumberRawUnfoldedRatio_P0      = new TGraph(g_Antiproton_number_raw_High_P0->GetN());
    TGraph *g_HighAntiprotonNumberRawUnfoldedRatio_P1      = new TGraph(g_Antiproton_number_raw_High_P1->GetN());
    TGraph *g_HighAntiprotonNumberRawUnfoldedRatio_P2      = new TGraph(g_Antiproton_number_raw_High_P2->GetN());
    TGraph *g_HighAntiprotonNumberRawUnfoldedRatio_P4      = new TGraph(g_Antiproton_number_raw_High_P4->GetN());

    for (int i=0; i<g_Antiproton_number_unfolded_Low_Over_RigidityBin->GetN(); i++){
        g_LowAntiprotonNumberRawUnfoldedRatio->SetPoint(i, g_Antiproton_number_unfolded_Low->GetX()[i], g_Antiproton_number_unfolded_Low->GetY()[i]/g_antiproton_number_raw_Low->GetY()[i]);
    }
    for (int i=0; i<g_Antiproton_number_unfolded_Intermediate_Over_RigidityBin->GetN(); i++){
        g_IntermediateAntiprotonNumberRawUnfoldedRatio->SetPoint(i, g_Antiproton_number_unfolded_Intermediate->GetX()[i], g_Antiproton_number_unfolded_Intermediate->GetY()[i]/g_Antiproton_number_raw_Intermediate->GetY()[i]);  
    }
    for (int i=0; i<g_Antiproton_number_unfolded_High_Over_RigidityBin->GetN(); i++){ 
        g_HighAntiprotonNumberRawUnfoldedRatio->SetPoint(i, g_Antiproton_number_unfolded_High->GetX()[i], g_Antiproton_number_unfolded_High->GetY()[i]/g_Antiproton_number_raw_High->GetY()[i]);
        g_HighAntiprotonNumberRawUnfoldedRatio_P0VGGNN->SetPoint(i, g_Antiproton_number_unfolded_High_P0VGGNN->GetX()[i], g_Antiproton_number_unfolded_High_P0VGGNN->GetY()[i]/g_Antiproton_number_raw_High_P0VGGNN->GetY()[i]);
        g_HighAntiprotonNumberRawUnfoldedRatio_P1     ->SetPoint(i, g_Antiproton_number_unfolded_High_P1->GetX()[i]     , g_Antiproton_number_unfolded_High_P1->GetY()[i]     /g_Antiproton_number_raw_High_P1->GetY()[i]);
        g_HighAntiprotonNumberRawUnfoldedRatio_P2     ->SetPoint(i, g_Antiproton_number_unfolded_High_P2->GetX()[i]     , g_Antiproton_number_unfolded_High_P2->GetY()[i]     /g_Antiproton_number_raw_High_P2->GetY()[i]);
        g_HighAntiprotonNumberRawUnfoldedRatio_P4     ->SetPoint(i, g_Antiproton_number_unfolded_High_P4->GetX()[i]     , g_Antiproton_number_unfolded_High_P4->GetY()[i]     /g_Antiproton_number_raw_High_P4->GetY()[i]);
    }


    g_LowAntiprotonNumberRawUnfoldedRatio         ->Draw("AP");
    g_IntermediateAntiprotonNumberRawUnfoldedRatio->Draw("same P");
    g_HighAntiprotonNumberRawUnfoldedRatio        ->Draw("same P");

    g_HighAntiprotonNumberRawUnfoldedRatio_P0VGGNN ->Draw("same P");
    g_HighAntiprotonNumberRawUnfoldedRatio_P1      ->Draw("same P");
    g_HighAntiprotonNumberRawUnfoldedRatio_P2      ->Draw("same P");
    g_HighAntiprotonNumberRawUnfoldedRatio_P4      ->Draw("same P");

    g_LowAntiprotonNumberRawUnfoldedRatio         ->SetTitle("");
    g_LowAntiprotonNumberRawUnfoldedRatio->SetMarkerStyle(15);
    g_LowAntiprotonNumberRawUnfoldedRatio->SetMarkerColor(2);
    g_IntermediateAntiprotonNumberRawUnfoldedRatio->SetMarkerStyle(15);
    g_IntermediateAntiprotonNumberRawUnfoldedRatio->SetMarkerColor(2);
    g_HighAntiprotonNumberRawUnfoldedRatio        ->SetMarkerStyle(15);
    g_HighAntiprotonNumberRawUnfoldedRatio        ->SetMarkerColor(2);

    g_HighAntiprotonNumberRawUnfoldedRatio_P0VGGNN ->SetMarkerStyle(15);
    g_HighAntiprotonNumberRawUnfoldedRatio_P0VGGNN ->SetMarkerColor(3);
    g_HighAntiprotonNumberRawUnfoldedRatio_P1      ->SetMarkerStyle(15);
    g_HighAntiprotonNumberRawUnfoldedRatio_P1      ->SetMarkerColor(4);
    g_HighAntiprotonNumberRawUnfoldedRatio_P2      ->SetMarkerStyle(15);
    g_HighAntiprotonNumberRawUnfoldedRatio_P2      ->SetMarkerColor(5);
    g_HighAntiprotonNumberRawUnfoldedRatio_P4      ->SetMarkerStyle(15);
    g_HighAntiprotonNumberRawUnfoldedRatio_P4      ->SetMarkerColor(6);

    TLegend *legend = new TLegend(0.20, 0.15, 0.60, 0.42);

    legend->AddEntry(g_HighAntiprotonNumberRawUnfoldedRatio_P0VGGNN   , "P0","p");
    legend->AddEntry(g_HighAntiprotonNumberRawUnfoldedRatio_P1        , "P1","p");
    legend->AddEntry(g_HighAntiprotonNumberRawUnfoldedRatio_P2        , "P2","p");
    legend->AddEntry(g_HighAntiprotonNumberRawUnfoldedRatio_P4        , "P4","p");
    legend->SetTextSize(0.04);
    legend->SetTextFont(62);
    legend->Draw();


    TAxis * xaxis = g_LowAntiprotonNumberRawUnfoldedRatio->GetXaxis();
    TAxis * yaxis = g_LowAntiprotonNumberRawUnfoldedRatio->GetYaxis();
    xaxis->SetLimits(1, 600); //1, 600
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("Raw Counts / Unfolded Counts");
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelFont(62);
    yaxis->SetLabelSize(0.05);

    gPad->SetLogx();

    c_number_ratio.SaveAs( (string("NumberRawUnfoldRatioPlot_") + issversion + string(".pdf")).c_str());    

}


void Plot_Antiproton_and_Proton_Numbers_LowRange(TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber, TGraph *g_Antiproton_number_Low, TGraph *g_Proton_number_Low, std::string issversion){

    //// Plot Antiproton numbers
    TCanvas c_number("c_number","c_number",1000,500);

    g_PhysicsReport_pbarNumber      ->Draw("AP");
    g_Antiproton_number_Low->Draw("same P");

    TAxis * xaxis_number = g_PhysicsReport_pbarNumber->GetXaxis();
    TAxis * yaxis_number = g_PhysicsReport_pbarNumber->GetYaxis();
    xaxis_number->SetLimits(0.8, 6); //1, 600
    yaxis_number->SetRangeUser(1, 90000); //1000, 90000
    //yaxis_number->SetMoreLogLabels();
    xaxis_number->SetTitle("Rigidity (GV)");
    yaxis_number->SetTitle("Antiproton Numbers");
    xaxis_number->SetTitleFont(62);
    yaxis_number->SetTitleFont(62);
    xaxis_number->SetTitleSize(0.045);
    yaxis_number->SetTitleSize(0.045);
    xaxis_number->SetLabelFont(62);
    xaxis_number->SetLabelSize(0.05);
    yaxis_number->SetLabelFont(62);
    yaxis_number->SetLabelSize(0.05);

    g_Antiproton_number_Low->SetMarkerStyle(15);
    g_Antiproton_number_Low->SetMarkerColor(6);
    g_Antiproton_number_Low->SetMarkerSize(1.5);
    g_PhysicsReport_pbarNumber->SetMarkerStyle(15);
    g_PhysicsReport_pbarNumber->SetMarkerColor(0);  //3
    g_PhysicsReport_pbarNumber->SetMarkerSize(0);

    g_PhysicsReport_pbarNumber->SetTitle("");
    gPad->SetLogy();
    //gPad->SetLogx();

    c_number.SaveAs( (string("NumberPlot_LowRange_") + issversion + string(".pdf")).c_str());


    //// Plot Proton Numbers
    TCanvas c_P_number("","",1000,500);

    g_Proton_number_Low->Draw("AP");

    TAxis * xaxis_P_number = g_Proton_number_Low->GetXaxis();
    TAxis * yaxis_P_number = g_Proton_number_Low->GetYaxis();
    xaxis_P_number->SetLimits(0.8, 6);      // 1   , 600
    yaxis_P_number->SetRangeUser(1, 1000000000); // 1000, 90000

    xaxis_P_number->SetTitle("Rigidity (GV)");
    yaxis_P_number->SetTitle("Proton Numbers");
    xaxis_P_number->SetTitleFont(62);
    yaxis_P_number->SetTitleFont(62);
    xaxis_P_number->SetTitleSize(0.045);
    yaxis_P_number->SetTitleSize(0.045);
    xaxis_P_number->SetLabelFont(62);
    xaxis_P_number->SetLabelSize(0.05);
    yaxis_P_number->SetLabelFont(62);
    yaxis_P_number->SetLabelSize(0.05);

    g_Proton_number_Low->SetMarkerStyle(15);
    g_Proton_number_Low->SetMarkerColor(6);
    g_Proton_number_Low->SetMarkerSize(1.5);

    g_Proton_number_Low->SetTitle("");
    gPad->SetLogy();

    c_P_number.SaveAs( (string("ProtonNumberPlot_LowRange_") + issversion + string(".pdf")).c_str());


    // Both Proton And Pbar Numbers
    TCanvas c_PAndPbar_number("","",1000,500);

    g_Proton_number_Low->Draw("AP");
    g_Antiproton_number_Low->Draw("same P");

    TAxis * xaxis_PAndPbar_number = g_Proton_number_Low->GetXaxis();
    TAxis * yaxis_PAndPbar_number = g_Proton_number_Low->GetYaxis();
    xaxis_PAndPbar_number->SetLimits(0.8, 6);      // 1   , 600
    yaxis_PAndPbar_number->SetRangeUser(1, 1000000000); // 1000, 90000
    xaxis_PAndPbar_number->SetTitle("|R| / GV");
    yaxis_PAndPbar_number->SetTitle("N");
    xaxis_PAndPbar_number->SetTitleFont(62);
    yaxis_PAndPbar_number->SetTitleFont(62);
    xaxis_PAndPbar_number->SetTitleSize(0.045);
    yaxis_PAndPbar_number->SetTitleSize(0.045);
    xaxis_PAndPbar_number->SetLabelFont(62);
    xaxis_PAndPbar_number->SetLabelSize(0.05);
    yaxis_PAndPbar_number->SetLabelFont(62);
    yaxis_PAndPbar_number->SetLabelSize(0.05);

    g_Proton_number_Low->SetMarkerStyle(15);
    g_Proton_number_Low->SetMarkerColor(4);
    g_Proton_number_Low->SetMarkerSize(1.5);

    g_Antiproton_number_Low->SetMarkerStyle(15);
    g_Antiproton_number_Low->SetMarkerColor(6);
    g_Antiproton_number_Low->SetMarkerSize(1.5);

    g_Proton_number_Low->SetTitle("");
    gPad->SetLogy();

    TLegend *legend_BothNumbers_Low = new TLegend(0.7, 0.2, 0.89, 0.4);
    legend_BothNumbers_Low->AddEntry(g_Proton_number_Low    , "Proton"    , "p");
    legend_BothNumbers_Low->AddEntry(g_Antiproton_number_Low, "Antiproton", "p");
    legend_BothNumbers_Low->SetTextSize(0.04);
    legend_BothNumbers_Low->SetTextFont(62);
    legend_BothNumbers_Low->SetBorderSize(0);
    legend_BothNumbers_Low->Draw();

    c_PAndPbar_number.SaveAs( (string("NumberPbarAndProtonPlot_LowRange_") + issversion + string(".pdf")).c_str());
  
}


void Plot_Antiproton_and_Proton_Numbers_IntermediateRange(TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber, TGraph *g_Antiproton_number_Intermediate, TGraph *g_Proton_number_Intermediate, std::string issversion){


    //// Plot Antiproton Numbers
    TCanvas c_number("c_number","c_number",1000,500);

    g_PhysicsReport_pbarNumber               ->Draw("AP");
    g_Antiproton_number_Intermediate->Draw("same P");

    TAxis * xaxis_number = g_PhysicsReport_pbarNumber->GetXaxis();
    TAxis * yaxis_number = g_PhysicsReport_pbarNumber->GetYaxis();
    xaxis_number->SetLimits(2, 19); //1, 600
    yaxis_number->SetRangeUser(1, 90000); //1000, 90000
    //yaxis_number->SetMoreLogLabels();
    xaxis_number->SetTitle("Rigidity (GV)");
    yaxis_number->SetTitle("Antiproton Numbers");
    xaxis_number->SetTitleFont(62);
    yaxis_number->SetTitleFont(62);
    xaxis_number->SetTitleSize(0.045);
    yaxis_number->SetTitleSize(0.045);
    xaxis_number->SetLabelFont(62);
    xaxis_number->SetLabelSize(0.05);
    yaxis_number->SetLabelFont(62);
    yaxis_number->SetLabelSize(0.05);

    g_Antiproton_number_Intermediate->SetMarkerStyle(15);
    g_Antiproton_number_Intermediate->SetMarkerColor(6);
    g_Antiproton_number_Intermediate->SetMarkerSize(1.5);
    g_PhysicsReport_pbarNumber->SetMarkerStyle(15);
    g_PhysicsReport_pbarNumber->SetMarkerColor(0);//3
    g_PhysicsReport_pbarNumber->SetMarkerSize(0);

    g_PhysicsReport_pbarNumber->SetTitle("");
    gPad->SetLogy();

    c_number.SaveAs( (string("NumberPlot_IntermediateRange_") + issversion + string(".pdf")).c_str());


    //// Both numbers
    TCanvas c_PAndPbar_number("","",1000,500);

    g_Proton_number_Intermediate     ->Draw("AP");
    g_Antiproton_number_Intermediate ->Draw("same P");

    TAxis * xaxis_PAndPbar_number = g_Proton_number_Intermediate->GetXaxis();
    TAxis * yaxis_PAndPbar_number = g_Proton_number_Intermediate->GetYaxis();
    xaxis_PAndPbar_number->SetLimits(2, 19);      // 1   , 600
    yaxis_PAndPbar_number->SetRangeUser(1, 1000000000); // 1000, 90000
    xaxis_PAndPbar_number->SetTitle("|R| / GV");
    yaxis_PAndPbar_number->SetTitle("N");
    xaxis_PAndPbar_number->SetTitleFont(62);
    yaxis_PAndPbar_number->SetTitleFont(62);
    xaxis_PAndPbar_number->SetTitleSize(0.045);
    yaxis_PAndPbar_number->SetTitleSize(0.045);
    xaxis_PAndPbar_number->SetLabelFont(62);
    xaxis_PAndPbar_number->SetLabelSize(0.05);
    yaxis_PAndPbar_number->SetLabelFont(62);
    yaxis_PAndPbar_number->SetLabelSize(0.05);

    g_Proton_number_Intermediate->SetMarkerStyle(15);
    g_Proton_number_Intermediate->SetMarkerColor(4);
    g_Proton_number_Intermediate->SetMarkerSize(1.5);

    g_Antiproton_number_Intermediate->SetMarkerStyle(15);
    g_Antiproton_number_Intermediate->SetMarkerColor(6);
    g_Antiproton_number_Intermediate->SetMarkerSize(1.5);

    g_Proton_number_Intermediate->SetTitle("");
    gPad->SetLogy();

    TLegend *legend_BothNumbers_Intermediate = new TLegend(0.7, 0.2, 0.89, 0.4);
    legend_BothNumbers_Intermediate->AddEntry(g_Proton_number_Intermediate    , "Proton"    , "p");
    legend_BothNumbers_Intermediate->AddEntry(g_Antiproton_number_Intermediate, "Antiproton", "p");
    legend_BothNumbers_Intermediate->SetTextSize(0.04);
    legend_BothNumbers_Intermediate->SetTextFont(62);
    legend_BothNumbers_Intermediate->SetBorderSize(0);
    legend_BothNumbers_Intermediate->Draw();

    c_PAndPbar_number.SaveAs( (string("NumberPbarAndProtonPlot_IntermediateRange_") + issversion + string(".pdf")).c_str());

}

void  Plot_Antiproton_and_Proton_Numbers_HighRange(TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber, TGraph *g_Antiproton_number_unfolded_High, std::string issversion){

    TCanvas c_number("c_number","c_number",1000,500);

    //g_PhysicsReport_pbarNumber->Draw("AP");
    //g_Antiproton_number_unfolded_High->Draw("same P");

    g_Antiproton_number_unfolded_High->Draw("AP");

    //TAxis * xaxis_number = g_PhysicsReport_pbarNumber->GetXaxis();
    //TAxis * yaxis_number = g_PhysicsReport_pbarNumber->GetYaxis();
    TAxis * xaxis_number = g_Antiproton_number_unfolded_High->GetXaxis();
    TAxis * yaxis_number = g_Antiproton_number_unfolded_High->GetYaxis();

    xaxis_number->SetLimits(10, 600); //1, 600
    yaxis_number->SetRangeUser(1, 90000); //1000, 90000
    xaxis_number->SetMoreLogLabels();
    //yaxis_number->SetMoreLogLabels();
    xaxis_number->SetTitle("Rigidity (GV)");
    yaxis_number->SetTitle("Antiproton Numbers");
    xaxis_number->SetTitleFont(62);
    yaxis_number->SetTitleFont(62);
    xaxis_number->SetTitleSize(0.045);
    yaxis_number->SetTitleSize(0.045);
    xaxis_number->SetLabelFont(62);
    xaxis_number->SetLabelSize(0.05);
    yaxis_number->SetLabelFont(62);
    yaxis_number->SetLabelSize(0.05);
    xaxis_number->SetTitleOffset(1.5);

    gPad->SetBottomMargin(0.15);

    g_Antiproton_number_unfolded_High->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High->SetMarkerColor(6);
    g_Antiproton_number_unfolded_High->SetMarkerSize(0.9);
    g_PhysicsReport_pbarNumber->SetMarkerStyle(15);
    g_PhysicsReport_pbarNumber->SetMarkerColor(0);//3
    g_PhysicsReport_pbarNumber->SetMarkerSize(0);

    g_PhysicsReport_pbarNumber->SetTitle("");
    g_Antiproton_number_unfolded_High->SetTitle("");
    gPad->SetLogx();
    gPad->SetLogy();

    c_number.SaveAs( (string("NumberPlot_HighRange_") + issversion + string(".pdf")).c_str());

}

void Plot_Antiproton_and_Proton_Numbers_Ratio(TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber, TGraph *g_Antiproton_number_unfolded_Low, TGraph *g_Proton_number_unfolded_Low, TGraph *g_Antiproton_number_unfolded_Intermediate, TGraph *g_Proton_number_unfolded_Intermediate, TGraph *g_Antiproton_number_unfolded_High, TGraph *g_Antiproton_number_unfolded_High_P0, TGraph *g_Antiproton_number_unfolded_High_P1, TGraph *g_Antiproton_number_unfolded_High_P2, TGraph *g_Antiproton_number_unfolded_High_P3, TGraph *g_Antiproton_number_unfolded_High_P4, TGraph *g_Antiproton_number_unfolded_High_P5, TGraph *g_Proton_number_unfolded_High, std::string issversion, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin){


    TCanvas c_number("c_number","c_number",1000,500);

    //// Calculate PbarNumberRatio to Reference
    TGraph *g_ExpectedPbarNumber_fromPhysicsReport = new TGraph(*g_PhysicsReport_pbarNumber);

    // Calculate Expected PbarNumbers from Reference.
    for (int q = 0; q < g_PhysicsReport_pbarNumber->GetN(); q++){
        //g_ExpectedPbarNumber_fromPhysicsReport->SetPoint(q, g_PhysicsReport_pbarNumber->GetX()[q], g_PhysicsReport_pbarNumber->GetY()[q]/6.5*9);
        g_ExpectedPbarNumber_fromPhysicsReport->SetPoint(q, g_PhysicsReport_pbarNumber->GetX()[q], g_PhysicsReport_pbarNumber->GetY()[q]);
    }        

    // Calculate PbarNumber Ratio to Expection
    TGraph *g_PbarNumberRatio_withExpectedFromPhysicsReport = new TGraph(*g_PhysicsReport_pbarNumber);

    for (int q = 0; q < g_Antiproton_number_unfolded_Low->GetN(); q++){
        g_PbarNumberRatio_withExpectedFromPhysicsReport->SetPoint( q, g_Antiproton_number_unfolded_Low->GetX()[q], g_Antiproton_number_unfolded_Low->GetY()[q] / g_ExpectedPbarNumber_fromPhysicsReport->GetY()[q] );
    }


    for (int q = 0; q < g_Antiproton_number_unfolded_Intermediate->GetN(); q++){
        g_PbarNumberRatio_withExpectedFromPhysicsReport->SetPoint( q+g_Antiproton_number_unfolded_Low->GetN(), g_Antiproton_number_unfolded_Intermediate->GetX()[q], g_Antiproton_number_unfolded_Intermediate->GetY()[q]/g_ExpectedPbarNumber_fromPhysicsReport->GetY()[q+g_Antiproton_number_unfolded_Low->GetN()]);
    }


    for (int q = 0; q < g_Antiproton_number_unfolded_High->GetN(); q++){
        g_PbarNumberRatio_withExpectedFromPhysicsReport->SetPoint( q+g_Antiproton_number_unfolded_Low->GetN()+g_Antiproton_number_unfolded_Intermediate->GetN(), g_Antiproton_number_unfolded_High->GetX()[q], g_Antiproton_number_unfolded_High->GetY()[q]/g_ExpectedPbarNumber_fromPhysicsReport->GetY()[q+g_Antiproton_number_unfolded_Low->GetN()+g_Antiproton_number_unfolded_Intermediate->GetN()]);
    }


    g_PbarNumberRatio_withExpectedFromPhysicsReport->Draw("AP");

    TAxis * xaxis_number = g_PbarNumberRatio_withExpectedFromPhysicsReport->GetXaxis();
    TAxis * yaxis_number = g_PbarNumberRatio_withExpectedFromPhysicsReport->GetYaxis();
    xaxis_number->SetLimits(1, 600);
    yaxis_number->SetRangeUser(0.1, 10);
    xaxis_number->SetTitle("Rigidity (GV)");
    yaxis_number->SetTitle("Antiproton Numbers Ratio (This analysis / PhyRep)");
    xaxis_number->SetTitleFont(62);
    yaxis_number->SetTitleFont(62);
    xaxis_number->SetTitleSize(0.045);
    yaxis_number->SetTitleSize(0.045);
    xaxis_number->SetLabelFont(62);
    xaxis_number->SetLabelSize(0.05);
    yaxis_number->SetLabelFont(62);
    yaxis_number->SetLabelSize(0.05);

    g_PbarNumberRatio_withExpectedFromPhysicsReport->SetMarkerStyle(15);
    g_PbarNumberRatio_withExpectedFromPhysicsReport->SetMarkerColor(2);

    gPad->SetLogy();
    gPad->SetLogx();

    c_number.SaveAs( (string("NumberRatioPlot_") + issversion + string(".pdf")).c_str());
}


void Plot_AntiprotonNumbersAllPattern(TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber, TGraph *g_Antiproton_number_unfolded_Low, TGraph *g_Proton_number_unfolded_Low, TGraph *g_Antiproton_number_unfolded_Intermediate, TGraph *g_Proton_number_unfolded_Intermediate, TGraph *g_Antiproton_number_unfolded_High, TGraph *g_Antiproton_number_unfolded_High_P0, TGraph *g_Antiproton_number_unfolded_High_P1, TGraph *g_Antiproton_number_unfolded_High_P2, TGraph *g_Antiproton_number_unfolded_High_P3, TGraph *g_Antiproton_number_unfolded_High_P4, TGraph *g_Antiproton_number_unfolded_High_P5, TGraph *g_Proton_number_unfolded_High, std::string issversion){

    TCanvas c_number("c_number","c_number",1000,500);

    g_Antiproton_number_unfolded_High_P0->Draw("AP");
    g_Antiproton_number_unfolded_High_P1->Draw("same P");
    g_Antiproton_number_unfolded_High_P2->Draw("same P");
    //g_Antiproton_number_unfolded_High_P3->Draw("same P");
    g_Antiproton_number_unfolded_High_P4->Draw("same P");
    //g_Antiproton_number_unfolded_High_P5->Draw("same P");

    TAxis * xaxis_number = g_Antiproton_number_unfolded_High_P0->GetXaxis();
    TAxis * yaxis_number = g_Antiproton_number_unfolded_High_P0->GetYaxis();
    xaxis_number->SetLimits(10, 600);
    yaxis_number->SetRangeUser(1, 90000);
    xaxis_number->SetTitle("Rigidity (GV)");
    yaxis_number->SetTitle("Antiproton Numbers");
    xaxis_number->SetTitleFont(62);
    yaxis_number->SetTitleFont(62);
    xaxis_number->SetTitleSize(0.045);
    yaxis_number->SetTitleSize(0.045);
    xaxis_number->SetLabelFont(62);
    xaxis_number->SetLabelSize(0.05);
    yaxis_number->SetLabelFont(62);
    yaxis_number->SetLabelSize(0.05);

    g_Antiproton_number_unfolded_High_P0->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P0->SetMarkerColor(7);
    g_Antiproton_number_unfolded_High_P1->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P1->SetMarkerColor(41);
    g_Antiproton_number_unfolded_High_P2->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P2->SetMarkerColor(8);
    g_Antiproton_number_unfolded_High_P3->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P3->SetMarkerColor(5);
    g_Antiproton_number_unfolded_High_P4->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P4->SetMarkerColor(28);
    g_Antiproton_number_unfolded_High_P5->SetMarkerStyle(15);
    g_Antiproton_number_unfolded_High_P5->SetMarkerColor(50);

    g_Antiproton_number_unfolded_High_P0->SetTitle("");

    gPad->SetLogy();
    gPad->SetLogx();

    TLegend *legend_number = new TLegend(0.20, 0.15, 0.60, 0.42);
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P0, "This analysis (HighRange)(P0)","p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P1, "This analysis (HighRange)(P1)", "p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P2, "This analysis (HighRange)(P2)", "p");
    //legend_number->AddEntry(g_Antiproton_number_unfolded_High_P3, "This analysis (HighRange)(P3)", "p");
    legend_number->AddEntry(g_Antiproton_number_unfolded_High_P4, "This analysis (HighRange)(P4)", "p");
    //legend_number->AddEntry(g_Antiproton_number_unfolded_High_P5, "This analysis (HighRange)(P5)", "p");
    legend_number->SetTextSize(0.04);
    legend_number->SetTextFont(62);
    legend_number->Draw();

    TLine line1(38.9, 1, 38.9, 90000);
    line1.Draw("same");
    TLine line2(147, 1, 147, 90000);
    line2.Draw("same");
    TLine line3(175, 1, 175, 90000);
    line3.Draw("same");
    line1.SetLineStyle(2);
    line2.SetLineStyle(2);

    c_number.SaveAs( (string("NumberPlot_AllPatterns") + issversion + string(".pdf")).c_str());
}


void Plot_AntiprotonNumbersOverRigidityBinWidth_High(
TGraph *g_Antiproton_number_High_P0, TGraph *g_Antiproton_number_High_P1, TGraph *g_Antiproton_number_High_P2, TGraph *g_Antiproton_number_High_P3, TGraph *g_Antiproton_number_High_P4, TGraph *g_Antiproton_number_High_P5, 
TH1D *h_Antiproton_number_High_P0, TH1D *h_Antiproton_number_High_P1, TH1D *h_Antiproton_number_High_P2, TH1D *h_Antiproton_number_High_P3, TH1D *h_Antiproton_number_High_P4, TH1D *h_Antiproton_number_High_P5, 
TGraph *g_Proton_number_High_P0, TGraph *g_Proton_number_High_P1, TGraph *g_Proton_number_High_P2, TGraph *g_Proton_number_High_P3, TGraph *g_Proton_number_High_P4, TGraph *g_Proton_number_High_P5,
std::string issversion, int HighRangeRemovedAtBegin){

    //// index for used range
    int index_38  = 13 - HighRangeRemovedAtBegin;
    int index_147 = 27 - HighRangeRemovedAtBegin;
    int index_175 = 28 - HighRangeRemovedAtBegin;

    //// Calculating Antiproton and Proton number/RigidityBinWidth
    TGraph *g_Antiproton_number_High_P0_Over_RigidityBin_used = new TGraph(g_Antiproton_number_High_P0->GetN());
    TGraph *g_Antiproton_number_High_P1_Over_RigidityBin_used = new TGraph(index_147);
    TGraph *g_Antiproton_number_High_P2_Over_RigidityBin_used = new TGraph(index_175);
    TGraph *g_Antiproton_number_High_P4_Over_RigidityBin_used = new TGraph(index_38);
    TGraph *g_Antiproton_number_High_P1_Over_RigidityBin_notused = new TGraph(g_Antiproton_number_High_P0->GetN()-index_147);
    TGraph *g_Antiproton_number_High_P2_Over_RigidityBin_notused = new TGraph(g_Antiproton_number_High_P0->GetN()-index_175);
    TGraph *g_Antiproton_number_High_P4_Over_RigidityBin_notused = new TGraph(g_Antiproton_number_High_P0->GetN()-index_38);

    TGraph *g_Proton_number_High_P0_Over_RigidityBin_used = new TGraph(g_Proton_number_High_P0->GetN());
    TGraph *g_Proton_number_High_P1_Over_RigidityBin_used = new TGraph(index_147);
    TGraph *g_Proton_number_High_P2_Over_RigidityBin_used = new TGraph(index_175);
    TGraph *g_Proton_number_High_P4_Over_RigidityBin_used = new TGraph(index_38);
    TGraph *g_Proton_number_High_P1_Over_RigidityBin_notused = new TGraph(g_Proton_number_High_P0->GetN()-index_147);
    TGraph *g_Proton_number_High_P2_Over_RigidityBin_notused = new TGraph(g_Proton_number_High_P0->GetN()-index_175);
    TGraph *g_Proton_number_High_P4_Over_RigidityBin_notused = new TGraph(g_Proton_number_High_P0->GetN()-index_38);


    //Used Area 
    for (int i=0; i<g_Antiproton_number_High_P0->GetN(); i++){
        g_Antiproton_number_High_P0_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_High_P0->GetX()[i], g_Antiproton_number_High_P0->GetY()[i]/h_Antiproton_number_High_P0->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_High_P0_Over_RigidityBin_used    ->SetPoint(i, g_Proton_number_High_P0->GetX()[i]    , g_Proton_number_High_P0->GetY()[i]/h_Antiproton_number_High_P0->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }
    for (int i=0; i<index_147; i++){
        g_Antiproton_number_High_P1_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_High_P1->GetX()[i], g_Antiproton_number_High_P1->GetY()[i]/h_Antiproton_number_High_P1->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_High_P1_Over_RigidityBin_used->SetPoint(i, g_Proton_number_High_P1->GetX()[i], g_Proton_number_High_P1->GetY()[i]/h_Antiproton_number_High_P1->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }
    for (int i=0; i<index_175; i++){
        g_Antiproton_number_High_P2_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_High_P2->GetX()[i], g_Antiproton_number_High_P2->GetY()[i]/h_Antiproton_number_High_P2->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_High_P2_Over_RigidityBin_used->SetPoint(i, g_Proton_number_High_P2->GetX()[i], g_Proton_number_High_P2->GetY()[i]/h_Antiproton_number_High_P2->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }
    for (int i=0; i<index_38; i++){
        g_Antiproton_number_High_P4_Over_RigidityBin_used->SetPoint(i, g_Antiproton_number_High_P4->GetX()[i], g_Antiproton_number_High_P4->GetY()[i]/h_Antiproton_number_High_P4->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_High_P4_Over_RigidityBin_used->SetPoint(i, g_Proton_number_High_P4->GetX()[i], g_Proton_number_High_P4->GetY()[i]/h_Antiproton_number_High_P4->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }

    //Not Used Area
    for (int i=index_147; i<g_Antiproton_number_High_P1->GetN(); i++){
        g_Antiproton_number_High_P1_Over_RigidityBin_notused->SetPoint(i, g_Antiproton_number_High_P1->GetX()[i], g_Antiproton_number_High_P1->GetY()[i]/h_Antiproton_number_High_P1->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_High_P1_Over_RigidityBin_notused->SetPoint(i, g_Proton_number_High_P1->GetX()[i], g_Proton_number_High_P1->GetY()[i]/h_Antiproton_number_High_P1->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }
    for (int i=index_175; i<g_Antiproton_number_High_P2->GetN(); i++){
        g_Antiproton_number_High_P2_Over_RigidityBin_notused->SetPoint(i, g_Antiproton_number_High_P2->GetX()[i], g_Antiproton_number_High_P2->GetY()[i]/h_Antiproton_number_High_P2->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_High_P2_Over_RigidityBin_notused->SetPoint(i, g_Proton_number_High_P2->GetX()[i], g_Proton_number_High_P2->GetY()[i]/h_Antiproton_number_High_P2->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }
    for (int i=index_38; i<g_Antiproton_number_High_P4->GetN(); i++){
        g_Antiproton_number_High_P4_Over_RigidityBin_notused->SetPoint(i, g_Antiproton_number_High_P4->GetX()[i], g_Antiproton_number_High_P4->GetY()[i]/h_Antiproton_number_High_P4->GetBinWidth(i+1+HighRangeRemovedAtBegin));
        g_Proton_number_High_P4_Over_RigidityBin_notused->SetPoint(i, g_Proton_number_High_P4->GetX()[i], g_Proton_number_High_P4->GetY()[i]/h_Antiproton_number_High_P4->GetBinWidth(i+1+HighRangeRemovedAtBegin));
    }

    // Plot 
    TCanvas c_number_OverRigidityBinWidth("c_number_OverRigidityBinWidth","c_number_OverRigidityBinWidth",1000,500);

    g_Antiproton_number_High_P0_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Antiproton_number_High_P0_Over_RigidityBin_used->SetMarkerColor(46);
    g_Antiproton_number_High_P0_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Antiproton_number_High_P1_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Antiproton_number_High_P1_Over_RigidityBin_used->SetMarkerColor(41);
    g_Antiproton_number_High_P1_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Antiproton_number_High_P1_Over_RigidityBin_notused->SetMarkerStyle(4);
    g_Antiproton_number_High_P1_Over_RigidityBin_notused->SetMarkerColor(41);
    g_Antiproton_number_High_P1_Over_RigidityBin_notused->SetMarkerSize(0.9);
    g_Antiproton_number_High_P2_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Antiproton_number_High_P2_Over_RigidityBin_used->SetMarkerColor(38);
    g_Antiproton_number_High_P2_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Antiproton_number_High_P2_Over_RigidityBin_notused->SetMarkerStyle(4);
    g_Antiproton_number_High_P2_Over_RigidityBin_notused->SetMarkerColor(38);
    g_Antiproton_number_High_P2_Over_RigidityBin_notused->SetMarkerSize(0.9);
    g_Antiproton_number_High_P4_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Antiproton_number_High_P4_Over_RigidityBin_used->SetMarkerColor(49);
    g_Antiproton_number_High_P4_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Antiproton_number_High_P4_Over_RigidityBin_notused->SetMarkerStyle(4);
    g_Antiproton_number_High_P4_Over_RigidityBin_notused->SetMarkerColor(49);
    g_Antiproton_number_High_P4_Over_RigidityBin_notused->SetMarkerSize(0.9);

    g_Proton_number_High_P0_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Proton_number_High_P0_Over_RigidityBin_used->SetMarkerColor(8);
    g_Proton_number_High_P0_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Proton_number_High_P1_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Proton_number_High_P1_Over_RigidityBin_used->SetMarkerColor(9);
    g_Proton_number_High_P1_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Proton_number_High_P1_Over_RigidityBin_notused->SetMarkerStyle(4);
    g_Proton_number_High_P1_Over_RigidityBin_notused->SetMarkerColor(9);
    g_Proton_number_High_P1_Over_RigidityBin_notused->SetMarkerSize(0.9);
    g_Proton_number_High_P2_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Proton_number_High_P2_Over_RigidityBin_used->SetMarkerColor(12);
    g_Proton_number_High_P2_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Proton_number_High_P2_Over_RigidityBin_notused->SetMarkerStyle(4);
    g_Proton_number_High_P2_Over_RigidityBin_notused->SetMarkerColor(12);
    g_Proton_number_High_P2_Over_RigidityBin_notused->SetMarkerSize(0.9);
    g_Proton_number_High_P4_Over_RigidityBin_used->SetMarkerStyle(15);
    g_Proton_number_High_P4_Over_RigidityBin_used->SetMarkerColor(30);
    g_Proton_number_High_P4_Over_RigidityBin_used->SetMarkerSize(0.9);
    g_Proton_number_High_P4_Over_RigidityBin_notused->SetMarkerStyle(4);
    g_Proton_number_High_P4_Over_RigidityBin_notused->SetMarkerColor(30);
    g_Proton_number_High_P4_Over_RigidityBin_notused->SetMarkerSize(0.9);

    g_Antiproton_number_High_P0_Over_RigidityBin_used->Draw("AP");
    g_Antiproton_number_High_P1_Over_RigidityBin_used->Draw("same P");
    g_Antiproton_number_High_P2_Over_RigidityBin_used->Draw("same P");
    g_Antiproton_number_High_P4_Over_RigidityBin_used->Draw("same P");
    g_Proton_number_High_P0_Over_RigidityBin_used->Draw("same P");
    g_Proton_number_High_P1_Over_RigidityBin_used->Draw("same P");
    g_Proton_number_High_P2_Over_RigidityBin_used->Draw("same P");
    g_Proton_number_High_P4_Over_RigidityBin_used->Draw("same P");
    //g_Antiproton_number_High_P1_Over_RigidityBin_notused->Draw("same P");
    //g_Antiproton_number_High_P2_Over_RigidityBin_notused->Draw("same P");
    //g_Antiproton_number_High_P4_Over_RigidityBin_notused->Draw("same P");

    g_Antiproton_number_High_P0_Over_RigidityBin_used->SetTitle("");

    TAxis * xaxis_number_Over_RigidityBin = g_Antiproton_number_High_P0_Over_RigidityBin_used->GetXaxis();
    TAxis * yaxis_number_Over_RigidityBin = g_Antiproton_number_High_P0_Over_RigidityBin_used->GetYaxis();
    xaxis_number_Over_RigidityBin->SetLimits(10, 600);
    yaxis_number_Over_RigidityBin->SetRangeUser(0.1, 15000000000); //0.1, 90000
    xaxis_number_Over_RigidityBin->SetTitle("|R| / GV");
    yaxis_number_Over_RigidityBin->SetTitle("N / #Delta R");
    xaxis_number_Over_RigidityBin->SetTitleFont(62);
    yaxis_number_Over_RigidityBin->SetTitleFont(62);
    xaxis_number_Over_RigidityBin->SetTitleSize(0.045);
    yaxis_number_Over_RigidityBin->SetTitleSize(0.045);
    xaxis_number_Over_RigidityBin->SetLabelFont(62);
    xaxis_number_Over_RigidityBin->SetLabelSize(0.05);
    yaxis_number_Over_RigidityBin->SetLabelFont(62);
    yaxis_number_Over_RigidityBin->SetLabelSize(0.05);
    xaxis_number_Over_RigidityBin->SetMoreLogLabels();

    xaxis_number_Over_RigidityBin->SetTitleOffset(1.5);
    gPad->SetBottomMargin(0.15);

    TLine line1_OverRigidityBin(38.9, 0.1, 38.9, 10000000000);
    line1_OverRigidityBin.Draw("same");
    TLine line2_OverRigidityBin(147, 0.1, 147, 10000000);
    line2_OverRigidityBin.Draw("same");
    TLine line3_OverRigidityBin(175, 0.1, 175, 10000000);
    line3_OverRigidityBin.Draw("same");
    line1_OverRigidityBin.SetLineStyle(2);
    line2_OverRigidityBin.SetLineStyle(2);
    line3_OverRigidityBin.SetLineStyle(2);

    TLegend *legend_number_OverRigidityBin = new TLegend(0.40, 0.7, 0.85, 0.89);
    legend_number_OverRigidityBin->AddEntry(g_Antiproton_number_High_P0_Over_RigidityBin_used, "FullSpan (Antiproton)"  ,"p"); // Pattern 0 
    legend_number_OverRigidityBin->AddEntry(g_Proton_number_High_P0_Over_RigidityBin_used    , "FullSpan (Proton)"  ,"p");
    legend_number_OverRigidityBin->AddEntry(g_Antiproton_number_High_P1_Over_RigidityBin_used, "Inner + L1 (Antiproton)","p"); // Pattern 1
    legend_number_OverRigidityBin->AddEntry(g_Proton_number_High_P1_Over_RigidityBin_used    , "Inner + L1 (Proton)","p");
    legend_number_OverRigidityBin->AddEntry(g_Antiproton_number_High_P2_Over_RigidityBin_used, "Inner + L9 (Antiproton)","p"); // Pattern 2
    legend_number_OverRigidityBin->AddEntry(g_Proton_number_High_P2_Over_RigidityBin_used    , "Inner + L9 (Proton)","p");
    legend_number_OverRigidityBin->AddEntry(g_Antiproton_number_High_P4_Over_RigidityBin_used, "Inner Only (Antiproton)","p"); // Pattern 4
    legend_number_OverRigidityBin->AddEntry(g_Proton_number_High_P4_Over_RigidityBin_used    , "Inner Only (Proton)","p");
    legend_number_OverRigidityBin->SetTextSize(0.04);
    legend_number_OverRigidityBin->SetTextFont(62);
    legend_number_OverRigidityBin->SetBorderSize(0);
    legend_number_OverRigidityBin->SetNColumns(2);
    legend_number_OverRigidityBin->Draw();

    gPad->SetLogy();
    gPad->SetLogx();

    c_number_OverRigidityBinWidth.SaveAs( (string("NumberPlot_AllPatternsOverRigidityBinWidth") + issversion + string(".pdf")).c_str());
}


void Plot_AntiprotonNumbersToRrferencesRatio_InHighRange(TGraph *g_PhysicsReport_pbarNumber, TGraph *g_Published_pbarNumber, TGraph *g_Antiproton_number_unfolded_Low, TGraph *g_Proton_number_unfolded_Low, TGraph *g_Antiproton_number_unfolded_Intermediate, TGraph *g_Proton_number_unfolded_Intermediate, TGraph *g_Antiproton_number_unfolded_High, TGraph *g_Antiproton_number_unfolded_High_P0, TGraph *g_Antiproton_number_unfolded_High_P1, TGraph *g_Antiproton_number_unfolded_High_P2, TGraph *g_Antiproton_number_unfolded_High_P4, TGraph *g_Proton_number_unfolded_High, std::string issversion){

    TCanvas c_number("c_number","c_number",1000,500);

    //// Calculate PbarNumberRatio to Reference
    TGraph *g_ExpectedPbarNumber_fromPhysicsReport = new TGraph(*g_PhysicsReport_pbarNumber);

    // Calculate Expected PbarNumbers from Reference.
    for (int q = 0; q < g_PhysicsReport_pbarNumber->GetN(); q++){
        g_ExpectedPbarNumber_fromPhysicsReport->SetPoint(q, g_PhysicsReport_pbarNumber->GetX()[q], g_PhysicsReport_pbarNumber->GetY()[q]/6.5*9);
    }

    // Calculate PbarNumber Ratio to Expection
    int LowAndInterNumber = g_ExpectedPbarNumber_fromPhysicsReport->GetN() - g_Antiproton_number_unfolded_High->GetN();
    TGraph *g_PbarNumberRatio_withExpectedFromPhysicsReport = new TGraph(*g_Antiproton_number_unfolded_High);
    for (int q = 0; q < g_Antiproton_number_unfolded_High->GetN(); q++){
        if (g_Antiproton_number_unfolded_High->GetX()[q] != g_ExpectedPbarNumber_fromPhysicsReport->GetX()[q+LowAndInterNumber]){
           cout<< "X binning don't match, check below:" << endl;
           cout << g_Antiproton_number_unfolded_High->GetX()[q] << "_" << g_ExpectedPbarNumber_fromPhysicsReport->GetX()[q+LowAndInterNumber] <<endl; 
        }
        g_PbarNumberRatio_withExpectedFromPhysicsReport->SetPoint( q, g_Antiproton_number_unfolded_High->GetX()[q], g_Antiproton_number_unfolded_High->GetY()[q] / g_ExpectedPbarNumber_fromPhysicsReport->GetY()[q+LowAndInterNumber] );
    }

    //// Plot PbarNumberRatio
    g_PbarNumberRatio_withExpectedFromPhysicsReport->Draw("AP");

    TAxis * xaxis_number = g_PbarNumberRatio_withExpectedFromPhysicsReport->GetXaxis();
    TAxis * yaxis_number = g_PbarNumberRatio_withExpectedFromPhysicsReport->GetYaxis();
    xaxis_number->SetLimits(10,600);
    //yaxis_number->SetRangeUser(1,90000);
    xaxis_number->SetTitle("Rigidity (GV)");
    yaxis_number->SetTitle("Antiproton Numbers Ratio");
    xaxis_number->SetTitleFont(62);
    yaxis_number->SetTitleFont(62);
    xaxis_number->SetTitleSize(0.045);
    yaxis_number->SetTitleSize(0.045);
    xaxis_number->SetLabelFont(62);
    xaxis_number->SetLabelSize(0.05);
    yaxis_number->SetLabelFont(62);
    yaxis_number->SetLabelSize(0.05);
    gPad->SetLogx();
    xaxis_number->SetMoreLogLabels();

    g_PbarNumberRatio_withExpectedFromPhysicsReport->SetMarkerStyle(15);
    g_PbarNumberRatio_withExpectedFromPhysicsReport->SetMarkerColor(2);

    g_PbarNumberRatio_withExpectedFromPhysicsReport->SetTitle("");
    /*
    TLegend *legend_number = new TLegend(0.25,0.15,0.78,0.36);
    legend_number->AddEntry(g_PbarNumberRatio_withExpectedFromPhysicsReport              , "PbarRatioToPhysicsReport", "p");
    legend_number->SetTextSize(0.04);
    legend_number->SetTextFont(62);
    legend_number->Draw();
    */

    //TLine line1(38.9, 0, 38.9, 1);
    //line1.Draw("same");
    //TLine line2(147, 0, 147, 1);
    //line2.Draw("same");
    //TLine line3(175, 0, 175, 1);
    //line3.Draw("same");
    //line1.SetLineStyle(2);
    //line2.SetLineStyle(2);
    //line3.SetLineStyle(2);

    c_number.SaveAs( (string("AntiprotonNumberRatioPlot_") + issversion + string(".pdf")).c_str());
}


// Acceptance Plots
void Plot_EffectiveAcceptanceInThreeRanges(TGraphErrors *Effective_Acceptance_High_B1042_P0, TGraphErrors *Effective_Acceptance_High_B1042_P1, TGraphErrors *Effective_Acceptance_High_B1042_P2, TGraphErrors *Effective_Acceptance_High_B1042_P4, TGraphErrors *Effective_Acceptance_High, TGraphErrors *Effective_Acceptance_Low, TGraphErrors *Effective_Acceptance_Intermediate, TGraphErrors *Effective_Acceptance_Intermediate_inHigh){
    //Note: Effective_Acceptance_Intermediate_inHigh only has 2 points, others are nan.

    // Fit with more complelx functions
    TF1  *functionP0 = new TF1("function1","[0]*log(log(log(x)))+[1]",14,525);
    Effective_Acceptance_High->Fit(functionP0, "", "", 14, 525);
    TF1  *functionP1 = new TF1("function1","[0]*log(log(log(x)))+[1]",14,525);
    Effective_Acceptance_High_B1042_P1->Fit(functionP1, "", "", 14, 525);   
    TF1  *functionP2 = new TF1("function1","[0]*log(log(log(x)))+[1]",14,525);
    Effective_Acceptance_High_B1042_P2->Fit(functionP2, "", "", 14, 525);
    TF1  *functionP4 = new TF1("function1","[0]*log(log(log(x)))+[1]",14,525);
    Effective_Acceptance_High_B1042_P4->Fit(functionP4, "", "", 14, 525); 

    TCanvas c_acceptance("c_acceptance","c_acceptance",1000,500);

    // remove fit in acceptance
    //Effective_Acceptance_Low                ->GetListOfFunctions()->Remove(Effective_Acceptance_Low                ->GetFunction("function1"));
    //Effective_Acceptance_Intermediate       ->GetListOfFunctions()->Remove(Effective_Acceptance_Intermediate       ->GetFunction("function1"));
    //Effective_Acceptance_High               ->GetListOfFunctions()->Remove(Effective_Acceptance_High               ->GetFunction("function1"));
    //Effective_Acceptance_Intermediate_inHigh->GetListOfFunctions()->Remove(Effective_Acceptance_Intermediate_inHigh->GetFunction("function1"));

    // Fit Line color
    Effective_Acceptance_Low          ->GetFunction("function1")->SetLineColor(2);
    Effective_Acceptance_Intermediate ->GetFunction("function1")->SetLineColor(2);
    Effective_Acceptance_High         ->GetFunction("function1")->SetLineColor(8);
    Effective_Acceptance_High_B1042_P1->GetFunction("function1")->SetLineColor(9);
    Effective_Acceptance_High_B1042_P2->GetFunction("function1")->SetLineColor(12);
    Effective_Acceptance_High_B1042_P4->GetFunction("function1")->SetLineColor(30);

    // Fit Line range
    Effective_Acceptance_Low         ->GetFunction("function1")->SetRange(1.0 , 4.43);
    Effective_Acceptance_Intermediate->GetFunction("function1")->SetRange(4.43, 15.3);
    Effective_Acceptance_High        ->GetFunction("function1")->SetRange(15.3, 525);

    Effective_Acceptance_Low                ->Draw("AP");
    Effective_Acceptance_Intermediate       ->Draw("same P");
    //Effective_Acceptance_Intermediate_inHigh->Draw("same P");
    Effective_Acceptance_High               ->Draw("same P");

    //Effective_Acceptance_High_B1042_P0 ->Draw("same P");
    Effective_Acceptance_High_B1042_P1 ->Draw("same P");
    Effective_Acceptance_High_B1042_P2 ->Draw("same P");
    Effective_Acceptance_High_B1042_P4 ->Draw("same P");

    // remove statistics box 
    TPaveStats *pavestats = (TPaveStats*)Effective_Acceptance_High->FindObject("stats");
    TPaveStats *pavestats_P1 = (TPaveStats*)Effective_Acceptance_High_B1042_P1->FindObject("stats");
    TPaveStats *pavestats_P2 = (TPaveStats*)Effective_Acceptance_High_B1042_P2->FindObject("stats");
    TPaveStats *pavestats_P4 = (TPaveStats*)Effective_Acceptance_High_B1042_P4->FindObject("stats");
    pavestats->SetX1NDC(2);   // new x start position
    pavestats->SetX2NDC(2.1); // new y start position
    pavestats->SetY1NDC(2);   // new x end position
    pavestats->SetY2NDC(2.1); // new y end position
    pavestats_P1->SetX1NDC(2);  
    pavestats_P1->SetX2NDC(2.1);
    pavestats_P1->SetY1NDC(2);
    pavestats_P1->SetY2NDC(2.1);
    pavestats_P2->SetX1NDC(2);
    pavestats_P2->SetX2NDC(2.1);
    pavestats_P2->SetY1NDC(2);
    pavestats_P2->SetY2NDC(2.1);
    pavestats_P4->SetX1NDC(2);
    pavestats_P4->SetX2NDC(2.1);
    pavestats_P4->SetY1NDC(2);
    pavestats_P4->SetY2NDC(2.1);


    TAxis * xaxis_acceptance = Effective_Acceptance_Low->GetXaxis();
    TAxis * yaxis_acceptance = Effective_Acceptance_Low->GetYaxis();
    xaxis_acceptance->SetLimits(1.0, 600);
    yaxis_acceptance->SetRangeUser(1, 1.35);
    xaxis_acceptance->SetTitle("|R|_{true} / (GV)");
    yaxis_acceptance->SetTitle("A_{p} / A_{p}"); // SetTitle("#frac{A_{p}}{A_{#bar{p}}}");
    xaxis_acceptance->SetTitleFont(62);
    yaxis_acceptance->SetTitleFont(62);
    xaxis_acceptance->SetTitleSize(0.045);
    yaxis_acceptance->SetTitleSize(0.045);
    xaxis_acceptance->SetLabelFont(62);
    xaxis_acceptance->SetLabelSize(0.05);
    yaxis_acceptance->SetLabelFont(62);
    yaxis_acceptance->SetLabelSize(0.05);
    xaxis_acceptance->SetLabelOffset(0);
    yaxis_acceptance->SetLabelOffset(0.02);
    xaxis_acceptance->SetTitleOffset(1.2);
    yaxis_acceptance->SetTitleOffset(0);

    TLatex latex;
    latex.SetNDC(1);
    latex.SetTextSize(0.055);
    latex.SetTextAngle(90);
    latex.SetTextAlign(13);  //align at top
    latex.DrawLatex(0.063, 0.88, "#minus");

    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.13);
    gStyle->SetOptFit(0);

    Effective_Acceptance_Low                ->SetMarkerStyle(15);
    Effective_Acceptance_Low                ->SetMarkerColor(2);
    Effective_Acceptance_Low                ->SetMarkerSize(0.9);

    Effective_Acceptance_Intermediate       ->SetMarkerStyle(15);
    Effective_Acceptance_Intermediate       ->SetMarkerColor(2);
    Effective_Acceptance_Intermediate       ->SetMarkerSize(0.9);

    Effective_Acceptance_Intermediate_inHigh->SetMarkerStyle(15);
    Effective_Acceptance_Intermediate_inHigh->SetMarkerColor(0);
    Effective_Acceptance_Intermediate_inHigh->SetMarkerSize(0.9);

    Effective_Acceptance_High               ->SetMarkerStyle(15);
    Effective_Acceptance_High               ->SetMarkerColor(8);
    Effective_Acceptance_High               ->SetMarkerSize(0.9);
    Effective_Acceptance_High_B1042_P1->SetMarkerStyle(15);
    Effective_Acceptance_High_B1042_P1->SetMarkerColor(9);
    Effective_Acceptance_High_B1042_P1->SetMarkerSize(0.9);
    Effective_Acceptance_High_B1042_P2->SetMarkerStyle(15);
    Effective_Acceptance_High_B1042_P2->SetMarkerColor(12);
    Effective_Acceptance_High_B1042_P2->SetMarkerSize(0.9);
    Effective_Acceptance_High_B1042_P4->SetMarkerStyle(15);
    Effective_Acceptance_High_B1042_P4->SetMarkerColor(30);
    Effective_Acceptance_High_B1042_P4->SetMarkerSize(0.9);

    Effective_Acceptance_Low->SetTitle("");
    gPad->SetLogx();
    gStyle->SetOptStat(0);

    TLegend *legend_acceptance = new TLegend(0.5, 0.65, 0.8, 0.85);
    legend_acceptance->AddEntry(Effective_Acceptance_Low            , "InnerCentral"     , "p");
    //legend_acceptance->AddEntry(Effective_Acceptance_Intermediate , "intermediate range", "p");
    legend_acceptance->AddEntry(Effective_Acceptance_High           , "FullSpan"        , "p");
    legend_acceptance->AddEntry(Effective_Acceptance_High_B1042_P1  , "Inner + L1"        , "p");
    legend_acceptance->AddEntry(Effective_Acceptance_High_B1042_P2  , "Inner + L9"        , "p");
    legend_acceptance->AddEntry(Effective_Acceptance_High_B1042_P4  , "Inner Only"        , "p");
    legend_acceptance->SetTextSize(0.04);
    legend_acceptance->SetTextFont(62);
    legend_acceptance->SetBorderSize(0);
    legend_acceptance->Draw();

    c_acceptance.SaveAs( (string("AcceptancePlot") + string(".pdf")).c_str());
}



// Ratio Plots  (overlapped range included)
void Plot_Ratio_ThisAnalysisOnly_Overlapped(TGraphErrors *gPublishedRatio, TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, TGraphErrors *g_HighResult, std::string issversionname){
    TCanvas c1("c1","c1",1000,500);
    //TPad *pad1 = new TPad("pad1", "pad1", 0.03, 0.62, 0.50, 0.92, 32);
    TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.0, 1.0, 1.0);
    pad1->Draw("");
    pad1->cd();

    //// Remove the first two points in interrmediate range for illustration.
    g_IntermediateResult->RemovePoint(0);
    g_IntermediateResult->RemovePoint(0);

    //// Plot
    g_HighResult->Draw("AP");
    //gPublishedRatio->Draw("same P");
    g_LowResult->Draw("same P");
    g_IntermediateResult->Draw("same P");
    g_HighResult->Draw("same P"); // Replot it ,to make it on top of Published Ratio.


    TAxis * xaxis = g_HighResult->GetXaxis();
    TAxis * yaxis = g_HighResult->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    //yaxis->SetTitle("Antiproton to Proton ratio");
    //yaxis->SetTitle("#frac{#bar{p}}{p}");
    yaxis->SetTitle("#bar{p}/p");
    //yaxis->CenterTitle(true);
    xaxis->SetLimits(1.0,600);
    yaxis->SetRangeUser(0.000000001, 0.00028);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.07);
    yaxis->SetTitleSize(0.07);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.07);
    yaxis->SetLabelSize(0.07);
    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(0.9);

    //gPad->SetLogy();
    gPad->SetLogx();
    xaxis->SetMoreLogLabels();
    g_HighResult->SetTitle("");

    g_HighResult->SetMarkerStyle(15);
    g_HighResult->SetMarkerColor(28);
    g_HighResult->SetLineColor(28);
    g_LowResult->SetMarkerStyle(15);
    g_LowResult->SetMarkerColor(6);
    g_LowResult->SetLineColor(6);
    g_IntermediateResult->SetMarkerStyle(15);
    g_IntermediateResult->SetMarkerColor(4);
    g_IntermediateResult->SetLineColor(4);
    gPublishedRatio->SetMarkerStyle(15);
    gPublishedRatio->SetMarkerSize(0.7);
    gPublishedRatio->SetMarkerColor(1);
    gPublishedRatio->SetLineWidth(1);

    TLegend *legend1 = new TLegend(0.42, 0.22, 0.82, 0.42);
    legend1->AddEntry(g_LowResult         , "Low Rigidity Range"         , "p");
    legend1->AddEntry(g_IntermediateResult, "Intermediate Rigidity Range", "p");
    legend1->AddEntry(g_HighResult        , "High Rigidity Range"        , "p");
    //legend1->AddEntry(gPublishedRatio,"PRL paper 2016","lpf");
    legend1->SetTextSize(0.05);
    legend1->SetTextFont(62);
    legend1->SetBorderSize(0);
    legend1->Draw();

    pad1->SetBottomMargin(0.18);
    pad1->SetLeftMargin(0.12);

    c1.SaveAs( (string("fullratio_") + issversionname + string(".pdf")).c_str());
}


void Plot_RawRatioAndUnfoldedRatio_Compare(std::string issversion, TGraphErrors *g_Raw_Unfold_Compare_P0VGGNN, TGraphErrors *g_Raw_Unfold_Compare_P1, TGraphErrors *g_Raw_Unfold_Compare_P2, TGraphErrors *g_Raw_Unfold_Compare_P4, TGraphErrors *g_RawLowResult, TGraphErrors *g_LowResult, TGraphErrors *g_RawIntermediateResult, TGraphErrors *g_IntermediateResult, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin){

    //// index for used range
    int index_38  = 13;
    int index_147 = 27;
    int index_175 = 28;

    TGraphErrors *g_Raw_Unfold_Compare_P0VGGNN_used = new TGraphErrors(g_Raw_Unfold_Compare_P0VGGNN->GetN());
    TGraphErrors *g_Raw_Unfold_Compare_P1_used      = new TGraphErrors(index_147);
    TGraphErrors *g_Raw_Unfold_Compare_P2_used      = new TGraphErrors(index_175);
    TGraphErrors *g_Raw_Unfold_Compare_P4_used      = new TGraphErrors(index_38);

    //Used Area
    for (int i=0; i<g_Raw_Unfold_Compare_P0VGGNN->GetN(); i++){
        g_Raw_Unfold_Compare_P0VGGNN_used->SetPoint     (i, g_Raw_Unfold_Compare_P0VGGNN->GetX()[i], g_Raw_Unfold_Compare_P0VGGNN->GetY()[i]);
        g_Raw_Unfold_Compare_P0VGGNN_used->SetPointError(i, 0                                      , g_Raw_Unfold_Compare_P0VGGNN->GetErrorY(i));
    }
    for (int i=0; i<index_147; i++){
        g_Raw_Unfold_Compare_P1_used     ->SetPoint     (i, g_Raw_Unfold_Compare_P1->GetX()[i], g_Raw_Unfold_Compare_P1->GetY()[i]);
        g_Raw_Unfold_Compare_P1_used     ->SetPointError(i, 0                                 , g_Raw_Unfold_Compare_P1->GetErrorY(i)); 
    }
    for (int i=0; i<index_175; i++){
        g_Raw_Unfold_Compare_P2_used     ->SetPoint     (i, g_Raw_Unfold_Compare_P2->GetX()[i], g_Raw_Unfold_Compare_P2->GetY()[i]);
        g_Raw_Unfold_Compare_P2_used     ->SetPointError(i, 0                                 , g_Raw_Unfold_Compare_P2->GetErrorY(i));
    }
    for (int i=0; i<index_38; i++){
        g_Raw_Unfold_Compare_P4_used     ->SetPoint     (i, g_Raw_Unfold_Compare_P4->GetX()[i], g_Raw_Unfold_Compare_P4->GetY()[i]);
        g_Raw_Unfold_Compare_P4_used     ->SetPointError(i, 0                                 , g_Raw_Unfold_Compare_P4->GetErrorY(i));
    }

    g_Raw_Unfold_Compare_P4_used->SetPoint(index_38-1, g_Raw_Unfold_Compare_P4->GetX()[index_38-1], g_Raw_Unfold_Compare_P4->GetY()[index_38-2]);
    
    //// Calculate Raw and Unfolded Ratio In: Low And Intermediate Range
    // Low Range
    TGraphErrors *g_Raw_Unfold_Compare_Low = new TGraphErrors(g_LowResult->GetN()-3);
    for (int i = 0; i < 6; i=i+2) {
        g_Raw_Unfold_Compare_Low ->SetPoint     (i, (g_LowResult->GetX()[i]+g_LowResult->GetX()[i+1])/2, ( ((g_RawLowResult->GetY()[i] - g_LowResult->GetY()[i])/g_RawLowResult->GetY()[i]) + ((g_RawLowResult->GetY()[i+1] - g_LowResult->GetY()[i+1])/g_RawLowResult->GetY()[i+1]) )/2 );
        g_Raw_Unfold_Compare_Low ->SetPointError(i, 0                     ,(sqrt( pow(g_RawLowResult->GetErrorY(i)/g_LowResult->GetY()[i], 2) + pow(g_RawLowResult->GetY()[i]*g_LowResult->GetErrorY(i)/pow(g_LowResult->GetY()[i],2), 2) ) ) / 1.414 );
    }
    for (int i = 6; i < g_LowResult->GetN(); ++i) {
        g_Raw_Unfold_Compare_Low ->SetPoint     (i, g_LowResult->GetX()[i], (g_RawLowResult->GetY()[i] - g_LowResult->GetY()[i])/g_RawLowResult->GetY()[i] );
        g_Raw_Unfold_Compare_Low ->SetPointError(i, 0                     , (sqrt( pow(g_RawLowResult->GetErrorY(i)/g_LowResult->GetY()[i], 2) + pow(g_RawLowResult->GetY()[i]*g_LowResult->GetErrorY(i)/pow(g_LowResult->GetY()[i],2), 2) ) ) );
    }
    /*
    for (int i = 0; i < g_LowResult->GetN(); ++i) {
        g_Raw_Unfold_Compare_Low ->SetPoint     (i, g_LowResult->GetX()[i], (g_RawLowResult->GetY()[i] - g_LowResult->GetY()[i])/g_RawLowResult->GetY()[i] );
        g_Raw_Unfold_Compare_Low ->SetPointError(i, 0                     , (sqrt( pow(g_RawLowResult->GetErrorY(i)/g_LowResult->GetY()[i], 2) + pow(g_RawLowResult->GetY()[i]*g_LowResult->GetErrorY(i)/pow(g_LowResult->GetY()[i],2), 2) ) ) );
    }
    */
    // Intermediate Range
    TGraphErrors *g_Raw_Unfold_Compare_Intermediate = new TGraphErrors(g_RawIntermediateResult->GetN());
    for (int i = 0; i < g_RawIntermediateResult->GetN(); ++i) {
        g_Raw_Unfold_Compare_Intermediate ->SetPoint     (i, g_IntermediateResult->GetX()[i], (g_RawIntermediateResult->GetY()[i] - g_IntermediateResult->GetY()[i])/g_RawIntermediateResult->GetY()[i] );
        g_Raw_Unfold_Compare_Intermediate ->SetPointError(i, 0                              , (sqrt( pow(g_RawIntermediateResult->GetErrorY(i)/g_IntermediateResult->GetY()[i], 2) + pow(g_RawIntermediateResult->GetY()[i]*g_IntermediateResult->GetErrorY(i)/pow(g_IntermediateResult->GetY()[i],2), 2) ) ) );
    }

    
    for (int i=0; i < LowRangeRemoveAtEnd; i++){
        g_Raw_Unfold_Compare_Low->RemovePoint(g_Raw_Unfold_Compare_Low->GetN()-1);        
    }
    for (int i=0; i < IntermediateRangeRemoveAtBegin; i++){
        g_Raw_Unfold_Compare_Intermediate->RemovePoint(0);
    }
    for (int i=0; i < IntermediateRangeRemoveAtEnd; i++){
        g_Raw_Unfold_Compare_Intermediate->RemovePoint(g_Raw_Unfold_Compare_Intermediate->GetN()-1);
    }
    for (int i=0; i < HighRangeRemovedAtBegin; i++){
        g_Raw_Unfold_Compare_P0VGGNN_used->RemovePoint(0);
        g_Raw_Unfold_Compare_P1_used     ->RemovePoint(0);
        g_Raw_Unfold_Compare_P2_used     ->RemovePoint(0);
        g_Raw_Unfold_Compare_P4_used     ->RemovePoint(0);
    }
    
    //// Plot

    TCanvas c1("c1","c1",1000,500);
    TPad *p1 = new TPad("p1", "", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1->Draw();
    p1->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gPad->SetLogx();

    g_Raw_Unfold_Compare_P0VGGNN_used->Draw("AP");
    g_Raw_Unfold_Compare_P1_used->Draw("same P");
    g_Raw_Unfold_Compare_P2_used->Draw("same P");
    g_Raw_Unfold_Compare_P4_used->Draw("same P");

    g_Raw_Unfold_Compare_Low         ->Draw("same P");
    g_Raw_Unfold_Compare_Intermediate->Draw("same P");

    g_Raw_Unfold_Compare_P0VGGNN_used->SetTitle("");

    g_Raw_Unfold_Compare_P0VGGNN_used->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P0VGGNN_used->SetMarkerColor(7);
    g_Raw_Unfold_Compare_P1_used->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P1_used->SetMarkerColor(41);
    g_Raw_Unfold_Compare_P2_used->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P2_used->SetMarkerColor(8);
    g_Raw_Unfold_Compare_P4_used->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_P4_used->SetMarkerColor(28);
    g_Raw_Unfold_Compare_Low->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_Low->SetMarkerColor(2);
    g_Raw_Unfold_Compare_Intermediate->SetMarkerStyle(15);
    g_Raw_Unfold_Compare_Intermediate->SetMarkerColor(4);

    TAxis * xaxis = g_Raw_Unfold_Compare_P0VGGNN_used->GetXaxis();
    TAxis * yaxis = g_Raw_Unfold_Compare_P0VGGNN_used->GetYaxis();
    xaxis->SetTitle("Rigidity (GV)");
    yaxis->SetTitle("#frac{Raw Ratio - Unfolded Ratio}{Unfolded Ratio}");
    xaxis->SetMoreLogLabels();
    xaxis->SetNoExponent();
    xaxis->SetLimits(1.0, 600);
    yaxis->SetRangeUser(-0.5, 0.8);

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

    //gPad->SetLogx();
    //xaxis->SetMoreLogLabels();

    //gPad->Update();
    //c1.Update();

    TLegend *legend = new TLegend(0.27, 0.60, 0.65, 0.88);
    legend->AddEntry(g_Raw_Unfold_Compare_Low         , "Low Rigidity Range"              ,"p");
    legend->AddEntry(g_Raw_Unfold_Compare_Intermediate, "Intermediate Rigidity Range"     ,"p");
    legend->AddEntry(g_Raw_Unfold_Compare_P0VGGNN_used, "High Rigidity Range (FullSpan)"  ,"p"); // Pattern 0
    legend->AddEntry(g_Raw_Unfold_Compare_P1_used     , "High Rigidity Range (Inner + L1)","p"); // Pattern 1
    legend->AddEntry(g_Raw_Unfold_Compare_P2_used     , "High Rigidity Range (Inner + L9)","p"); // Pattern 2
    legend->AddEntry(g_Raw_Unfold_Compare_P4_used     , "High Rigidity Range (Inner Only)","p"); // Pattern 4
    legend->SetTextSize(0.04);
    legend->SetTextFont(62);
    legend->Draw();

    c1.SaveAs((string("RawUnfoldedRatioCompare") + string("_") + issversion + string(".pdf")).c_str());
}


// Ratio Plots (Compare with two refenrences)  (overlapped range NOT included)
void Plot_Ratio_CompareWithReference_NotOverlapped(TGraphErrors *gPublishedRatio, TGraphErrors *gPhyrReortRatio, TGraphErrors *g_HighResult, TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, std::string issversionname, std::string issversion){

    vector <std::string> Referencelegendname = {"PRL paper 2016", "AMS02 published (May 2011 to Nov 2017)"};
    vector <std::string> Referencefilename   = {"PRLpaper2016"  , "PhysicsReport2021"  };

    // Loop for two reference
    for (int plotindex=0; plotindex<2; plotindex++){
        TGraphErrors *gReferenceRatio;
        if (plotindex == 0){
            gReferenceRatio = gPublishedRatio;}
        else if (plotindex == 1){
            gReferenceRatio = gPhyrReortRatio;}

        TCanvas c_ratio("c_ratio","c_ratio",1000,500);
        // Workaround Option to stwich off the residual plot (step 1):
        //TPad *pad_ratio      = new TPad("pad_ratio"   , "The pad1 ratio"   , 0.0, 0.0  , 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
        TPad *pad_ratio    = new TPad("pad_ratio"   , "The pad1 ratio"   , 0.0, 0.245, 1.0, 1.0 , 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
        TPad *pad_residual = new TPad("pad_residual", "The pad2 residaul", 0.0, 0.0  , 1.0, 0.35, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)

        pad_ratio->Draw();
        // Workaround Option to stwich off the residual plot (step 2):
        pad_residual->Draw();
        
        // upper plot: ratio plot
        pad_ratio->cd();
        //gStyle->SetPadBorderMode(0);
        //gStyle->SetFrameBorderMode(0);

        gReferenceRatio->Draw("AP"); // Option: "Z" is to remove the small horizontal line in vertical error bar.

        g_HighResult        ->Draw("same P");
        //g_HighResult_P0   ->Draw("same P");
        //g_HighResult_P1   ->Draw("same P");
        //g_HighResult_P2   ->Draw("same P");
        //g_HighResult_P4   ->Draw("same P");
        g_LowResult         ->Draw("same P");
        g_IntermediateResult->Draw("same P");
        g_HighResult        ->Draw("same P"); // Replot it ,to make it on top of Published Ratio.

        TAxis * xaxis_ratio = gReferenceRatio->GetXaxis();
        TAxis * yaxis_ratio = gReferenceRatio->GetYaxis();
        yaxis_ratio->SetTitle("#bar{p}/p");
        xaxis_ratio->SetTitle("|R| / (GV)");
        xaxis_ratio->SetLimits(1.0,600);
        yaxis_ratio->SetRangeUser(0.000000001, 0.00032);
        xaxis_ratio->SetTitleFont(62);
        yaxis_ratio->SetTitleFont(62);
        xaxis_ratio->SetTitleSize(0.045);
        yaxis_ratio->SetTitleSize(0.065);
        xaxis_ratio->SetLabelFont(62);
        yaxis_ratio->SetLabelFont(62);
        xaxis_ratio->SetLabelSize(0.05);
        yaxis_ratio->SetLabelSize(0.05);

        xaxis_ratio->SetTitleOffset(1.3);
        yaxis_ratio->SetTitleOffset(0.6);
        pad_ratio->SetBottomMargin(0.18);
        pad_ratio->SetLeftMargin(0.1);

        // Workaround Option to stwich off the residual plot (step 3):
        xaxis_ratio->SetLabelOffset(999); // remove x axis
        xaxis_ratio->SetLabelSize(0);     // remove x label

        gPad->SetLogx();
        xaxis_ratio->SetMoreLogLabels();
        g_HighResult->SetTitle("");
        gReferenceRatio->SetTitle("");

        g_HighResult->SetMarkerStyle(15);
        g_HighResult->SetMarkerColor(2);
        g_HighResult->SetMarkerSize(0.9);
        g_LowResult->SetMarkerStyle(15);
        g_LowResult->SetMarkerColor(2);
        g_LowResult->SetLineColor(2);
        g_LowResult->SetMarkerSize(0.9);
        g_IntermediateResult->SetMarkerStyle(15);
        g_IntermediateResult->SetMarkerColor(2);
        g_IntermediateResult->SetMarkerSize(0.9);
        gReferenceRatio->SetMarkerStyle(24);
        gReferenceRatio->SetMarkerSize(0.9);
        gReferenceRatio->SetMarkerColor(1);
        gReferenceRatio->SetLineWidth(1);
        gReferenceRatio->SetLineColor(1);

        double xlow, xhigh, ylow, yhigh;
        if (issversion == "pass7.8"){
            xlow  = 0.22; 
            xhigh = 0.7;
            ylow  = 0.7;
            yhigh = 0.9;
        }
        else if (issversion == "PhyRep2021"){
            xlow  = 0.22;
            xhigh = 0.67;
            ylow  = 0.70;
            yhigh = 0.89;
        }
        else if (issversion == "2016paper"){
            xlow  = 0.22;
            xhigh = 0.7;
            ylow  = 0.7;
            yhigh = 0.9;
        }
        TLegend *legend_ratio = new TLegend(xlow, ylow, xhigh, yhigh);
        //legend_ratio->AddEntry(g_LowResult, (string("This analysis(X bin shift ") + to_string_with_precision((1-binshift)*100,0) + string("%)")).c_str(),"lp");
        if (issversion == "pass7.8"){
            legend_ratio->AddEntry(g_LowResult, (string("This analysis (All Time Range)") ).c_str(),"p");}
        else if (issversion == "2016paper"){
            legend_ratio->AddEntry(g_LowResult, (string("This analysis (PRL Time Range)") ).c_str(),"p");}
        else if (issversion == "PhyRep2021"){
            legend_ratio->AddEntry(g_LowResult, (string("This analysis (May 2011 - Nov 2017)") ).c_str(),"p");}
        //legend_ratio->AddEntry(g_HighResult_P0, "Patttern 0", "lp");
        //legend_ratio->AddEntry(g_HighResult_P1, "Patttern 1", "lp");
        //legend_ratio->AddEntry(g_HighResult_P2, "Patttern 2", "lp");
        //legend_ratio->AddEntry(g_HighResult_P4, "Patttern 4", "lp");

        legend_ratio->AddEntry(gReferenceRatio, (Referencelegendname.at(plotindex)).c_str(), "p");
        legend_ratio->SetTextSize(0.04);
        legend_ratio->SetTextFont(62);
        legend_ratio->SetBorderSize(0);
        legend_ratio->Draw();

        TLine line_break_upper(1.51, 0, 1.51, 0.00032);
        line_break_upper.SetLineStyle(2);
        line_break_upper.Draw("same");
        
        // lower plot
        pad_residual->cd();
        pad_residual->SetBottomMargin(0.3);
        pad_residual->SetLeftMargin(0.1);

        TGraphErrors *Residual_Low = new TGraphErrors(g_LowResult->GetN());
        TGraphErrors *Residual_Intermediate = new TGraphErrors(g_IntermediateResult->GetN());
        TGraphErrors *Residual_High = new TGraphErrors(g_HighResult->GetN());

        for (int e = 0; e < g_LowResult->GetN(); e++){
            Residual_Low->SetPoint(e, g_LowResult->GetX()[e], (g_LowResult->GetY()[e] - gReferenceRatio->GetY()[e])/gReferenceRatio->GetY()[e]);
            //cout<< "Residual_Low: " << (g_LowResult->GetY()[e] - gReferenceRatio->GetY()[e])/gReferenceRatio->GetY()[e] << endl;
        }

        int IntermediateStartIndex = g_LowResult->GetN(); // here, g_LowResult has already removed last overlapped points.
        for (int e = 0; e < g_IntermediateResult->GetN(); e++){
            Residual_Intermediate->SetPoint(e, g_IntermediateResult->GetX()[e], (g_IntermediateResult->GetY()[e] - gReferenceRatio->GetY()[e+IntermediateStartIndex]) / gReferenceRatio->GetY()[e+IntermediateStartIndex]);
            //cout<< "Residual_Intermediate: " << (g_IntermediateResult->GetY()[e] - gReferenceRatio->GetY()[e+IntermediateStartIndex]) / gReferenceRatio->GetY()[e+IntermediateStartIndex] <<endl;
        }

        int HighStartIndex = g_LowResult->GetN() + g_IntermediateResult->GetN();
        if ( (issversion == "pass7.8") || (issversion == "PhyRep2021") ){
            for (int e = 0; e < g_HighResult->GetN(); e++){
                Residual_High->SetPoint(e, g_HighResult->GetX()[e], (g_HighResult->GetY()[e] - gReferenceRatio->GetY()[e+HighStartIndex]) / gReferenceRatio->GetY()[e+HighStartIndex]);
            //cout<< "Residual_High: " << (g_HighResult->GetY()[e] - gReferenceRatio->GetY()[e+HighStartIndex]) / gReferenceRatio->GetY()[e+HighStartIndex] <<endl;
            }
        }
        else if (issversion == "2016paper"){
            for (int e = 0; e < g_HighResult->GetN()-2; e++){
                Residual_High->SetPoint(e, g_HighResult->GetX()[e], (g_HighResult->GetY()[e] - gReferenceRatio->GetY()[e+HighStartIndex]) / gReferenceRatio->GetY()[e+HighStartIndex]);
            }
        }

        Residual_Low->Draw("AP");

        TAxis * xaxis_residuallow = Residual_Low->GetXaxis();
        TAxis * yaxis_residuallow = Residual_Low->GetYaxis();
        xaxis_residuallow->SetLimits(1.0,600);
        yaxis_residuallow->SetRangeUser(-0.3, 0.5);
        Residual_Low->SetMarkerStyle(15);
        Residual_Low->SetMarkerColor(1);
        Residual_Low->SetMarkerSize(0.9);
        Residual_Intermediate->SetMarkerStyle(15);
        Residual_Intermediate->SetMarkerColor(1);
        Residual_Intermediate->SetMarkerSize(0.9);
        Residual_High->SetMarkerStyle(15);
        Residual_High->SetMarkerColor(1);
        Residual_High->SetMarkerSize(0.9);
        gPad->SetLogx();

        Residual_Low->SetTitle("");
        xaxis_residuallow->SetTitle("|R| / (GV)");
        if (plotindex == 0){
            yaxis_residuallow->SetTitle("#frac{This - PRL}{PRL}");
        }
        else if (plotindex == 1){
            yaxis_residuallow->SetTitle("#frac{This - Published}{Published}");
        }
        xaxis_residuallow->SetMoreLogLabels();
        xaxis_residuallow->SetTitleFont(62);
        yaxis_residuallow->SetTitleFont(62);
        xaxis_residuallow->SetTitleSize(0.12);
        yaxis_residuallow->SetTitleSize(0.11);
        xaxis_residuallow->SetLabelFont(62);
        yaxis_residuallow->SetLabelFont(62);
        xaxis_residuallow->SetLabelSize(0.12);
        yaxis_residuallow->SetLabelSize(0.1);
        xaxis_residuallow->SetTickSize(0.1);
        xaxis_residuallow->SetTitleOffset(1.2);
        yaxis_residuallow->SetTitleOffset(0.4);

        Residual_Intermediate->Draw("same P");
        Residual_High->Draw("same P");

      
        TLine line(1, 0, 600, 0);
        line.Draw("same");

        TLine line2(1, -0.1, 600, -0.1);
        line2.SetLineStyle(2);
        //line2.Draw("same");

        TLine line3(1, 0.1, 600, 0.1);
        line3.SetLineStyle(2);
        //line3.Draw("same");

        TLine line_break_lower(1.51, -0.3, 1.51, 0.5);
        line_break_lower.SetLineStyle(2);
        line_break_lower.Draw("same");
        
        c_ratio.SaveAs( (string("fullratio_NOT_overlapped_" + Referencefilename.at(plotindex) + "_") + issversionname + string(".pdf")).c_str());
    }
}


// Ratio All Patterns Plots (Compare with two refenrences)  (overlapped range NOT included)
void Plot_Ratio_AllTrackerPatterns_CompareWithReference_NotOverlapped(TGraphErrors *gPublishedRatio, TGraphErrors *gPhyrReortRatio, TGraphErrors *g_HighResult, TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, TGraphErrors *g_HighResult_P0, TGraphErrors *g_HighResult_P1, TGraphErrors *g_HighResult_P2, TGraphErrors *g_HighResult_P3, TGraphErrors *g_HighResult_P4, TGraphErrors *g_HighResult_P5, std::string issversionname){

    vector <std::string> Referencelegendname = {"PRL paper 2016", "Physics Report 2021"};
    vector <std::string> Referencefilename   = {"PRLpaper2016"  , "PhysicsReport2021"};

    for (int plotindex=0; plotindex<2; plotindex++){
        TGraphErrors *gReferenceRatio;
        if (plotindex == 0){
            gReferenceRatio = gPublishedRatio;
        }
        else if (plotindex == 1){
            gReferenceRatio = gPhyrReortRatio;
        }

        TCanvas c_ratio("c_ratio","c_ratio",1000,500);
        TPad *pad_ratio = new TPad("pad_ratio", "The pad1 ratio", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)

        pad_ratio->Draw();

        pad_ratio->cd();
        gStyle->SetPadBorderMode(0);
        gStyle->SetFrameBorderMode(0);

        gReferenceRatio->Draw("AP Z");

        g_HighResult_P0->Draw("same PZ");
        g_HighResult_P1->Draw("same PZ");
        g_HighResult_P2->Draw("same PZ");
        g_HighResult_P4->Draw("same PZ");

        g_LowResult->Draw("same PZ");
        g_IntermediateResult->Draw("same PZ");
        g_HighResult->Draw("same PZ"); // Replot it ,to make it on top of Published Ratio.

        TAxis * xaxis_ratio = gReferenceRatio->GetXaxis();
        TAxis * yaxis_ratio = gReferenceRatio->GetYaxis();
        //yaxis_ratio->SetTitle("Antiproton to Proton ratio");
        yaxis_ratio->SetTitle("#bar{p}/p");
        //yaxis_ratio->CenterTitle(true);
        xaxis_ratio->SetTitle("Rigidity (GV)");
        xaxis_ratio->SetLimits(1.0,600);
        yaxis_ratio->SetRangeUser(0.000000001, 0.00032);
        yaxis_ratio->SetTitleFont(62);
        xaxis_ratio->SetTitleSize(0.07);
        yaxis_ratio->SetTitleSize(0.07);
        xaxis_ratio->SetLabelFont(62);
        yaxis_ratio->SetLabelFont(62);
        xaxis_ratio->SetLabelSize(0.07);
        yaxis_ratio->SetLabelSize(0.07);
        xaxis_ratio->SetTitleOffset(1.3);
        yaxis_ratio->SetTitleOffset(0.9);
        pad_ratio->SetBottomMargin(0.18);
        pad_ratio->SetLeftMargin(0.18);

        gPad->SetLogx();
        xaxis_ratio->SetMoreLogLabels();
        g_HighResult->SetTitle("");
        gReferenceRatio->SetTitle("");

        g_HighResult->SetMarkerStyle(15);
        g_HighResult->SetMarkerColor(2);
        g_HighResult->SetMarkerSize(0.9);
        g_LowResult->SetMarkerStyle(15);
        g_LowResult->SetMarkerColor(2);
        g_LowResult->SetLineColor(2);
        g_LowResult->SetMarkerSize(0.9);
        g_IntermediateResult->SetMarkerStyle(15);
        g_IntermediateResult->SetMarkerColor(2);
        g_IntermediateResult->SetMarkerSize(0.9);
        gReferenceRatio->SetMarkerStyle(24);
        gReferenceRatio->SetMarkerSize(0.9);
        gReferenceRatio->SetMarkerColor(1);
        gReferenceRatio->SetLineWidth(1);
        gReferenceRatio->SetLineColor(1);

        g_HighResult_P0->SetMarkerStyle(15);
        g_HighResult_P0->SetMarkerColor(3);
        g_HighResult_P0->SetMarkerSize(0.9);
        g_HighResult_P1->SetMarkerStyle(15);
        g_HighResult_P1->SetMarkerColor(4);
        g_HighResult_P1->SetMarkerSize(0.9);
        g_HighResult_P2->SetMarkerStyle(15);
        g_HighResult_P2->SetMarkerColor(46);
        g_HighResult_P2->SetMarkerSize(0.9);
        g_HighResult_P4->SetMarkerStyle(15);
        g_HighResult_P4->SetMarkerColor(6);
        g_HighResult_P4->SetMarkerSize(0.9);


        TLegend *legend_ratio = new TLegend(0.22, 0.68, 0.5, 0.89);
        legend_ratio->AddEntry(g_LowResult, (string("All Patterns") ).c_str(),"lp");
        legend_ratio->AddEntry(g_HighResult_P0, "Pattern 0", "lp");
        legend_ratio->AddEntry(g_HighResult_P1, "Pattern 1", "lp");
        legend_ratio->AddEntry(g_HighResult_P2, "Pattern 2", "lp");
        legend_ratio->AddEntry(g_HighResult_P4, "Pattern 4", "lp");
        legend_ratio->AddEntry(gReferenceRatio, (Referencelegendname.at(plotindex)).c_str(), "lp");
        legend_ratio->SetTextSize(0.04);
        legend_ratio->SetTextFont(62);
        legend_ratio->Draw();

        c_ratio.SaveAs( (string("fullratio_AllPatterns_NOT_overlapped_" + Referencefilename.at(plotindex) + "_") + issversionname + string(".pdf")).c_str());
    }
}


// Ratio Plots (PRL paper and Physics Report)
void Plot_Ratio_CompareBetweenReferences(TGraphErrors *gPhyrReortRatio, TGraphErrors *gPublishedRatio){
    TCanvas c_twoRef("c_twoRef","c_twoRef",1000,500);
    TPad *pad_twoRefratio = new TPad("pad_twoRefratio", "The pad1 ratio", 0.0, 0.245, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    TPad *pad_twoRefresidual = new TPad("pad_twoRefresidual", "The pad2 residaul", 0.0, 0.0, 1.0, 0.3, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)

    pad_twoRefratio->Draw();
    pad_twoRefresidual->Draw();

    // upper plot: ratio plot
    pad_twoRefratio->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    gPhyrReortRatio->Draw("AP");
    gPublishedRatio->Draw("same P");

    TAxis * xaxis_twoRefratio = gPhyrReortRatio->GetXaxis();
    TAxis * yaxis_twoRefratio = gPhyrReortRatio->GetYaxis();
    xaxis_twoRefratio->SetTitle("");
    //yaxis_twoRefratio->SetTitle("Antiproton to Proton ratio");
    yaxis_twoRefratio->SetTitle("#bar{p}/p");
    //yaxis_twoRefratio->CenterTitle(true);
    xaxis_twoRefratio->SetLimits(1.0,600);
    yaxis_twoRefratio->SetRangeUser(0.000000001,0.00028);
    yaxis_twoRefratio->SetTitleFont(62);
    xaxis_twoRefratio->SetTitleSize(0);
    yaxis_twoRefratio->SetTitleSize(0.05);
    xaxis_twoRefratio->SetLabelFont(62);
    yaxis_twoRefratio->SetLabelFont(62);
    xaxis_twoRefratio->SetLabelSize(0);
    yaxis_twoRefratio->SetLabelSize(0.07);
    xaxis_twoRefratio->SetTitleOffset(0.8);
    yaxis_twoRefratio->SetTitleOffset(0.8);

    xaxis_twoRefratio->SetLabelOffset(999); // remove x axis
    xaxis_twoRefratio->SetLabelSize(0); // remove x label

    gPad->SetLogx();
    xaxis_twoRefratio->SetMoreLogLabels();
    gPhyrReortRatio->SetTitle("");

    gPhyrReortRatio->SetMarkerStyle(15);
    gPhyrReortRatio->SetMarkerColor(1);
    gPublishedRatio->SetMarkerStyle(15);
    gPublishedRatio->SetMarkerColor(2);

    TLegend *legend_twoRefratio = new TLegend(0.55,0.15,0.87,0.3);
    legend_twoRefratio->AddEntry(gPhyrReortRatio, "Physics Report 2021","lp");
    legend_twoRefratio->AddEntry(gPublishedRatio, "PRL paper 2016", "lp");
    legend_twoRefratio->SetTextSize(0.05);
    legend_twoRefratio->SetTextFont(62);
    legend_twoRefratio->Draw();

    // lower plot
    pad_twoRefresidual->cd();
    pad_twoRefresidual->SetBottomMargin(0.3);
    pad_twoRefresidual->SetLeftMargin(0.1);

    TGraphErrors *Residual_twoRef = new TGraphErrors(gPublishedRatio->GetN());

    for (int e = 0; e < gPublishedRatio->GetN()-2; e++){
        Residual_twoRef->SetPoint(e, gPublishedRatio->GetX()[e], (gPhyrReortRatio->GetY()[e] - gPublishedRatio->GetY()[e])/gPublishedRatio->GetY()[e]);
    }

    Residual_twoRef->Draw("AP");

    TAxis * xaxis_residuallow_twoRef = Residual_twoRef->GetXaxis();
    TAxis * yaxis_residuallow_twoRef = Residual_twoRef->GetYaxis();
    xaxis_residuallow_twoRef->SetLimits(1.0,600);
    //yaxis_residuallow_twoRef->SetRangeUser(-0.5,1.0);
    //yaxis_residuallow_twoRef->SetRangeUser(-0.2,0.2);
    Residual_twoRef->SetMarkerStyle(15);
    Residual_twoRef->SetMarkerColor(1);
    gPad->SetLogx();

    Residual_twoRef->SetTitle("");
    xaxis_residuallow_twoRef->SetTitle("|R| / (GV)");
    yaxis_residuallow_twoRef->SetTitle("#frac{PhyRep - PRL}{PRL}");

    xaxis_residuallow_twoRef->SetMoreLogLabels();
    xaxis_residuallow_twoRef->SetTitleFont(62);
    yaxis_residuallow_twoRef->SetTitleFont(62);
    xaxis_residuallow_twoRef->SetTitleSize(0.1);
    yaxis_residuallow_twoRef->SetTitleSize(0.1);
    xaxis_residuallow_twoRef->SetLabelFont(62);
    yaxis_residuallow_twoRef->SetLabelFont(62);
    xaxis_residuallow_twoRef->SetLabelSize(0.15);
    yaxis_residuallow_twoRef->SetLabelSize(0.1);
    yaxis_residuallow_twoRef->SetTitleOffset(0.4);
    xaxis_residuallow_twoRef->SetTickSize(0.1);

    TLine line(1, 0, 600, 0);
    line.Draw("same");
    TLine line2(1, -0.1, 600, -0.1);
    line2.SetLineStyle(2);
    line2.Draw("same");
    TLine line3(1, 0.1, 600, 0.1);
    line3.SetLineStyle(2);
    line3.Draw("same");

    c_twoRef.SaveAs( (string("fullratio_PhyRepPRLCompare"  )  + string(".pdf")).c_str());
}


// Plot ratio in this analysis only
void Plot_Ratio_ThisAnalysisOnly(TGraphErrors *g_HighResult, TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, std::string issversion, double Scaler){

    TCanvas c1("", "", 1000, 500);

    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1->Draw();
    p1->cd();

    g_HighResult        ->Draw("AP");
    g_LowResult         ->Draw("same P");
    g_IntermediateResult->Draw("same P");

    g_HighResult        ->SetMarkerStyle(15);
    g_HighResult        ->SetMarkerColor(2);
    g_HighResult        ->SetMarkerSize(0.9);
    g_LowResult         ->SetMarkerStyle(15);
    g_LowResult         ->SetMarkerColor(2);
    g_LowResult         ->SetMarkerSize(0.9);
    g_IntermediateResult->SetMarkerStyle(15);
    g_IntermediateResult->SetMarkerColor(2);
    g_IntermediateResult->SetMarkerSize(0.9);

    /*
    for (int p = 0; p < g_LowResult->GetN(); p++) {
        cout<< "low here:" << g_LowResult->GetX()[p] << endl;
        cout<< "low here:" << g_LowResult->GetY()[p] << endl;
    }
    for (int p = 0; p < g_IntermediateResult->GetN(); p++) {
        cout<< "intermediate here:" << g_IntermediateResult->GetX()[p] << endl;
        cout<< "intermediate here:" << g_IntermediateResult->GetY()[p] << endl;
    }
    for (int p = 0; p < g_HighResult->GetN(); p++) {
        cout<< "high here:" << g_HighResult->GetX()[p] << endl;
        cout<< "high here:" << g_HighResult->GetY()[p] << endl;
    }
    */

    TAxis * xaxis = g_HighResult->GetXaxis();
    TAxis * yaxis = g_HighResult->GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#bar{p}/p");
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.055);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);

    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(1.3);

    p1->SetBottomMargin(0.2);
    p1->SetLeftMargin  (0.15);

    g_HighResult->SetTitle("");
    //yaxis->SetMoreLogLabels();
    xaxis->SetLabelOffset(0.02);
 
    // Plot
    xaxis->SetLimits   (1.0    , 600);
    yaxis->SetRangeUser(0.00001, 0.00028);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    xaxis->SetLabelOffset(0);
    c1.SaveAs( (string("fullratio_ThisAnallysisOnly_") + issversion + string("_LogXLogY.pdf")).c_str() );

    xaxis->SetLimits   (1.0        , 600);
    yaxis->SetRangeUser(0.000000001, 0.00028);
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    xaxis->SetMoreLogLabels();
    xaxis->SetLabelOffset(0);
    c1.SaveAs( (string("fullratio_ThisAnallysisOnly_") + issversion + string(".pdf")).c_str() );

    xaxis->SetLimits   (1.0    , 550);
    yaxis->SetRangeUser(0.00001, 0.00028);
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    // Straight line Fit
    int R_StartFit = 40;
    TF1 *linear = new TF1("linear", "[0]*log(x)+[1]", R_StartFit, 525);
    g_HighResult->Fit(linear, "", "", R_StartFit, 525);
    TF1 * f1 = g_HighResult->GetFunction("linear");
    f1->SetLineColor(kRed);
    f1->SetRange(R_StartFit, 525);
    f1->Draw("same");
    cout<< "Fit From:"        << R_StartFit << " GV"             <<endl;
    cout<< "Scaler:"          << Scaler                          <<endl;
    cout<< "Chi2:"            << f1->GetChisquare()              <<endl;
    cout<< "NDF:"             << f1->GetNDF()                    <<endl;
    cout<< "Chi2/NDF:"        << f1->GetChisquare()/f1->GetNDF() <<endl; 
    cout<< "Fit probability:" << f1->GetProb()                   <<endl; //Return the fit probability
    //Create a TGraphErrors to hold the confidence intervals//
    TGraphErrors *grint = new TGraphErrors( 525-R_StartFit );  //465 for from 60 GV, 485 for from 40 GV.   
    for (int i=R_StartFit; i<526; i=i+1){
        grint->SetPoint(i-R_StartFit, i, 0);
    }
    // Compute the confidence intervals at the x points of the created graph
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint, 0.68);
    //Now the "grint" graph contains function values as its y-coordinates and confidence intervals as the errors on these coordinates
    //Draw the graph, the function and the confidence intervals
    c1.cd(1);
    grint->SetLineColor(-9);
    grint->SetMarkerColor(-9);
    grint->SetFillColor(kRed-9);
    grint->SetLineWidth(0);
    grint->SetFillStyle(1001);
    grint->Draw("same P3");
    g_HighResult->Draw("same P");
    
    c1.SaveAs( (string("fullratio_ThisAnallysisOnly_") + issversion + string("_LogXLinearY.pdf")).c_str() );

    gPad->SetLogy(1);
    gPad->SetLogx(0);
    xaxis->SetLimits   (1.0    , 550);
    yaxis->SetRangeUser(0.00001, 0.00035);

    c1.SaveAs( (string("fullratio_ThisAnallysisOnly_") + issversion + string("_LinearXLogY.pdf")).c_str() );


}


void Plot_Ratio_OnlyStaErr(TGraphErrors *g_HighResult, TGraphErrors *g_LowResult, TGraphErrors *g_IntermediateResult, std::string issversion, int LowRangeRemoveAtEnd, int IntermediateRangeRemoveAtBegin, int IntermediateRangeRemoveAtEnd, int HighRangeRemovedAtBegin){

    TCanvas c1("", "", 1000, 500);

    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1->Draw();
    p1->cd();

    //// Remove the first two points in interrmediate range for illustration.
    TGraphErrors *g_IntermediateResult_FirstTwoRemoved  = new TGraphErrors(*g_IntermediateResult);
    g_IntermediateResult_FirstTwoRemoved->RemovePoint(0);
    g_IntermediateResult_FirstTwoRemoved->RemovePoint(0);

    g_HighResult                        ->Draw("AP");
    g_LowResult                         ->Draw("same P");
    g_IntermediateResult_FirstTwoRemoved->Draw("same P");

    g_HighResult->SetMarkerStyle(15);
    g_HighResult->SetMarkerColor(28);
    g_HighResult->SetLineColor(28);
    g_HighResult->SetMarkerSize(0.9);

    g_LowResult->SetMarkerStyle(15);
    g_LowResult->SetMarkerColor(6);
    g_LowResult->SetLineColor(6);
    g_LowResult->SetMarkerSize(0.9);

    g_IntermediateResult_FirstTwoRemoved->SetMarkerStyle(15);
    g_IntermediateResult_FirstTwoRemoved->SetMarkerColor(4);
    g_IntermediateResult_FirstTwoRemoved->SetLineColor(4);
    g_IntermediateResult_FirstTwoRemoved->SetMarkerSize(0.9);

    TAxis * xaxis = g_HighResult->GetXaxis();
    TAxis * yaxis = g_HighResult->GetYaxis();

    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#bar{p}/p");
    //yaxis->CenterTitle(true);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);

    xaxis->SetLabelOffset(-0.01);
    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(0.9);

    p1->SetBottomMargin(0.18);
    p1->SetLeftMargin  (0.12);

    g_HighResult->SetTitle("");

    TLatex latex;
    latex.SetNDC(1);
    latex.SetTextSize(0.04);
    latex.SetTextAngle(0);
    latex.SetTextAlign(13);  //align at top
    latex.DrawLatex(0.2, 0.8, "Statistical uncertainty only");

    TLegend *legend1 = new TLegend(0.42, 0.22, 0.82, 0.42);
    legend1->AddEntry(g_LowResult                         , "Low Rigidity Range"         , "p");
    legend1->AddEntry(g_IntermediateResult_FirstTwoRemoved, "Intermediate Rigidity Range", "p");
    legend1->AddEntry(g_HighResult                        , "High Rigidity Range"        , "p");
    legend1->SetTextSize(0.04);
    legend1->SetTextFont(62);
    legend1->SetBorderSize(0);
    legend1->Draw();

    xaxis->SetLimits   (1.0    , 550);
    yaxis->SetRangeUser(0.00001, 0.00035);
    gPad->SetLogy(1);
    gPad->SetLogx(0);
    c1.SaveAs( (string("fullratio_StaErrOnly_ThisAnallysisOnly_") + issversion + string("_LinearXLogY.pdf")).c_str() );

    xaxis->SetLimits   (1.0    , 600);
    yaxis->SetRangeUser(0.00001, 0.00028);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    c1.SaveAs( (string("fullratio_StaErrOnly_ThisAnallysisOnly_") + issversion + string("_LogXLogY.pdf")).c_str() );

    xaxis->SetLimits   (1.0        , 600);
    yaxis->SetRangeUser(0.000000001, 0.00028);
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    xaxis->SetMoreLogLabels();
    c1.SaveAs( (string("fullratio_StaErrOnly_ThisAnallysisOnly_") + issversion + string(".pdf")).c_str() );


    // Not overlapped
    TGraphErrors *g_LowResult_NotOverLap           = new TGraphErrors(*g_LowResult);
    TGraphErrors *g_IntermediateResult_NotOverLap  = new TGraphErrors(*g_IntermediateResult);
    TGraphErrors *g_HighResult_NotOverLap          = new TGraphErrors(*g_HighResult);

    for (int i=0; i<LowRangeRemoveAtEnd; i++){
        g_LowResult_NotOverLap->RemovePoint(g_LowResult_NotOverLap->GetN()-1);
    }
    for (int i=0; i<IntermediateRangeRemoveAtBegin; i++){
        g_IntermediateResult_NotOverLap->RemovePoint(0);
    }
    for (int i=0; i<IntermediateRangeRemoveAtEnd; i++){
        g_IntermediateResult_NotOverLap->RemovePoint(g_IntermediateResult_NotOverLap->GetN()-1);
    }
    for (int i=0; i<HighRangeRemovedAtBegin; i++){
        g_HighResult_NotOverLap->RemovePoint(0);
    }

    g_HighResult_NotOverLap        ->Draw("AP");
    g_LowResult_NotOverLap         ->Draw("same P");
    g_IntermediateResult_NotOverLap->Draw("same P");
    
    g_HighResult_NotOverLap        ->SetMarkerColor(2);
    g_HighResult_NotOverLap        ->SetLineColor(2);
    g_HighResult_NotOverLap        ->SetMarkerSize(0.9);
    g_LowResult_NotOverLap         ->SetMarkerColor(2);
    g_LowResult_NotOverLap         ->SetLineColor(2);
    g_LowResult_NotOverLap         ->SetMarkerSize(0.9);
    g_IntermediateResult_NotOverLap->SetMarkerColor(2);
    g_IntermediateResult_NotOverLap->SetLineColor(2);   
    g_IntermediateResult_NotOverLap->SetMarkerSize(0.9); 

    TLatex latex2;
    latex2.SetNDC(1);
    latex2.SetTextSize(0.04);
    latex2.SetTextAngle(0);
    latex2.SetTextAlign(13);  //align at top
    latex2.DrawLatex(0.2, 0.8, "Statistical uncertainty only");

    xaxis->SetLimits   (1.0    , 550);
    yaxis->SetRangeUser(0.00001, 0.00035);
    gPad->SetLogy(1);
    gPad->SetLogx(0);
    c1.SaveAs( (string("fullratio_StaErrOnly_NotOverlapped_ThisAnallysisOnly_") + issversion + string("_LinearXLogY.pdf")).c_str() );

    xaxis->SetLimits   (1.0    , 600);
    yaxis->SetRangeUser(0.00001, 0.00028);
    gPad->SetLogy(1);
    gPad->SetLogx(1);
    c1.SaveAs( (string("fullratio_StaErrOnly_NotOverlapped_ThisAnallysisOnly_") + issversion + string("_LogXLogY.pdf")).c_str() );

    xaxis->SetLimits   (1.0        , 600);
    yaxis->SetRangeUser(0.000000001, 0.00028);
    gPad->SetLogy(0);
    gPad->SetLogx(1);
    xaxis->SetMoreLogLabels();
    c1.SaveAs( (string("fullratio_StaErrOnly_NotOverlapped_ThisAnallysisOnly_") + issversion + string(".pdf")).c_str() );

}


// Plot Statictical Error (absolute value)
void Plot_Statictical_Error_Absolute(TGraph *g_PublishedRatioStatisticErrorPRL, TGraphErrors *g_StatisticalError_Low, TGraph *g_StatisticError_Intermediate, TGraph *g_StatisticalError_High, TGraph *g_PhysicsReportRatioStatisticError, std::string issversionname){

    TCanvas c_StaErr("c_StaErr","c_StaErr",1000,500);
    TPad *pad_StaErr = new TPad("pad_StaErr","pad_StaErr",0.03,0.62,0.50,0.92,32);

    g_PublishedRatioStatisticErrorPRL->SetMarkerStyle(15);
    g_PublishedRatioStatisticErrorPRL->SetMarkerColor(0);
    g_PublishedRatioStatisticErrorPRL->SetMarkerSize(0.9);
    g_StatisticalError_Low->SetMarkerStyle(15);
    g_StatisticalError_Low->SetMarkerColor(2);
    g_StatisticalError_Low->SetMarkerSize(0.9);
    g_StatisticError_Intermediate->SetMarkerStyle(15);
    g_StatisticError_Intermediate->SetMarkerColor(2);
    g_StatisticError_Intermediate->SetMarkerSize(0.9);
    g_StatisticalError_High->SetMarkerStyle(15);
    g_StatisticalError_High->SetMarkerColor(2);
    g_StatisticalError_High->SetMarkerSize(0.9);
    g_PhysicsReportRatioStatisticError->SetMarkerStyle(15);
    g_PhysicsReportRatioStatisticError->SetMarkerColor(4);
    g_PhysicsReportRatioStatisticError->SetMarkerSize(0.9);

    gPad->SetLogx();
    //gPad->SetLogy();

    TAxis * xaxis_StaErr = g_PublishedRatioStatisticErrorPRL->GetXaxis();
    TAxis * yaxis_StaErr = g_PublishedRatioStatisticErrorPRL->GetYaxis();
    xaxis_StaErr->SetTitle("|R| / (GV)");
    yaxis_StaErr->SetTitle("Absolute statistical uncertainty");
    xaxis_StaErr->SetLimits(1.0,600);
    yaxis_StaErr->SetRangeUser(0, 0.00008);
    xaxis_StaErr->SetTitleFont(62);
    yaxis_StaErr->SetTitleFont(62);
    xaxis_StaErr->SetTitleSize(0.045);
    yaxis_StaErr->SetTitleSize(0.045);
    xaxis_StaErr->SetLabelFont(62);
    yaxis_StaErr->SetLabelFont(62);
    xaxis_StaErr->SetLabelSize(0.05);
    yaxis_StaErr->SetLabelSize(0.05);
    xaxis_StaErr->SetTitleOffset(0.8);
    yaxis_StaErr->SetTitleOffset(1.1);

    g_PublishedRatioStatisticErrorPRL ->Draw("AP");
    g_StatisticalError_Low            ->Draw("same P");
    g_StatisticError_Intermediate     ->Draw("same P");
    g_StatisticalError_High           ->Draw("same P");
    //g_PhysicsReportRatioStatisticError->Draw("same P");

    cout<< "StatisticError:" << endl;
    for (int p = 0; p < g_StatisticalError_Low->GetN(); p++) {
        cout<< "low here:" << g_StatisticalError_Low->GetX()[p] << endl;
        cout<< "low here:" << g_StatisticalError_Low->GetY()[p] << endl;
    }    
    

    for (int p = 0; p < g_StatisticError_Intermediate->GetN(); p++) {
        cout<< "intermediate here:" << g_StatisticError_Intermediate->GetX()[p] << endl;
        cout<< "intermediate here:" << g_StatisticError_Intermediate->GetY()[p] << endl;
    }

    for (int p = 0; p < g_StatisticalError_High->GetN(); p++) {
        cout<< "high here:" << g_StatisticalError_High->GetX()[p] << endl;
        cout<< "high here:" << g_StatisticalError_High->GetY()[p] << endl;
    }
    

    g_PublishedRatioStatisticErrorPRL->SetTitle("");

    //TLegend *legend_staerr_absolute = new TLegend(0.35,0.7,0.65,0.85);
    //legend_staerr_absolute->AddEntry(g_StatisticalError_Low, (string("This analysis")).c_str(),"lpf");
    //legend_staerr_absolute->AddEntry(g_PublishedRatioStatisticErrorPRL,"PRL paper 2016","lpf");
    //legend_staerr_absolute->AddEntry(g_PhysicsReportRatioStatisticError,"Physics Report 2021","lpf");
    //legend_staerr_absolute->SetTextSize(0.05);
    //legend_staerr_absolute->SetTextFont(62);
    //legend_staerr_absolute->Draw();

    c_StaErr.SaveAs( (string("StaErr_") + issversionname + string("_LinearY.pdf")).c_str());
    gPad->SetLogy();
    c_StaErr.SaveAs( (string("StaErr_") + issversionname + string("_LogY.pdf")).c_str());

}


// Plot Statictical Error (relative value)
void Plot_Statictical_Error_Relative(TH1D h_PhysicsReportStatisticRatioRelativeError, TH1D h_PublishedRatioStatisticRelativeErrorPRL, TH1D *h_StatisticalRelError_Low, TH1D *h_StatisticalRelError_Intermediate, TH1D *h_Statistic_Error_Relative_High, std::vector<double> PublishedPRLBinEdge, std::string issversion, std::string issversionname){

    TCanvas c_StaRelErr("c_StaRelErr","c_StaRelErr",1000,500);
    TPad *pad_StaRelErr = new TPad("pad_StaRelErr","pad_StaRelErr", 0.0, 0.0, 1.0, 1.0);

    pad_StaRelErr->Draw("");
    pad_StaRelErr->cd();

    h_PhysicsReportStatisticRatioRelativeError.Draw("same HIST ][");
    h_StatisticalRelError_Low                ->Draw("same HIST ][");
    h_StatisticalRelError_Intermediate       ->Draw("same HIST ][");
    h_Statistic_Error_Relative_High          ->Draw("same HIST ][");
    //h_PhysicsReportStatisticRatioRelativeError.Draw("same HIST ][");
    //h_PublishedRatioStatisticRelativeErrorPRL .Draw("same HIST ][");

    /*
    TLine *l_PRLRelativeStaErrorPRL = new TLine(PublishedPRLBinEdge.back(), 0, PublishedPRLBinEdge.back(), h_PublishedRatioStatisticRelativeErrorPRL.GetBinContent(h_PublishedRatioStatisticRelativeErrorPRL.GetNbinsX()));
    l_PRLRelativeStaErrorPRL->SetLineWidth(3);
    l_PRLRelativeStaErrorPRL->Draw();

    if (issversion == "2016paper"){
        TLine *l_2016HighRelativeStaErrorPRL = new TLine(PublishedPRLBinEdge.back(), 0, PublishedPRLBinEdge.back(), h_Statistic_Error_Relative_High->GetBinContent(h_Statistic_Error_Relative_High->GetNbinsX()));
        l_2016HighRelativeStaErrorPRL->SetLineWidth(3);
        l_2016HighRelativeStaErrorPRL->SetLineColor(2);
        l_2016HighRelativeStaErrorPRL->Draw();
    }
    */

    h_PhysicsReportStatisticRatioRelativeError.SetStats(0);

    h_PublishedRatioStatisticRelativeErrorPRL.SetLineColor(0); //1
    h_PhysicsReportStatisticRatioRelativeError.SetLineColor(0); //4
    h_PublishedRatioStatisticRelativeErrorPRL.SetLineWidth(0); //3
    h_PhysicsReportStatisticRatioRelativeError.SetLineWidth(0); //3

    h_StatisticalRelError_Low         ->SetLineColor(2);
    h_StatisticalRelError_Intermediate->SetLineColor(2);
    h_Statistic_Error_Relative_High   ->SetLineColor(2);
    h_StatisticalRelError_Low         ->SetLineWidth(3);
    h_StatisticalRelError_Intermediate->SetLineWidth(3);
    h_Statistic_Error_Relative_High   ->SetLineWidth(3);

    TAxis * xaxis_StaRelErr = h_PhysicsReportStatisticRatioRelativeError.GetXaxis();
    TAxis * yaxis_StaRelErr = h_PhysicsReportStatisticRatioRelativeError.GetYaxis();
    xaxis_StaRelErr->SetTitle("|R| / (GV)");
    yaxis_StaRelErr->SetTitle("Relative statistical uncertainty (%)");
    xaxis_StaRelErr->SetLimits(1, 600);
    yaxis_StaRelErr->SetRangeUser(0, 35);
    xaxis_StaRelErr->SetTitleFont(62);
    yaxis_StaRelErr->SetTitleFont(62);
    xaxis_StaRelErr->SetTitleSize(0.045);
    yaxis_StaRelErr->SetTitleSize(0.045);
    xaxis_StaRelErr->SetLabelFont(62);
    yaxis_StaRelErr->SetLabelFont(62);
    xaxis_StaRelErr->SetLabelSize(0.05);
    yaxis_StaRelErr->SetLabelSize(0.05);
    xaxis_StaRelErr->SetMoreLogLabels();

    xaxis_StaRelErr->SetTitleOffset(1.3);
    yaxis_StaRelErr->SetTitleOffset(1.0);
    pad_StaRelErr->SetBottomMargin(0.13);
    pad_StaRelErr->SetLeftMargin(0.15);

    /*
    TLegend *legend_staerr_relative = new TLegend(0.2,0.7,0.75,0.85);
    if (issversion == "pass7.8"){
        legend_staerr_relative->AddEntry(h_StatisticalRelError_Low, (string("This analysis (All time range)")).c_str(),"lpf");
    }
    else if (issversion == "2016paper"){
        legend_staerr_relative->AddEntry(h_StatisticalRelError_Low, (string("This analysis (PRL paper time range)")).c_str(),"lpf");
    }
    else if (issversion == "PhyRep2021"){
        legend_staerr_relative->AddEntry(h_StatisticalRelError_Low, (string("This analysis (Phy Rep time range)")).c_str(),"lpf");
    }
    legend_staerr_relative->AddEntry(&h_PublishedRatioStatisticRelativeErrorPRL,"PRL paper 2016","lpf");
    legend_staerr_relative->AddEntry(&h_PhysicsReportStatisticRatioRelativeError,"Physics Report 2021","lpf");
    legend_staerr_relative->SetTextSize(0.05);
    legend_staerr_relative->SetTextFont(62);
    legend_staerr_relative->Draw();
    */

    gPad->SetLogx();
    c_StaRelErr.SaveAs( (string("StaRelErr_") + issversionname + string(".pdf")).c_str());
}


// Plot Absolute Systematic error
void Plot_Systematical_Error_Absolute_Acc(TGraph *g_System_ACC_Low, TGraph *g_System_ACC_Intermediate, TGraph *g_System_ACC_High){
    TCanvas c1("c1","c1",1000,500);

    g_System_ACC_High         ->Draw("AP");
    g_System_ACC_Intermediate ->Draw("same P");
    g_System_ACC_Low          ->Draw("same P");

    g_System_ACC_High->SetMarkerStyle(15);
    g_System_ACC_High->SetMarkerColor(8); 
    g_System_ACC_High->SetMarkerSize(0.9);
    g_System_ACC_Intermediate->SetMarkerStyle(15);
    g_System_ACC_Intermediate->SetMarkerColor(8);
    g_System_ACC_Intermediate->SetMarkerSize(0.8); 
    g_System_ACC_Low->SetMarkerStyle(15);
    g_System_ACC_Low->SetMarkerColor(8); 
    g_System_ACC_Low->SetMarkerSize(0.9);

    g_System_ACC_High->SetTitle("");

    TAxis * xaxis = g_System_ACC_High->GetXaxis();
    TAxis * yaxis = g_System_ACC_High->GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#splitline{Absolute systematic uncertainty}{           (From Acceptance)}");
    xaxis->SetLimits(1.0, 600);  // 1.0, 600
    yaxis->SetRangeUser(0,0.00010); //0,0.00008
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetTitleOffset(0.8);
    yaxis->SetTitleOffset(1.1);

    xaxis->SetTitleOffset(1.3);
    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.13);

    gPad->SetLogx();
    c1.SaveAs( (string("SysErr_Acceptance") + string(".pdf")).c_str());

    gPad->SetLogy();
    c1.SaveAs( (string("SysErr_Acceptance_LogY") + string(".pdf")).c_str());
}


void Plot_Systematical_Error_Absolute_FitRange(TGraphErrors *g_SystematicError_TRD_Intermediate, TGraphErrors *g_SystematicError_Shape_Low, std::string issversionname){

    TCanvas c1("","",1000,500);

    g_SystematicError_TRD_Intermediate->Draw("AP");
    g_SystematicError_Shape_Low       ->Draw("same P");

    g_SystematicError_TRD_Intermediate->SetMarkerStyle(15);
    g_SystematicError_TRD_Intermediate->SetMarkerColor(4);
    g_SystematicError_TRD_Intermediate->SetMarkerSize(0.9); 
    g_SystematicError_Shape_Low       ->SetMarkerStyle(15);
    g_SystematicError_Shape_Low       ->SetMarkerColor(6);  
    g_SystematicError_Shape_Low       ->SetMarkerSize(0.9);

    g_SystematicError_TRD_Intermediate->SetTitle("");

    TAxis * xaxis = g_SystematicError_TRD_Intermediate->GetXaxis();
    TAxis * yaxis = g_SystematicError_TRD_Intermediate->GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#splitline{Absolute systematic uncertainty}{           (From Fit Range)}");
    xaxis->SetLimits(1.0, 20); 
    yaxis->SetRangeUser(0,0.00008); 
    xaxis->SetMoreLogLabels();
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);

    xaxis->SetTitleOffset(1.3);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.16);
    yaxis->SetTitleOffset(1.1);

    TLegend *legend = new TLegend(0.25, 0.75, 0.45, 0.85);  // 0.25, 0.7, 0.65, 0.9
    legend->AddEntry(g_SystematicError_TRD_Intermediate , "Intermediate rigidity range", "p");
    legend->AddEntry(g_SystematicError_Shape_Low        , "Low rigidity range"         , "p");
    legend->SetTextSize(0.04);
    legend->SetTextFont(62);
    legend->SetBorderSize(0);
    legend->Draw();


    gPad->SetLogx();
    c1.SaveAs( (string("SysErr_FitRange_") + issversionname + string(".pdf")).c_str());

    gPad->SetLogy();
    c1.SaveAs( (string("SysErr_FitRange_LogY_") + issversionname + string(".pdf")).c_str());

}


void Plot_Systematical_Error_Absolute_CC(TGraph * g_System_CC_High, std::string issversionname){

    TCanvas c1("","",1000,500);

    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0);
    p1->Draw("");
    p1->cd();

    g_System_CC_High->Draw("AP");

    g_System_CC_High->SetMarkerStyle(15);
    g_System_CC_High->SetMarkerColor(6);
    g_System_CC_High->SetMarkerSize(0.9);

    g_System_CC_High->SetTitle("");

    TAxis * xaxis = g_System_CC_High->GetXaxis();
    TAxis * yaxis = g_System_CC_High->GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#splitline{Absolute systematic uncertainty}{      (From Charge Confusion)}");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0,0.00008);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetMoreLogLabels();

    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(1.3);
    p1->SetBottomMargin(0.13);
    p1->SetLeftMargin(0.15);

    gPad->SetLogx();
    c1.SaveAs( (string("SysErr_CC_") + issversionname + string(".pdf")).c_str());

    gPad->SetLogy();
    c1.SaveAs( (string("SysErr_CC_LogY_") + issversionname + string(".pdf")).c_str());

}


void Plot_Systematical_Error_Relative_CC( TH1D *h_System_CC_relative_High, std::string issversionname){

    TCanvas c1("","",1000,500);

    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0);
    p1->Draw("");
    p1->cd();

    h_System_CC_relative_High->Draw("");

    h_System_CC_relative_High->SetTitle("");

    h_System_CC_relative_High->SetLineColor(6);
    h_System_CC_relative_High->SetLineWidth(3);

    TAxis * xaxis = h_System_CC_relative_High->GetXaxis();
    TAxis * yaxis = h_System_CC_relative_High->GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#splitline{Relative systematic uncertainty (%)}{      (From Charge Confusion)}");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0, 50);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetMoreLogLabels();

    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(1.3);
    p1->SetBottomMargin(0.13);
    p1->SetLeftMargin(0.15);

    c1.SaveAs( (string("SysRelErr_CC_") + issversionname + string(".pdf")).c_str());

    gPad->SetLogx();
    c1.SaveAs( (string("SysRelErr_CC_") + issversionname + string("_LogX.pdf")).c_str());

}


void Plot_Systematical_Error_Relative_ACC(TH1D h_SysUncertaintyRel_ACCratio, std::string issversionname){

    TCanvas c1("","",1000,500);

    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0);
    p1->Draw("");
    p1->cd();

    h_SysUncertaintyRel_ACCratio.Draw("");

    h_SysUncertaintyRel_ACCratio.SetTitle("");

    h_SysUncertaintyRel_ACCratio.SetLineColor(8);
    h_SysUncertaintyRel_ACCratio.SetLineWidth(3);

    TAxis * xaxis = h_SysUncertaintyRel_ACCratio.GetXaxis();
    TAxis * yaxis = h_SysUncertaintyRel_ACCratio.GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#splitline{Relative systematic uncertainty (%)}{           (From Acceptance)}");
    xaxis->SetLimits(10, 600);
    yaxis->SetRangeUser(0, 10);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetMoreLogLabels();

    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(1.3);
    p1->SetBottomMargin(0.13);
    p1->SetLeftMargin(0.15);

    gPad->SetLogx();
    c1.SaveAs( (string("SysRelErr_ACC_") + issversionname + string(".pdf")).c_str());

}


void Plot_Systematical_Error_Relative_FitRange( TH1D *h_SystematicError_Shape_Low, TH1D *h_SystematicRelativeError_TRD, TGraph *g_SystematicError_Shape_Low, std::string issversionname){

    TCanvas c1("","",1000,500);

    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0);
    p1->Draw("");
    p1->cd();

    g_SystematicError_Shape_Low  ->Draw("");
    h_SystematicError_Shape_Low  ->Draw("same HIST ][");
    h_SystematicRelativeError_TRD->Draw("same HIST ][");

    h_SystematicError_Shape_Low->SetTitle("");

    g_SystematicError_Shape_Low->SetLineColor(0);
    g_SystematicError_Shape_Low->SetMarkerSize(0);
    h_SystematicError_Shape_Low->SetLineColor(45);
    h_SystematicError_Shape_Low->SetLineWidth(3);
    h_SystematicRelativeError_TRD      ->SetLineColor(45);
    h_SystematicRelativeError_TRD      ->SetLineWidth(3);

    TAxis * xaxis = g_SystematicError_Shape_Low->GetXaxis();
    TAxis * yaxis = g_SystematicError_Shape_Low->GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("#splitline{Relative systematic uncertainty (%)}{           (From Fit Range)}");
    xaxis->SetLimits(1, 18);
    yaxis->SetRangeUser(0, 15);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    //xaxis->SetMoreLogLabels();

    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(1.3);
    p1->SetBottomMargin(0.13);
    p1->SetLeftMargin(0.15);

    gPad->Update();
    c1.Update();
    //gPad->SetLogx();
    c1.SaveAs( (string("SysRelErr_FitRange_") + issversionname + string(".pdf")).c_str());

}



void Plot_Systematical_Error_Absolute( TGraph *g_TotalSysError_High, TGraph *g_System_CC_High, TGraph *g_TotalSysError_Intermediate, TGraphErrors *g_SystematicError_TRD_Intermediate, TGraph *g_TotalSysError_Low, TGraphErrors *g_SystematicError_Shape_Low, TGraph *g_System_ACC_Low, TGraph *g_System_ACC_Intermediate, TGraph *g_System_ACC_High, TGraph *g_PublishedRatioSystematicErrorPRL, TGraph *g_PhysicsReportRatioSystematicError, std::string issversionname ){
    TCanvas c_SysAboErr("c_SysAboErr","c_SysAboErr",1000,500);

    // High
    g_System_CC_High      ->Draw("AP");
    //g_System_ACC_High   ->Draw("same P");
    g_TotalSysError_High->Draw("same P");
    // Intermediate
    //g_System_ACC_Intermediate         ->Draw("same P");
    //g_SystematicError_TRD_Intermediate->Draw("same P");
    g_TotalSysError_Intermediate      ->Draw("same P");
    // Low
    //g_System_ACC_Low           ->Draw("same P");
    //g_SystematicError_Shape_Low->Draw("same P");
    g_TotalSysError_Low        ->Draw("same P");
    // Reference
    //g_PublishedRatioSystematicErrorPRL ->Draw("same P");
    //g_PhysicsReportRatioSystematicError->Draw("same P");

   
    cout<< "Systematical Error:" << endl; 
    for (int p = 0; p < g_TotalSysError_Low->GetN(); p++) {
        cout<< "low here:" << g_TotalSysError_Low->GetX()[p] << endl;
        cout<< "low here:" << g_TotalSysError_Low->GetY()[p] << endl;
    }

    for (int p = 0; p < g_TotalSysError_Intermediate->GetN(); p++) {
        cout<< "Intermediate here:" << g_TotalSysError_Intermediate->GetX()[p] << endl;
        cout<< "Intermediate here:" << g_TotalSysError_Intermediate->GetY()[p] << endl;
    }

    for (int p = 0; p < g_TotalSysError_High->GetN(); p++) {
        cout<< "High here:" << g_TotalSysError_High->GetX()[p] << endl;
        cout<< "High here:" << g_TotalSysError_High->GetY()[p] << endl;
    }
    


    // Reference
    //g_PublishedRatioSystematicErrorPRL->SetMarkerStyle(15);
    //g_PublishedRatioSystematicErrorPRL->SetMarkerColor(1);
    g_PhysicsReportRatioSystematicError->SetMarkerStyle(15);
    g_PhysicsReportRatioSystematicError->SetMarkerColor(0); //4
    // High
    g_TotalSysError_High->SetMarkerStyle(15);
    g_TotalSysError_High->SetMarkerColor(4); 
    g_TotalSysError_High->SetMarkerSize(0.9);
    g_System_CC_High->SetMarkerStyle(15);
    g_System_CC_High->SetMarkerColor(0);//6 
    g_System_CC_High->SetMarkerSize(0.9);
    g_System_ACC_High->SetMarkerStyle(15);
    g_System_ACC_High->SetMarkerColor(0); //3
    g_System_ACC_High->SetMarkerSize(0.9);
    // Intermediate
    g_System_ACC_Intermediate->SetMarkerStyle(15);
    g_System_ACC_Intermediate->SetMarkerColor(0); //3
    g_System_ACC_Intermediate->SetMarkerSize(0.9);
    g_SystematicError_TRD_Intermediate->SetMarkerStyle(15);
    g_SystematicError_TRD_Intermediate->SetMarkerColor(0); //46
    g_SystematicError_TRD_Intermediate->SetMarkerSize(0.9);
    g_TotalSysError_Intermediate->SetMarkerStyle(15);
    g_TotalSysError_Intermediate->SetMarkerColor(4); 
    g_TotalSysError_Intermediate->SetMarkerSize(0.9);
    // Low
    g_System_ACC_Low->SetMarkerStyle(15);
    g_System_ACC_Low->SetMarkerColor(0); //3
    g_System_ACC_Low->SetMarkerSize(0.9);
    g_SystematicError_Shape_Low->SetMarkerStyle(15);
    g_SystematicError_Shape_Low->SetMarkerColor(0);  //7
    g_SystematicError_Shape_Low->SetMarkerSize(0.9);
    g_TotalSysError_Low->SetMarkerStyle(15);
    g_TotalSysError_Low->SetMarkerColor(4); 
    g_TotalSysError_Low->SetMarkerSize(0.9);

    g_System_CC_High->SetTitle("");

    TAxis * xaxis_AboSysErr = g_System_CC_High->GetXaxis();
    TAxis * yaxis_AboSysErr = g_System_CC_High->GetYaxis();
    xaxis_AboSysErr->SetTitle("|R| / (GV)");
    yaxis_AboSysErr->SetTitle("Absolute total systematic uncertainty"); 
    xaxis_AboSysErr->SetLimits(1.0, 600);  // 1.0, 600
    yaxis_AboSysErr->SetRangeUser(0,0.00010); //0,0.00008
    xaxis_AboSysErr->SetTitleFont(62);
    yaxis_AboSysErr->SetTitleFont(62);
    xaxis_AboSysErr->SetTitleSize(0.045);
    yaxis_AboSysErr->SetTitleSize(0.045);
    xaxis_AboSysErr->SetLabelFont(62);
    yaxis_AboSysErr->SetLabelFont(62);
    xaxis_AboSysErr->SetLabelSize(0.05);
    yaxis_AboSysErr->SetLabelSize(0.05);

    gPad->SetBottomMargin(0.2);
    xaxis_AboSysErr->SetTitleOffset(1.3);
    yaxis_AboSysErr->SetTitleOffset(1.1);

    TLegend *legend_AboSysErr = new TLegend(0.25, 0.75, 0.45, 0.85);  // 0.25, 0.7, 0.65, 0.9
    //legend_AboSysErr->AddEntry(g_PublishedRatioSystematicErrorPRL , "PRL paper 2016","lpf");
    //legend_AboSysErr->AddEntry(g_TotalSysError_High               , "This analysis(Total)","lpf");

    //legend_AboSysErr->AddEntry(g_System_CC_High                   , "CC","lpf");
    //legend_AboSysErr->AddEntry(g_System_ACC_Low                   , "Acceptance","lpf");
    //legend_AboSysErr->AddEntry(g_SystematicError_TRD_Intermediate , "Fit in Intermediate range","lpf");
    //legend_AboSysErr->AddEntry(g_SystematicError_Shape_Low        , "Fit in Low range","lpf");
    //legend_AboSysErr->AddEntry(g_PhysicsReportRatioSystematicError, "Totol Systemcatic Error in Physics Report 2021","lpf");
    legend_AboSysErr->SetTextSize(0.04);
    legend_AboSysErr->SetTextFont(62);
    legend_AboSysErr->SetBorderSize(0);
    //legend_AboSysErr->Draw();

    gPad->SetLogx();
    c_SysAboErr.SaveAs( (string("SysErr_Tot_") + issversionname + string(".pdf")).c_str());

    gPad->SetLogy();
    c_SysAboErr.SaveAs( (string("SysErr_Tot_LogY_") + issversionname + string(".pdf")).c_str());
}






// Plot Relatvie systematical error
void Plot_Relatvie_Systematical_Error(TH1D h_SysUncertaintyRel_ACCratio, TH1D h_PublishedRatioSystematicRelativeErrorPRL, TH1D h_PhysicsReportSystematicRatioRelativeError, TH1D *h_System_CC_relative_High, TH1D *h_TotalSysRelError_High, TH1D *h_TotalSysRelError_Intermediate, TH1D *h_TotalSysRelError_Low, std::vector<double>PublishedPRLBinEdge, std::string issversionname){

    TCanvas c_SysRelErr("c_SysRelErr","c_SysRelErr",1000,500);
    TPad *p1 = new TPad("", "", 0.0, 0.0, 1.0, 1.0);
    p1->Draw("");
    p1->cd();

    // Referernce
    h_PhysicsReportSystematicRatioRelativeError.Draw("HIST");
    //h_PublishedRatioSystematicRelativeErrorPRL .Draw("same HIST ][");
    // This analysis
    h_TotalSysRelError_High        ->Draw("same HIST ][");
    h_TotalSysRelError_Intermediate->Draw("same HIST ][");
    h_TotalSysRelError_Low         ->Draw("same HIST ][");

    // Font Style 
    h_PublishedRatioSystematicRelativeErrorPRL .SetLineColor(0);
    h_PhysicsReportSystematicRatioRelativeError.SetLineColor(0);
    h_TotalSysRelError_High        ->SetLineColor(4);
    h_TotalSysRelError_Intermediate->SetLineColor(4);
    h_TotalSysRelError_Low         ->SetLineColor(4);

    h_PublishedRatioSystematicRelativeErrorPRL .SetLineWidth(0);
    h_PhysicsReportSystematicRatioRelativeError.SetLineWidth(0);
    h_TotalSysRelError_High        ->SetLineWidth(3);
    h_TotalSysRelError_Intermediate->SetLineWidth(3);
    h_TotalSysRelError_Low         ->SetLineWidth(3);

    /*
    // Line to fix 450-525 GV difference
    TLine *l_PRLRelativeSysErrorPRL = new TLine(PublishedPRLBinEdge.back(), 0, PublishedPRLBinEdge.back(), h_PublishedRatioSystematicRelativeErrorPRL.GetBinContent(h_PublishedRatioSystematicRelativeErrorPRL.GetNbinsX()));
    l_PRLRelativeSysErrorPRL->SetLineWidth(3);
    l_PRLRelativeSysErrorPRL->Draw();
    */

    // Axis
    TAxis * xaxis = h_PhysicsReportSystematicRatioRelativeError.GetXaxis();
    TAxis * yaxis = h_PhysicsReportSystematicRatioRelativeError.GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("Relative total systematic uncertainty (%)");
    xaxis->SetLimits(1, 600);
    yaxis->SetRangeUser(0, 50);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetMoreLogLabels();

    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(1.0);
    p1->SetBottomMargin(0.13);
    p1->SetLeftMargin(0.15);

    /*
    // Legend
    TLegend *legend_SysErr_relative = new TLegend(0.25,0.65,0.55,0.85);
    //legend_SysErr_relative->AddEntry(&h_SysUncertaintyRel_ACCratio  , (string("This analysis(Acc)")).c_str(), "lpf");
    //legend_SysErr_relative->AddEntry(h_System_CC_relative_Highi     , (string("This analysis(CC)")).c_str() , "lpf");
    legend_SysErr_relative->AddEntry(h_TotalSysRelError_High          , (string("This analysis")).c_str()     , "lpf");
    //legend_SysErr_relative->AddEntry(h_TotalSysRelError_Intermediate, (string("This analysis(Intermediate,)")).c_str(), "lpf"); 
    //legend_SysErr_relative->AddEntry(h_TotalSysRelError_Low         , (string("This analysis(Low)")).c_str(), "lpf");
    legend_SysErr_relative->AddEntry(&h_PublishedRatioSystematicRelativeErrorPRL , "PRL paper 2016"     , "lpf");
    legend_SysErr_relative->AddEntry(&h_PhysicsReportSystematicRatioRelativeError, "Physics Report 2021", "lpf");
    legend_SysErr_relative->SetTextSize(0.05);
    legend_SysErr_relative->SetTextFont(62);
    legend_SysErr_relative->Draw();
    */

    gPad->SetLogx();
    c_SysRelErr.SaveAs( (string("SysRelErr_") + issversionname + string(".pdf")).c_str());
}



void Plot_StaSysRelErrCompare(TH1D h_PublishedRatioSystematicRelativeErrorPRL, TH1D h_PhysicsReportSystematicRatioRelativeError, TH1D *h_TotalSysRelError_High, TH1D *h_TotalSysRelError_Intermediate, TH1D *h_TotalSysRelError_Low, TH1D h_PhysicsReportStatisticRatioRelativeError, TH1D h_PublishedRatioStatisticRelativeErrorPRL, TH1D *h_StatisticalRelError_Low, TH1D *h_StatisticalRelError_Intermediate, TH1D *h_Statistic_Error_Relative_High, std::vector<double>PublishedPRLBinEdge, std::string issversionname, std::string issversion){

    TCanvas c_StaSysRelErrCompare("c_StaSysRelErrCompare", "c_StaSysRelErrCompare", 1000, 500);
    TPad *p1 = new TPad("","", 0, 0, 1.0, 1.0, 0);
    p1->Draw();
    p1->cd();

    // Plot
    h_PhysicsReportSystematicRatioRelativeError.Draw("HIST ][");
    h_TotalSysRelError_High        ->Draw("same HIST ][");
    h_TotalSysRelError_Intermediate->Draw("same HIST ][");
    h_TotalSysRelError_Low         ->Draw("same HIST ][");

    h_PhysicsReportStatisticRatioRelativeError.Draw("same HIST ][");
    h_StatisticalRelError_Low         ->Draw("same HIST ][");
    h_StatisticalRelError_Intermediate->Draw("same HIST ][");
    h_Statistic_Error_Relative_High   ->Draw("same HIST ][");

    // Font Style
    h_PhysicsReportSystematicRatioRelativeError.SetLineColor(28);
    h_TotalSysRelError_High        ->SetLineColor(4);
    h_TotalSysRelError_Intermediate->SetLineColor(4);
    h_TotalSysRelError_Low         ->SetLineColor(4);
    h_PhysicsReportStatisticRatioRelativeError.SetLineColor(6);
    h_StatisticalRelError_Low         ->SetLineColor(2);
    h_StatisticalRelError_Intermediate->SetLineColor(2);
    h_Statistic_Error_Relative_High   ->SetLineColor(2);

    h_PhysicsReportSystematicRatioRelativeError.SetLineWidth(3);
    h_PhysicsReportStatisticRatioRelativeError .SetLineWidth(3);
    h_TotalSysRelError_High        ->SetLineWidth(3);
    h_TotalSysRelError_Intermediate->SetLineWidth(3);
    h_TotalSysRelError_Low         ->SetLineWidth(3);
    h_StatisticalRelError_Low         ->SetLineWidth(3);
    h_StatisticalRelError_Intermediate->SetLineWidth(3);
    h_Statistic_Error_Relative_High   ->SetLineWidth(3);

    /*
    // Line to fix 450-525 GV difference
    TLine *l_PRLRelativeSysErrorPRL = new TLine(PublishedPRLBinEdge.back(), 0, PublishedPRLBinEdge.back(), h_PublishedRatioSystematicRelativeErrorPRL.GetBinContent(h_PublishedRatioSystematicRelativeErrorPRL.GetNbinsX()));
    l_PRLRelativeSysErrorPRL->SetLineWidth(3);
    l_PRLRelativeSysErrorPRL->Draw();
    */

    // Axis
    TAxis * xaxis = h_PhysicsReportSystematicRatioRelativeError.GetXaxis();
    TAxis * yaxis = h_PhysicsReportSystematicRatioRelativeError.GetYaxis();
    xaxis->SetTitle("|R| / (GV)");
    yaxis->SetTitle("Relative uncertainty (%)");
    xaxis->SetLimits(1.0, 525);
    yaxis->SetRangeUser(0, 60);
    xaxis->SetTitleFont(62);
    yaxis->SetTitleFont(62);
    xaxis->SetTitleSize(0.045);
    yaxis->SetTitleSize(0.045);
    xaxis->SetLabelFont(62);
    yaxis->SetLabelFont(62);
    xaxis->SetLabelSize(0.05);
    yaxis->SetLabelSize(0.05);
    xaxis->SetMoreLogLabels();

    gPad->SetLogx();
    h_PhysicsReportSystematicRatioRelativeError.SetStats(0);

    p1->SetBottomMargin(0.13);
    xaxis->SetTitleOffset(1.3);
    yaxis->SetTitleOffset(0.9);

    /*
    h_PhysicsReportSystematicRatioRelativeError.Draw("HIST ][");
    h_TotalSysRelError_High        ->Draw("same HIST ][");
    h_TotalSysRelError_Intermediate->Draw("same HIST ][");
    h_TotalSysRelError_Low         ->Draw("same HIST ][");

    h_PhysicsReportStatisticRatioRelativeError.Draw("same HIST ][");
    h_StatisticalRelError_Low         ->Draw("same HIST ][");
    h_StatisticalRelError_Intermediate->Draw("same HIST ][");
    h_Statistic_Error_Relative_High   ->Draw("same HIST ][");
    */

    // Legend
    TLegend *legend = new TLegend(0.25,0.65,0.55,0.85);
    legend->AddEntry(h_Statistic_Error_Relative_High             , "Statistical uncertainty (This analysis)"      , "lpf");
    legend->AddEntry(&h_PhysicsReportStatisticRatioRelativeError , "Statistical uncertainty (Physics Report 2021)", "lpf");
    legend->AddEntry(h_TotalSysRelError_High                     , "Systematic uncertainty (This analysis)"       , "lpf");
    legend->AddEntry(&h_PhysicsReportSystematicRatioRelativeError, "Systematic uncertainty (Physics Report 2021)" , "lpf");
    legend->SetTextSize(0.04);
    legend->SetTextFont(62);
    legend->SetBorderSize(0);
    legend->Draw();

    c_StaSysRelErrCompare.SaveAs( (string("StaSysRelErrCompare_") + issversionname + string(".pdf")).c_str());

}


// Plot total error (absolute)  (FIXME: Deal with overlap range for g_TotalError_Low, Intermediate, High )
void Plot_Total_Error_Absolute(TGraph *g_PublishedRatioToralError, TGraph *g_TotalError_Low, TGraph *g_TotalError_Intermediate, TGraph *g_TotalError_High, TGraph *g_PhysicsReportRatioError, std::string issversionname){
    TCanvas c_totalerr("c_totalerr","c_totalerr",1000,500);
    TPad *pad_totalerr = new TPad("pad_totalerr","pad_totalerr",0.03,0.62,0.50,0.92,32);

    g_PublishedRatioToralError->Draw("AP X");
    g_TotalError_Low->Draw("same P");
    g_TotalError_Intermediate->Draw("same P");
    g_TotalError_High->Draw("same P");
    g_PhysicsReportRatioError->Draw("same P");

    g_PublishedRatioToralError->SetMarkerStyle(15);
    g_PublishedRatioToralError->SetMarkerColor(1);
    g_TotalError_Low->SetMarkerStyle(15);
    g_TotalError_Low->SetMarkerColor(2);
    g_TotalError_Intermediate->SetMarkerStyle(15);
    g_TotalError_Intermediate->SetMarkerColor(2);
    g_TotalError_High->SetMarkerStyle(15);
    g_TotalError_High->SetMarkerColor(2);
    g_PhysicsReportRatioError->SetMarkerStyle(15);
    g_PhysicsReportRatioError->SetMarkerColor(4);

    g_PublishedRatioToralError->SetTitle("");

    TAxis * xaxis_TotErr = g_PublishedRatioToralError->GetXaxis();
    TAxis * yaxis_TotErr = g_PublishedRatioToralError->GetYaxis();
    xaxis_TotErr->SetTitle("Rigidity (GV)");
    yaxis_TotErr->SetTitle("Total Error");
    xaxis_TotErr->SetLimits(1.0,600);
    //yaxis_TotErr->SetRangeUser(0,30);
    xaxis_TotErr->SetTitleFont(62);
    yaxis_TotErr->SetTitleFont(62);
    xaxis_TotErr->SetTitleSize(0.05);
    yaxis_TotErr->SetTitleSize(0.05);
    xaxis_TotErr->SetLabelFont(62);
    yaxis_TotErr->SetLabelFont(62);
    xaxis_TotErr->SetLabelSize(0.07);
    yaxis_TotErr->SetLabelSize(0.07);
    xaxis_TotErr->SetTitleOffset(0.8);
    yaxis_TotErr->SetTitleOffset(0.8);

    TLegend *legend_TotErr = new TLegend(0.35,0.7,0.65,0.85);
    legend_TotErr->AddEntry(g_PublishedRatioToralError,"PRL paper 2016","lpf");
    legend_TotErr->AddEntry(g_TotalError_Low,"This analysis","lpf");
    legend_TotErr->AddEntry(g_PhysicsReportRatioError,"Physics Report 2021","lpf");
    legend_TotErr->SetTextSize(0.05);
    legend_TotErr->SetTextFont(62);
    legend_TotErr->Draw();
    gPad->SetLogx();

    c_totalerr.SaveAs( (string("TotalErr_") + issversionname + string(".pdf")).c_str());
}



// Plot total relative error (Compare with PRL and Physics Report)
void Plot_Total_Relative_Error_CompareReference(TH1D h_PhysicsReportRatioRelativeError, TH1D h_PublishedRatioRelativeErrorPRL, TH1D *h_TotalRelError_Low, TH1D *h_TotalRelError_Intermediate, TH1D *h_TotalRelError_High, std::vector<double> PublishedPRLBinEdge, std::string issversion, std::string issversionname){
    TCanvas c_totalerr_rel("c_totalerr_rel", "c_totalerr_rel", 1000, 500);
    TPad *pad_totalerr_rel = new TPad("pad_totalerr_rel","pad_totalerr_rel", 0.03, 0.62, 0.50, 0.92, 32);

    h_PhysicsReportRatioRelativeError.Draw("HIST ][");
    h_PublishedRatioRelativeErrorPRL .Draw("same HIST ][");

    h_TotalRelError_Low         ->Draw("same HIST ][");
    h_TotalRelError_Intermediate->Draw("same HIST ][");
    h_TotalRelError_High        ->Draw("same HIST ][");

    TLine *l_PRLRelativeTotErrorPRL = new TLine(PublishedPRLBinEdge.back(), 0, PublishedPRLBinEdge.back(), h_PublishedRatioRelativeErrorPRL.GetBinContent(h_PublishedRatioRelativeErrorPRL.GetNbinsX()));
    l_PRLRelativeTotErrorPRL->SetLineWidth(3);
    l_PRLRelativeTotErrorPRL->Draw();

    if (issversion == "2016paper"){
        TLine *l_2016TotalRelativeStaErrorPRL = new TLine(PublishedPRLBinEdge.back(), 0, PublishedPRLBinEdge.back(), h_TotalRelError_High->GetBinContent(h_TotalRelError_High->GetNbinsX()));
        l_2016TotalRelativeStaErrorPRL->SetLineWidth(3);
        l_2016TotalRelativeStaErrorPRL->SetLineColor(2);
        l_2016TotalRelativeStaErrorPRL->Draw();
    }

    h_PhysicsReportRatioRelativeError.SetStats(0);

    h_PublishedRatioRelativeErrorPRL.SetLineColor(1);
    h_PublishedRatioRelativeErrorPRL.SetLineWidth(3);
    h_PhysicsReportRatioRelativeError.SetLineColor(4);
    h_PhysicsReportRatioRelativeError.SetLineWidth(3);
    h_TotalRelError_Low->SetLineColor(2);
    h_TotalRelError_Low->SetLineWidth(3);
    h_TotalRelError_Intermediate->SetLineColor(2);
    h_TotalRelError_Intermediate->SetLineWidth(3);
    h_TotalRelError_High->SetLineColor(2);
    h_TotalRelError_High->SetLineWidth(3);

    TAxis * xaxis_TotRelErr = h_PhysicsReportRatioRelativeError.GetXaxis();
    TAxis * yaxis_TotRelErr = h_PhysicsReportRatioRelativeError.GetYaxis();
    xaxis_TotRelErr->SetTitle("Rigidity (GV)");
    yaxis_TotRelErr->SetTitle("Relative Total Error (%)");
    xaxis_TotRelErr->SetLimits(1.0,525);
    yaxis_TotRelErr->SetRangeUser(0,50);
    xaxis_TotRelErr->SetTitleFont(62);
    yaxis_TotRelErr->SetTitleFont(62);
    xaxis_TotRelErr->SetTitleSize(0.05);
    yaxis_TotRelErr->SetTitleSize(0.05);
    xaxis_TotRelErr->SetLabelFont(62);
    yaxis_TotRelErr->SetLabelFont(62);
    xaxis_TotRelErr->SetLabelSize(0.07);
    yaxis_TotRelErr->SetLabelSize(0.07);
    xaxis_TotRelErr->SetTitleOffset(0.8);
    yaxis_TotRelErr->SetTitleOffset(0.8);

    /*
    TLegend *legend_TotRelErr = new TLegend(0.20, 0.6, 0.70, 0.85);
    if (issversion == "pass7.8"){
        legend_TotRelErr->AddEntry(h_TotalRelError_Low, "This analysis (All time range)"      , "lpf");}
    else if (issversion == "2016paper"){
        legend_TotRelErr->AddEntry(h_TotalRelError_Low, "This analysis (PRL paper time range)", "lpf");}
    else if (issversion == "PhyRep2021"){
        legend_TotRelErr->AddEntry(h_TotalRelError_Low, "This analysis (Phy Rep time range)"  , "lpf");}
    */
    TLegend *legend_TotRelErr = new TLegend(0.20, 0.6, 0.55, 0.85);
    legend_TotRelErr->AddEntry(h_TotalRelError_Low               , "This analysis"      , "lpf");
    legend_TotRelErr->AddEntry(&h_PublishedRatioRelativeErrorPRL , "PRL paper 2016"     , "lpf");
    legend_TotRelErr->AddEntry(&h_PhysicsReportRatioRelativeError, "Physics Report 2021", "lpf");
    legend_TotRelErr->SetTextSize(0.05);
    legend_TotRelErr->SetTextFont(62);
    legend_TotRelErr->Draw();
    gPad->SetLogx();

    c_totalerr_rel.SaveAs( (string("TotalRelErr_") + issversionname + string(".pdf")).c_str());
}



// Plot total relative error (Break Down of my own error result)
void Plot_Total_Relative_Error_ThisAnalysis(TH1D h_PhysicsReportRatioRelativeError, TH1D *h_TotalRelError_High, TH1D *h_TotalRelError_Intermediate, TH1D *h_TotalRelError_Low, TH1D *h_TotalSysRelError_Low, TH1D *h_TotalSysRelError_Intermediate, TH1D *h_TotalSysRelError_High, TH1D *h_StatisticalRelError_Low, TH1D *h_StatisticalRelError_Intermediate, TH1D *h_Statistic_Error_Relative_High, std::string issversionname){

    TCanvas c_totalerr_rel_breakdown("c_totalerr_rel_breakdown","c_totalerr_rel_breakdown",1000,500);
    TPad *pad_totalerr_rel_breakdown = new TPad("pad_totalerr_rel_breakdown","pad_totalerr_rel_breakdown",0, 0, 1, 1, 0);

    h_PhysicsReportRatioRelativeError.Draw("HIST ][");
    h_TotalRelError_High            ->Draw("same HIST ][");
    h_TotalRelError_Intermediate    ->Draw("same HIST ][");
    h_TotalRelError_Low             ->Draw("same HIST ][");

    h_TotalSysRelError_Low         ->Draw("same HIST ][");
    h_TotalSysRelError_Intermediate->Draw("same HIST ][");
    h_TotalSysRelError_High        ->Draw("same HIST ][");

    h_StatisticalRelError_Low         ->Draw("same HIST ][");
    h_StatisticalRelError_Intermediate->Draw("same HIST ][");
    h_Statistic_Error_Relative_High   ->Draw("same HIST ][");

    h_PhysicsReportRatioRelativeError.SetStats(0);

    h_TotalRelError_Low->SetLineColor(1);
    h_TotalRelError_Low->SetLineWidth(3);
    h_TotalRelError_Intermediate->SetLineColor(1);
    h_TotalRelError_Intermediate->SetLineWidth(3);
    h_TotalRelError_High->SetLineColor(1);
    h_TotalRelError_High->SetLineWidth(3);

    h_PhysicsReportRatioRelativeError.SetLineColor(0);
    h_PhysicsReportRatioRelativeError.SetLineWidth(3);

    h_TotalSysRelError_Low->SetLineWidth(3);
    h_TotalSysRelError_Low->SetLineColor(4);
    h_TotalSysRelError_Intermediate->SetLineWidth(3);
    h_TotalSysRelError_Intermediate->SetLineColor(4);
    h_TotalSysRelError_High->SetLineWidth(3);
    h_TotalSysRelError_High->SetLineColor(4);

    h_StatisticalRelError_Low->SetLineWidth(3);
    h_StatisticalRelError_Low->SetLineColor(2);
    h_StatisticalRelError_Intermediate->SetLineWidth(3);
    h_StatisticalRelError_Intermediate->SetLineColor(2);
    h_Statistic_Error_Relative_High->SetLineWidth(3);
    h_Statistic_Error_Relative_High->SetLineColor(2);

    TAxis * xaxis_TotRelErr_BreakDown = h_PhysicsReportRatioRelativeError.GetXaxis();
    TAxis * yaxis_TotRelErr_BreakDown = h_PhysicsReportRatioRelativeError.GetYaxis();
    xaxis_TotRelErr_BreakDown->SetTitle("|R| / (GV)");
    yaxis_TotRelErr_BreakDown->SetTitle("Relative uncertainty (%)");
    xaxis_TotRelErr_BreakDown->SetLimits   (100, 525);
    yaxis_TotRelErr_BreakDown->SetRangeUser(0  , 60);
    xaxis_TotRelErr_BreakDown->SetTitleFont(62);
    yaxis_TotRelErr_BreakDown->SetTitleFont(62);
    xaxis_TotRelErr_BreakDown->SetTitleSize(0.045);
    yaxis_TotRelErr_BreakDown->SetTitleSize(0.045);
    xaxis_TotRelErr_BreakDown->SetLabelFont(62);
    yaxis_TotRelErr_BreakDown->SetLabelFont(62);
    xaxis_TotRelErr_BreakDown->SetLabelSize(0.05);
    yaxis_TotRelErr_BreakDown->SetLabelSize(0.05);
    xaxis_TotRelErr_BreakDown->SetMoreLogLabels();

    gPad->SetBottomMargin(0.18);
    xaxis_TotRelErr_BreakDown->SetTitleOffset(1.3);
    yaxis_TotRelErr_BreakDown->SetTitleOffset(0.8);

    TLegend *legend_TotRelErr_BreakDown = new TLegend(0.20, 0.6, 0.70, 0.85);
    legend_TotRelErr_BreakDown->AddEntry(h_TotalRelError_Low      , "Total uncertainty"      , "lpf");
    legend_TotRelErr_BreakDown->AddEntry(h_StatisticalRelError_Low, "Statistical uncertainty", "lpf");
    legend_TotRelErr_BreakDown->AddEntry(h_TotalSysRelError_Low   , "Systematic uncertainty" , "lpf");
    legend_TotRelErr_BreakDown->SetTextSize(0.04);
    legend_TotRelErr_BreakDown->SetTextFont(62);
    legend_TotRelErr_BreakDown->SetBorderSize(0);
    legend_TotRelErr_BreakDown->Draw();
    gPad->SetLogx();

    c_totalerr_rel_breakdown.SaveAs( (string("TotalRelErr_BreakDown_") + issversionname + string(".pdf")).c_str());
}



// Plot effective acceptance ratios (original, plus, minus)
void Plot_Effective_Acceptance_Ratios_With_Uncertainty_B1220(TGraphErrors *Effective_Acceptance_ratio, TGraphErrors *Effective_Acceptance_ratio_minus10, TGraphErrors *Effective_Acceptance_ratio_plus10){
    TCanvas c_acc("c_acc","c_acc",1000,500);

    Effective_Acceptance_ratio->SetMarkerStyle(15);
    Effective_Acceptance_ratio->SetMarkerColor(2);
    Effective_Acceptance_ratio_minus10->SetMarkerStyle(15);
    Effective_Acceptance_ratio_minus10->SetMarkerColor(3);
    Effective_Acceptance_ratio_plus10->SetMarkerStyle(15);
    Effective_Acceptance_ratio_plus10->SetMarkerColor(4);

    Effective_Acceptance_ratio_plus10->SetTitle("");
    TAxis * xaxis_accerr = Effective_Acceptance_ratio_plus10->GetXaxis();
    TAxis * yaxis_accerr = Effective_Acceptance_ratio_plus10->GetYaxis();
    xaxis_accerr->SetTitle("Rigidity (GV)");
    yaxis_accerr->SetTitle("Effective Acceptance ratio");
    xaxis_accerr->SetLimits(1.0,600);
    yaxis_accerr->SetRangeUser(1.0, 1.4);
    xaxis_accerr->SetTitleFont(62);
    yaxis_accerr->SetTitleFont(62);
    xaxis_accerr->SetTitleSize(0.05);
    yaxis_accerr->SetTitleSize(0.05);
    xaxis_accerr->SetLabelFont(62);
    yaxis_accerr->SetLabelFont(62);
    xaxis_accerr->SetLabelSize(0.07);
    yaxis_accerr->SetLabelSize(0.07);
    xaxis_accerr->SetTitleOffset(0.8);
    yaxis_accerr->SetTitleOffset(0.9);

    gPad->SetLogx();

    Effective_Acceptance_ratio_plus10 ->Draw("AP");
    Effective_Acceptance_ratio_minus10->Draw("same P");
    Effective_Acceptance_ratio        ->Draw("same P");

    TLegend *legend_accerr = new TLegend(0.35,0.6,0.8,0.8);
    // "#frac{A_{p}}{A_{#bar{p}}} (cross section varied +10%)"
    legend_accerr->AddEntry(Effective_Acceptance_ratio_plus10 , "A_{p}/A_{#bar{p}} (cross section varied +10%)","lpf");
    legend_accerr->AddEntry(Effective_Acceptance_ratio_minus10, "A_{p}/A_{#bar{p}} (cross section varied -10%)","lpf");
    legend_accerr->AddEntry(Effective_Acceptance_ratio        , "A_{p}/A_{#bar{p}}","lpf");
    legend_accerr->SetTextSize(0.04);
    legend_accerr->SetTextFont(62);
    legend_accerr->Draw();

    c_acc.SaveAs( (string("EffectiveAcceptanceRatio") + string(".pdf")).c_str());
}


void RemoveErrorbars_TGraphErrors(TGraphErrors *result){
    for (int i=0; i < result->GetN(); i++){
        result->SetPointError(i,0,0);
    }
}

void RemoveErrorbars_TGraphAsymmErrors(TGraphAsymmErrors *result){
    for (int i=0; i < result->GetN(); i++){
        result->SetPointError(i,0,0,0,0);
    }
}

void Plot_Effective_Acceptance_Ratios_With_Uncertainty_B1042(TGraphErrors *Effective_Acceptance_High_B1042, TGraphErrors *Effective_Acceptance_Low_B1042, TGraphErrors *Effective_Acceptance_Intermediate_B1042, TGraphErrors *Effective_Acceptance_Intermediate_inHigh_B1042, TGraphErrors *Effective_Acceptance_ratio_B1042_minus10_AllRange, TGraphErrors *Effective_Acceptance_ratio_B1042_plus10_AllRange, TGraphAsymmErrors *Effective_Acceptance_ratio_B1042_AllRange){

    TCanvas c_acc("c_acc","c_acc",1000,500);

    /*
    // remove fit in acceptance
    Effective_Acceptance_Low_B1042                ->GetListOfFunctions()->Remove(Effective_Acceptance_Low_B1042                ->GetFunction("function1"));
    Effective_Acceptance_Intermediate_B1042       ->GetListOfFunctions()->Remove(Effective_Acceptance_Intermediate_B1042       ->GetFunction("function1"));
    Effective_Acceptance_High_B1042               ->GetListOfFunctions()->Remove(Effective_Acceptance_High_B1042               ->GetFunction("function1"));
    Effective_Acceptance_Intermediate_inHigh_B1042->GetListOfFunctions()->Remove(Effective_Acceptance_Intermediate_inHigh_B1042->GetFunction("function1"));
    */

    //Effective_Acceptance_ratio_B1042_minus10_AllRange->Draw("A C");
    //Effective_Acceptance_ratio_B1042_plus10_AllRange ->Draw("same  C");
    //Effective_Acceptance_High_B1042                  ->Draw("same P");
    //Effective_Acceptance_Low_B1042                   ->Draw("same P");
    //Effective_Acceptance_Intermediate_B1042          ->Draw("same P");
    Effective_Acceptance_ratio_B1042_AllRange        ->Draw("A P3");
    Effective_Acceptance_ratio_B1042_AllRange->SetFillColor(90);
    Effective_Acceptance_ratio_B1042_AllRange->SetFillStyle(1001);

    TGraphAsymmErrors *Effective_Acceptance_ratio_B1042_AllRange_NoErrorBar = new TGraphAsymmErrors(*Effective_Acceptance_ratio_B1042_AllRange);

    RemoveErrorbars_TGraphAsymmErrors(Effective_Acceptance_ratio_B1042_AllRange_NoErrorBar);
    //RemoveErrorbars_TGraphAsymmErrors(Effective_Acceptance_ratio_B1042_AllRange);
    //RemoveErrorbars_TGraphErrors(Effective_Acceptance_High_B1042);
    //RemoveErrorbars_TGraphErrors(Effective_Acceptance_Low_B1042);
    //RemoveErrorbars_TGraphErrors(Effective_Acceptance_Intermediate_B1042);

    Effective_Acceptance_ratio_B1042_AllRange_NoErrorBar->Draw("same P");

    Effective_Acceptance_ratio_B1042_minus10_AllRange->SetMarkerStyle(15);
    Effective_Acceptance_ratio_B1042_minus10_AllRange->SetMarkerColor(3);
    Effective_Acceptance_ratio_B1042_minus10_AllRange->SetMarkerSize(0.9);
    Effective_Acceptance_ratio_B1042_minus10_AllRange->SetLineColor(3);
    Effective_Acceptance_ratio_B1042_plus10_AllRange->SetMarkerStyle(15);
    Effective_Acceptance_ratio_B1042_plus10_AllRange->SetMarkerColor(4);
    Effective_Acceptance_ratio_B1042_plus10_AllRange->SetMarkerSize(0.9); 
    Effective_Acceptance_ratio_B1042_plus10_AllRange->SetLineColor(4);   
    /*
    Effective_Acceptance_High_B1042->SetMarkerStyle(15);
    Effective_Acceptance_High_B1042->SetMarkerColor(2);
    Effective_Acceptance_High_B1042->SetMarkerSize(0.9);
    Effective_Acceptance_Low_B1042->SetMarkerStyle(15);
    Effective_Acceptance_Low_B1042->SetMarkerColor(2);
    Effective_Acceptance_Low_B1042->SetMarkerSize(0.9);
    Effective_Acceptance_Intermediate_B1042->SetMarkerStyle(15);
    Effective_Acceptance_Intermediate_B1042->SetMarkerColor(2);
    Effective_Acceptance_Intermediate_B1042->SetMarkerSize(0.9);
    */
    Effective_Acceptance_ratio_B1042_AllRange->SetMarkerStyle(15);
    Effective_Acceptance_ratio_B1042_AllRange->SetMarkerColor(2);
    Effective_Acceptance_ratio_B1042_AllRange->SetMarkerSize(0.9);
    Effective_Acceptance_ratio_B1042_AllRange_NoErrorBar->SetMarkerStyle(15);
    Effective_Acceptance_ratio_B1042_AllRange_NoErrorBar->SetMarkerColor(2);
    Effective_Acceptance_ratio_B1042_AllRange_NoErrorBar->SetMarkerSize(0.9);

    Effective_Acceptance_ratio_B1042_AllRange->SetTitle("");
    TAxis * xaxis_accerr = Effective_Acceptance_ratio_B1042_AllRange->GetXaxis();
    TAxis * yaxis_accerr = Effective_Acceptance_ratio_B1042_AllRange->GetYaxis();
    xaxis_accerr->SetTitle("|R| / (GV)");
    yaxis_accerr->SetTitle("A_{p} / A_{p}"); 
    xaxis_accerr->SetLimits(1.0, 600);
    yaxis_accerr->SetRangeUser(0.95, 1.4);
    xaxis_accerr->SetTitleFont(62);
    yaxis_accerr->SetTitleFont(62);
    xaxis_accerr->SetTitleSize(0.045);
    yaxis_accerr->SetTitleSize(0.045);
    xaxis_accerr->SetLabelFont(62);
    yaxis_accerr->SetLabelFont(62);
    xaxis_accerr->SetLabelSize(0.05);
    yaxis_accerr->SetLabelSize(0.05);
    xaxis_accerr->SetLabelOffset(0);
    yaxis_accerr->SetLabelOffset(0.02);
    xaxis_accerr->SetTitleOffset(1.2);
    yaxis_accerr->SetTitleOffset(0);

    TLatex latex;
    latex.SetNDC(1);
    latex.SetTextSize(0.053);
    latex.SetTextAngle(90);
    latex.SetTextAlign(13);  //align at top
    latex.DrawLatex(0.064, 0.877, "#minus");

    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.13);
    gPad->SetLogx();
    /*
    TLegend *legend_accerr = new TLegend(0.35,0.6,0.8,0.8);
    legend_accerr->AddEntry(Effective_Acceptance_ratio_plus10 , "A_{p}/A_{#bar{p}} (cross section varied +10%)","lpf");
    legend_accerr->AddEntry(Effective_Acceptance_ratio_minus10, "A_{p}/A_{#bar{p}} (cross section varied -10%)","lpf");
    legend_accerr->AddEntry(Effective_Acceptance_ratio        , "A_{p}/A_{#bar{p}}","lpf");
    legend_accerr->SetTextSize(0.04);
    legend_accerr->SetTextFont(62);
    legend_accerr->Draw();
    */
    c_acc.SaveAs( (string("EffectiveAcceptanceRatio_B1042") + string(".pdf")).c_str());


}













