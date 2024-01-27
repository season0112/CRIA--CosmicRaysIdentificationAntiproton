#include "Plot_Ratio_AllPatterns_HighRange.hh"
//// FIXME: Now use P0 (Full Span) AcceptanceCorrection !  
//// FIXME: Temporary use P0 to take place of mixedPattern. (Effective_Acceptance_Ratio_mixedPattern)

void Plot_Ratio_AllPatterns_HighRange(){

    std::string issversion = "pass7.8";
    //std::string issversion = "PhyRep2021";

    std::string suffix;
    if (issversion == "pass7.8"){
        suffix = "";}
    else if (issversion == "2016paper"){
        suffix = "_May2015";}
    else if (issversion == "PhyRep2021"){
        suffix = "_Nov2017";}


    //// Define Patterns Selection range
    int Range1_FirstIndex = 1;
    int Range1_LastIndex  = 12;

    int Range2_FirstIndex = 13;
    int Range2_LastIndex  = 27;

    int Range3_FirstIndex = 28;
    int Range3_LastIndex  = 28;

    int Range4_FirstIndex = 29;
    int Range4_LastIndex  = 32;    

    //// Load PhysicsReport result
    TFile *f_phyreport = new TFile( (string(getenv("MY_ANALYSIS")) + string("/ReferenceFiles/AntiprotonToProtonRatio/AntiprotonToProtonRatio_PhysicsReport/ssdc_canvas.root")).c_str() );
    TGraphErrors *gPhyrReortRatio = (TGraphErrors*)f_phyreport->Get("graph1");
    

    //// Load result in High range 
    string Ratio_Pattern0      = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_0") + suffix + string(".root");
    string Ratio_Pattern0VGGNN = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_0_VGG16NN") + suffix + string(".root");
    string Ratio_Pattern1      = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_1") + suffix + string(".root");
    string Ratio_Pattern2      = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_2") + suffix + string(".root");
    string Ratio_Pattern3      = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_3") + suffix + string(".root");
    string Ratio_Pattern4      = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_4") + suffix + string(".root");
    string Ratio_Pattern5      = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_5") + suffix + string(".root");
    string Ratio_PatternMinus1 = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_Pattern_-1") + suffix + string(".root");

    TFile *f_Ratio_Pattern0 = new TFile(Ratio_Pattern0.c_str());
    TGraphErrors *g_HighResult_P0                         = (TGraphErrors*)f_Ratio_Pattern0->Get("My_ratio");
    TGraph *g_StatisticalError_P0                         = (TGraph*)f_Ratio_Pattern0->Get("g_StatisticError"); 
    TH1D *h_Antiproton_number_unfolded_High_P0            = (TH1D*)f_Ratio_Pattern0->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_P0                = (TH1D*)f_Ratio_Pattern0->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_P0                           = (TH1D*)f_Ratio_Pattern0->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_P0                 = (TH1D*)f_Ratio_Pattern0->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_P0                     = (TH1D*)f_Ratio_Pattern0->Get("hProton_raw");
    TGraph *g_FitChi2_P0                                  = (TGraph*)f_Ratio_Pattern0->Get("g_fitchi2");
    TH1D *h_Statistic_Error_Relative_mixedPattern         = (TH1D*)f_Ratio_Pattern0->Get("h_Statistic_Error_Relative"); // This object is to be refilled to mixedPattern later. No worry here. 
    TGraphErrors *Effective_Acceptance_Ratio_mixedPattern = (TGraphErrors*)f_Ratio_Pattern0->Get("g_EffectiveAcceptanceRatio_withManualError"); // Temporary use P0 to take place of mixedPattern.
    TGraphErrors *Effective_Acceptance_Ratio_P0           = (TGraphErrors*)f_Ratio_Pattern0->Get("g_EffectiveAcceptanceRatio_withManualError");
    //TGraph *g_System_CC_High_mixedPattern                 = (TGraph*)f_Ratio_Pattern0->Get("gSystem_CC"); // Temporary use P0 to take place of mixedPattern.
    TGraphErrors *g_PRL2016_Ratio                         = (TGraphErrors*)f_Ratio_Pattern0->Get("g_PRL2016_Ratio"); 
    TGraph *gFitFunction_P0                               = (TGraph*)f_Ratio_Pattern0->Get("Fitfunction");
    TGraphErrors *g_Raw_Unfold_Compare_P0                 = (TGraphErrors*)f_Ratio_Pattern0->Get("g_Raw_Unfold_Compare");

    TFile *f_Ratio_Pattern0VGGNN = new TFile(Ratio_Pattern0VGGNN.c_str());
    TGraphErrors *g_HighResult_P0VGGNN                = (TGraphErrors*)f_Ratio_Pattern0VGGNN->Get("My_ratio");
    TGraph *g_StatisticalError_P0VGGNN                = (TGraph*)f_Ratio_Pattern0VGGNN->Get("g_StatisticError");
    TH1D *h_Antiproton_number_unfolded_High_P0VGGNN   = (TH1D*)f_Ratio_Pattern0VGGNN->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_P0VGGNN       = (TH1D*)f_Ratio_Pattern0VGGNN->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_P0VGGNN                  = (TH1D*)f_Ratio_Pattern0VGGNN->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_P0VGGNN        = (TH1D*)f_Ratio_Pattern0VGGNN->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_P0VGGNN            = (TH1D*)f_Ratio_Pattern0VGGNN->Get("hProton_raw");
    TGraph *g_FitChi2_P0VGGNN                         = (TGraph*)f_Ratio_Pattern0VGGNN->Get("g_fitchi2");
    TGraphErrors *Effective_Acceptance_Ratio_P0VGGNN  = (TGraphErrors*)f_Ratio_Pattern0VGGNN->Get("g_EffectiveAcceptanceRatio_withManualError");
    TGraph *gFitFunction_P0VGGNN                      = (TGraph*)f_Ratio_Pattern0VGGNN->Get("Fitfunction");
    TGraph *g_System_CC_High_mixedPattern             = (TGraph*)f_Ratio_Pattern0VGGNN->Get("gSystem_CC"); // Temporary use P0 to take place of mixedPattern.
    TGraphErrors *g_Raw_Unfold_Compare_P0VGGNN        = (TGraphErrors*)f_Ratio_Pattern0VGGNN->Get("g_Raw_Unfold_Compare");

    TFile *f_Ratio_Pattern1 = new TFile(Ratio_Pattern1.c_str());
    TGraphErrors *g_HighResult_P1                = (TGraphErrors*)f_Ratio_Pattern1->Get("My_ratio");
    TGraph *g_StatisticalError_P1                = (TGraph*)f_Ratio_Pattern1->Get("g_StatisticError");
    TH1D *h_Antiproton_number_unfolded_High_P1   = (TH1D*)f_Ratio_Pattern1->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_P1       = (TH1D*)f_Ratio_Pattern1->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_P1                  = (TH1D*)f_Ratio_Pattern1->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_P1        = (TH1D*)f_Ratio_Pattern1->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_P1            = (TH1D*)f_Ratio_Pattern1->Get("hProton_raw");
    TGraph *g_FitChi2_P1                         = (TGraph*)f_Ratio_Pattern1->Get("g_fitchi2");
    TGraphErrors *Effective_Acceptance_Ratio_P1  = (TGraphErrors*)f_Ratio_Pattern1->Get("g_EffectiveAcceptanceRatio_withManualError");
    TGraph *gFitFunction_P1                      = (TGraph*)f_Ratio_Pattern1->Get("Fitfunction");
    TGraphErrors *g_Raw_Unfold_Compare_P1        = (TGraphErrors*)f_Ratio_Pattern1->Get("g_Raw_Unfold_Compare");

    TFile *f_Ratio_Pattern2 = new TFile(Ratio_Pattern2.c_str());
    TGraphErrors *g_HighResult_P2                = (TGraphErrors*)f_Ratio_Pattern2->Get("My_ratio");
    TGraph *g_StatisticalError_P2                = (TGraph*)f_Ratio_Pattern2->Get("g_StatisticError");
    TH1D *h_Antiproton_number_unfolded_High_P2   = (TH1D*)f_Ratio_Pattern2->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_P2       = (TH1D*)f_Ratio_Pattern2->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_P2                  = (TH1D*)f_Ratio_Pattern2->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_P2        = (TH1D*)f_Ratio_Pattern2->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_P2            = (TH1D*)f_Ratio_Pattern2->Get("hProton_raw");
    TGraph *g_FitChi2_P2                         = (TGraph*)f_Ratio_Pattern2->Get("g_fitchi2");
    TGraphErrors *Effective_Acceptance_Ratio_P2  = (TGraphErrors*)f_Ratio_Pattern2->Get("g_EffectiveAcceptanceRatio_withManualError");
    TGraph *gFitFunction_P2                      = (TGraph*)f_Ratio_Pattern2->Get("Fitfunction");
    TGraphErrors *g_Raw_Unfold_Compare_P2        = (TGraphErrors*)f_Ratio_Pattern2->Get("g_Raw_Unfold_Compare");

    TFile *f_Ratio_Pattern3 = new TFile(Ratio_Pattern3.c_str());
    TGraphErrors *g_HighResult_P3                = (TGraphErrors*)f_Ratio_Pattern3->Get("My_ratio");
    TGraph *g_StatisticalError_P3                = (TGraph*)f_Ratio_Pattern3->Get("g_StatisticError");
    TH1D *h_Antiproton_number_unfolded_High_P3   = (TH1D*)f_Ratio_Pattern3->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_P3       = (TH1D*)f_Ratio_Pattern3->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_P3                  = (TH1D*)f_Ratio_Pattern3->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_P3        = (TH1D*)f_Ratio_Pattern3->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_P3            = (TH1D*)f_Ratio_Pattern3->Get("hProton_raw");
    TGraph *g_FitChi2_P3                         = (TGraph*)f_Ratio_Pattern3->Get("g_fitchi2");
    TGraphErrors *Effective_Acceptance_Ratio_P3  = (TGraphErrors*)f_Ratio_Pattern3->Get("g_EffectiveAcceptanceRatio_withManualError");
    TGraph *gFitFunction_P3                      = (TGraph*)f_Ratio_Pattern3->Get("Fitfunction");
    TGraphErrors *g_Raw_Unfold_Compare_P3        = (TGraphErrors*)f_Ratio_Pattern3->Get("g_Raw_Unfold_Compare");

    TFile *f_Ratio_Pattern4 = new TFile(Ratio_Pattern4.c_str());
    TGraphErrors *g_HighResult_P4                = (TGraphErrors*)f_Ratio_Pattern4->Get("My_ratio");
    TGraph *g_StatisticalError_P4                = (TGraph*)f_Ratio_Pattern4->Get("g_StatisticError");
    TH1D *h_Antiproton_number_unfolded_High_P4   = (TH1D*)f_Ratio_Pattern4->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_P4       = (TH1D*)f_Ratio_Pattern4->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_P4                  = (TH1D*)f_Ratio_Pattern4->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_P4        = (TH1D*)f_Ratio_Pattern4->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_P4            = (TH1D*)f_Ratio_Pattern4->Get("hProton_raw");
    TGraph *g_FitChi2_P4                         = (TGraph*)f_Ratio_Pattern4->Get("g_fitchi2");
    TGraphErrors *Effective_Acceptance_Ratio_P4  = (TGraphErrors*)f_Ratio_Pattern4->Get("g_EffectiveAcceptanceRatio_withManualError");
    TGraph *gFitFunction_P4                      = (TGraph*)f_Ratio_Pattern4->Get("Fitfunction");
    TGraphErrors *g_Raw_Unfold_Compare_P4        = (TGraphErrors*)f_Ratio_Pattern4->Get("g_Raw_Unfold_Compare");

    TFile *f_Ratio_Pattern5 = new TFile(Ratio_Pattern5.c_str());
    TGraphErrors *g_HighResult_P5                = (TGraphErrors*)f_Ratio_Pattern5->Get("My_ratio");
    TGraph *g_StatisticalError_P5                = (TGraph*)f_Ratio_Pattern5->Get("g_StatisticError");
    TH1D *h_Antiproton_number_unfolded_High_P5   = (TH1D*)f_Ratio_Pattern5->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_P5       = (TH1D*)f_Ratio_Pattern5->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_P5                  = (TH1D*)f_Ratio_Pattern5->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_P5        = (TH1D*)f_Ratio_Pattern5->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_P5            = (TH1D*)f_Ratio_Pattern5->Get("hProton_raw");
    TGraph *g_FitChi2_P5                         = (TGraph*)f_Ratio_Pattern5->Get("g_fitchi2");
    TGraphErrors *Effective_Acceptance_Ratio_P5  = (TGraphErrors*)f_Ratio_Pattern5->Get("g_EffectiveAcceptanceRatio_withManualError");
    TGraph *gFitFunction_P5                      = (TGraph*)f_Ratio_Pattern5->Get("Fitfunction");
    TGraphErrors *g_Raw_Unfold_Compare_P5        = (TGraphErrors*)f_Ratio_Pattern5->Get("g_Raw_Unfold_Compare");

    TFile *f_Ratio_PatternMinus1 = new TFile(Ratio_PatternMinus1.c_str());
    TGraphErrors *g_HighResult_PMinus1                = (TGraphErrors*)f_Ratio_PatternMinus1->Get("My_ratio");
    TGraph *g_StatisticalError_PMinus1                = (TGraph*)f_Ratio_PatternMinus1->Get("g_StatisticError");
    TH1D *h_Antiproton_number_unfolded_High_PMinus1   = (TH1D*)f_Ratio_PatternMinus1->Get("hAntiproton_unfolded");
    TH1D *h_Proton_number_unfolded_High_PMinus1       = (TH1D*)f_Ratio_PatternMinus1->Get("hProton_unfolded");
    TH1D *h_Delta_antiproton_PMinus1                  = (TH1D*)f_Ratio_PatternMinus1->Get("Delta_antiproton");
    TH1D *h_Antiproton_number_Raw_High_PMinus1        = (TH1D*)f_Ratio_PatternMinus1->Get("hAntiproton_raw");
    TH1D *h_Proton_number_Raw_High_PMinus1            = (TH1D*)f_Ratio_PatternMinus1->Get("hProton_raw");
    TGraph *g_FitChi2_PMinus1                         = (TGraph*)f_Ratio_PatternMinus1->Get("g_fitchi2");
    TGraphErrors *Effective_Acceptance_Ratio_PMinus1  = (TGraphErrors*)f_Ratio_PatternMinus1->Get("g_EffectiveAcceptanceRatio_withManualError");
    TGraph *gFitFunction_PMinus1                      = (TGraph*)f_Ratio_PatternMinus1->Get("Fitfunction");
    TGraphErrors *g_Raw_Unfold_Compare_PMinus1        = (TGraphErrors*)f_Ratio_PatternMinus1->Get("g_Raw_Unfold_Compare");    

    //// Calculate Unfolded Ratio without AcceptanceCorrection 
    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN = new TH1D(*h_Antiproton_number_unfolded_High_P0VGGNN);
    h_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN->Divide(h_Antiproton_number_unfolded_High_P0VGGNN, h_Proton_number_unfolded_High_P0VGGNN, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN);

    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_P0 = new TH1D(*h_Antiproton_number_unfolded_High_P0);
    h_Ratio_unfolded_without_AcceptanceCorrection_P0->Divide(h_Antiproton_number_unfolded_High_P0, h_Proton_number_unfolded_High_P0, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P0 = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_P0);

    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_P1 = new TH1D(*h_Antiproton_number_unfolded_High_P1);
    h_Ratio_unfolded_without_AcceptanceCorrection_P1->Divide(h_Antiproton_number_unfolded_High_P1, h_Proton_number_unfolded_High_P1, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P1 = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_P1);

    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_P2 = new TH1D(*h_Antiproton_number_unfolded_High_P2);
    h_Ratio_unfolded_without_AcceptanceCorrection_P2->Divide(h_Antiproton_number_unfolded_High_P2, h_Proton_number_unfolded_High_P2, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P2 = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_P2);

    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_P3 = new TH1D(*h_Antiproton_number_unfolded_High_P3);
    h_Ratio_unfolded_without_AcceptanceCorrection_P3->Divide(h_Antiproton_number_unfolded_High_P3, h_Proton_number_unfolded_High_P3, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P3 = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_P3);

    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_P4 = new TH1D(*h_Antiproton_number_unfolded_High_P4);
    h_Ratio_unfolded_without_AcceptanceCorrection_P4->Divide(h_Antiproton_number_unfolded_High_P4, h_Proton_number_unfolded_High_P4, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P4 = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_P4);

    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_P5 = new TH1D(*h_Antiproton_number_unfolded_High_P5);
    h_Ratio_unfolded_without_AcceptanceCorrection_P5->Divide(h_Antiproton_number_unfolded_High_P5, h_Proton_number_unfolded_High_P5, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_P5 = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_P5);
    
    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_PMinus1 = new TH1D(*h_Antiproton_number_unfolded_High_PMinus1);
    h_Ratio_unfolded_without_AcceptanceCorrection_PMinus1->Divide(h_Antiproton_number_unfolded_High_PMinus1, h_Proton_number_unfolded_High_PMinus1, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_PMinus1 = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_PMinus1);


    //// Pbar and Proton Numbers Patterns Selection in High range according to PRL paper, also Calculate StatisticalError for mixedPattern.
    TH1D *h_Antiproton_number_unfolded_High_mixedPattern = new TH1D(*h_Antiproton_number_unfolded_High_P0);
    TH1D *h_Proton_number_unfolded_High_mixedPattern     = new TH1D(*h_Antiproton_number_unfolded_High_P0);
    TH1D *h_Delta_antiproton_mixedPattern                = new TH1D(*h_Antiproton_number_unfolded_High_P0);
    TH1D *h_Antiproton_number_Raw_High_mixedPattern      = new TH1D(*h_Antiproton_number_unfolded_High_P0);
    TH1D *h_Proton_number_Raw_High_mixedPattern          = new TH1D(*h_Antiproton_number_unfolded_High_P0);
    TGraph *g_StatisticalError_mixedPattern = new TGraph(h_Antiproton_number_unfolded_High_P0);

    // MixedPattern_Selection(Also decide for P0, BDT or VGGNN is used.)
    MixedPattern_Selection(
    Range1_FirstIndex, Range1_LastIndex, Range2_FirstIndex, Range2_LastIndex, Range3_FirstIndex, Range3_LastIndex, Range4_FirstIndex, Range4_LastIndex, 
    h_Antiproton_number_unfolded_High_mixedPattern, h_Proton_number_unfolded_High_mixedPattern, 
    h_Antiproton_number_unfolded_High_P0, h_Antiproton_number_unfolded_High_P0VGGNN, h_Antiproton_number_unfolded_High_P1, h_Antiproton_number_unfolded_High_P2, h_Antiproton_number_unfolded_High_P3, h_Antiproton_number_unfolded_High_P4, h_Antiproton_number_unfolded_High_P5, h_Antiproton_number_unfolded_High_PMinus1, 
    h_Proton_number_unfolded_High_P0, h_Proton_number_unfolded_High_P0VGGNN, h_Proton_number_unfolded_High_P1, h_Proton_number_unfolded_High_P2, h_Proton_number_unfolded_High_P3, h_Proton_number_unfolded_High_P4, h_Proton_number_unfolded_High_P5, h_Proton_number_unfolded_High_PMinus1, 
    g_StatisticalError_mixedPattern, g_HighResult_P0, g_HighResult_P0VGGNN, g_HighResult_P1, g_HighResult_P2, g_HighResult_P3, g_HighResult_P4, g_HighResult_P5, g_HighResult_PMinus1, 
    h_Delta_antiproton_mixedPattern, h_Delta_antiproton_P0, h_Delta_antiproton_P0VGGNN, h_Delta_antiproton_P1, h_Delta_antiproton_P2, h_Delta_antiproton_P3, h_Delta_antiproton_P4, h_Delta_antiproton_P5, h_Delta_antiproton_PMinus1, 
    h_Proton_number_Raw_High_P0, h_Proton_number_Raw_High_P0VGGNN, h_Proton_number_Raw_High_P1, h_Proton_number_Raw_High_P2, h_Proton_number_Raw_High_P3, h_Proton_number_Raw_High_P4, h_Proton_number_Raw_High_P5, h_Proton_number_Raw_High_PMinus1, h_Proton_number_Raw_High_mixedPattern,
    h_Antiproton_number_Raw_High_P0, h_Antiproton_number_Raw_High_P0VGGNN, h_Antiproton_number_Raw_High_P1, h_Antiproton_number_Raw_High_P2, h_Antiproton_number_Raw_High_P3, h_Antiproton_number_Raw_High_P4, h_Antiproton_number_Raw_High_P5, h_Antiproton_number_Raw_High_PMinus1, h_Antiproton_number_Raw_High_mixedPattern, "VGGNN"); // last argument decide BDT or VGGNN is used. 


    //// Calculate Unfolded Ratio (mixedPattern) from Numbers.  --> No Acceptance Correction.
    TH1D *h_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern = new TH1D(*h_Antiproton_number_unfolded_High_P0);
    h_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->Divide(h_Antiproton_number_unfolded_High_mixedPattern, h_Proton_number_unfolded_High_mixedPattern, 1, 1);
    TGraphErrors *g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern = new TGraphErrors(h_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern);


    //// Calculate Unfolded Ratio (mixedPattern)                --> With Acceptance Correction.
    // FIXME: Now use P0 (Full Span) AcceptanceCorrection !
    TGraphErrors *g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern = new TGraphErrors(*g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern);
    for (int q = 0; q < g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetN(); q++){
        g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->SetPoint(q, g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetX()[q] ,g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetY()[q] * gFitFunction_P0->Eval(g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetX()[q]));
    }


    //// Reset Error for RatioFromUnfoldedNumbers (without AcceptanceCorrection)
    Reset_StatisticalError(g_HighResult_P0, g_HighResult_P0VGGNN, g_HighResult_P1, g_HighResult_P2, g_HighResult_P3, g_HighResult_P4, g_HighResult_P5, g_HighResult_PMinus1, 
    g_Ratio_unfolded_without_AcceptanceCorrection_P0, g_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN, g_Ratio_unfolded_without_AcceptanceCorrection_P1, g_Ratio_unfolded_without_AcceptanceCorrection_P2, g_Ratio_unfolded_without_AcceptanceCorrection_P3, g_Ratio_unfolded_without_AcceptanceCorrection_P4, g_Ratio_unfolded_without_AcceptanceCorrection_P5, g_Ratio_unfolded_without_AcceptanceCorrection_PMinus1, 
    g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern, g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern, g_StatisticalError_mixedPattern);


    //// Plot
    Plot_RawRatioAndUnfoldedRatio_Compare(issversion, g_Raw_Unfold_Compare_P0VGGNN, g_Raw_Unfold_Compare_P1, g_Raw_Unfold_Compare_P2, g_Raw_Unfold_Compare_P4);
     
    Plot_EffectiveAcceptanceRatio_AllPatterns(gFitFunction_P0, gFitFunction_P1, gFitFunction_P2, gFitFunction_P3, gFitFunction_P4, gFitFunction_P5, gFitFunction_PMinus1, Effective_Acceptance_Ratio_P0VGGNN, Effective_Acceptance_Ratio_P1, Effective_Acceptance_Ratio_P2, Effective_Acceptance_Ratio_P3, Effective_Acceptance_Ratio_P4, Effective_Acceptance_Ratio_P5, Effective_Acceptance_Ratio_PMinus1);
    Plot_EffectiveAcceptanceRatio_SinglePattern(Effective_Acceptance_Ratio_P0VGGNN, "Pattern0");
    Plot_EffectiveAcceptanceRatio_SinglePattern(Effective_Acceptance_Ratio_P1     , "Pattern1");
    Plot_EffectiveAcceptanceRatio_SinglePattern(Effective_Acceptance_Ratio_P2     , "Pattern2");
    Plot_EffectiveAcceptanceRatio_SinglePattern(Effective_Acceptance_Ratio_P4     , "Pattern4");

    Plot_Ratio_HighRange_AllPattern(gPhyrReortRatio, g_HighResult_P0, g_HighResult_P0VGGNN, g_HighResult_P1, g_HighResult_P2, g_HighResult_P3, g_HighResult_P4, g_HighResult_P5, g_HighResult_PMinus1, issversion, "WithAcceptanceCorrection", g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern);
    Plot_Ratio_HighRange_AllPattern(gPhyrReortRatio, g_Ratio_unfolded_without_AcceptanceCorrection_P0, g_Ratio_unfolded_without_AcceptanceCorrection_P0VGGNN, g_Ratio_unfolded_without_AcceptanceCorrection_P1, g_Ratio_unfolded_without_AcceptanceCorrection_P2, g_Ratio_unfolded_without_AcceptanceCorrection_P3, g_Ratio_unfolded_without_AcceptanceCorrection_P4, g_Ratio_unfolded_without_AcceptanceCorrection_P5, g_Ratio_unfolded_without_AcceptanceCorrection_PMinus1, issversion, "WithoutAcceptanceCorrection", g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern);

    Plot_Ratio_HighRange_mixedPattern(g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern, gPhyrReortRatio, g_PRL2016_Ratio, g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern);

    Plot_Ratio_HighRange_PatternSelection(gPhyrReortRatio, g_PRL2016_Ratio, Range1_FirstIndex, Range1_LastIndex, Range2_FirstIndex, Range2_LastIndex, Range3_FirstIndex, Range3_LastIndex, Range4_FirstIndex, Range4_LastIndex, g_HighResult_P0, g_HighResult_P0VGGNN, g_HighResult_P1, g_HighResult_P2, g_HighResult_P3, g_HighResult_P4, g_HighResult_P5, issversion, "WithAcceptanceCorrection", g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern); //Here the pattern selection applies to g_Ratio_unfolded_without(with)_AcceptanceCorrection_PX !!! (Namely RemovePoints)

    Plot_Ratio_HighRange_SinglePatternCompareReference(gPhyrReortRatio, g_HighResult_P0     , "Pattern0"     , issversion, "WithAcceptanceCorrection");
    Plot_Ratio_HighRange_SinglePatternCompareReference(gPhyrReortRatio, g_HighResult_P0VGGNN, "Pattern0VGGNN", issversion, "WithAcceptanceCorrection");
    Plot_Ratio_HighRange_SinglePatternCompareReference(gPhyrReortRatio, g_HighResult_P1     , "Pattern1"     , issversion, "WithAcceptanceCorrection");
    Plot_Ratio_HighRange_SinglePatternCompareReference(gPhyrReortRatio, g_HighResult_P2     , "Pattern2"     , issversion, "WithAcceptanceCorrection");
    Plot_Ratio_HighRange_SinglePatternCompareReference(gPhyrReortRatio, g_HighResult_P4     , "Pattern4"     , issversion, "WithAcceptanceCorrection");

    Plot_StatisticalError_ComparePatterns(g_StatisticalError_P0, g_StatisticalError_P0VGGNN, g_StatisticalError_P1, g_StatisticalError_P2, g_StatisticalError_P3, g_StatisticalError_P4, g_StatisticalError_P5, g_StatisticalError_mixedPattern, h_Antiproton_number_unfolded_High_P0, h_Antiproton_number_unfolded_High_P0VGGNN, h_Antiproton_number_unfolded_High_P1, h_Antiproton_number_unfolded_High_P2, h_Antiproton_number_unfolded_High_P3, h_Antiproton_number_unfolded_High_P4, h_Antiproton_number_unfolded_High_P5);
    //Plot_StatisticalErrorCompareWithSquareRootN(h_Antiproton_number_unfolded_High_P0, h_Antiproton_number_unfolded_High_P0VGGNN, h_Antiproton_number_unfolded_High_P1, h_Antiproton_number_unfolded_High_P2, h_Antiproton_number_unfolded_High_P3, h_Antiproton_number_unfolded_High_P4, h_Antiproton_number_unfolded_High_P5, h_Delta_antiproton_P0, h_Delta_antiproton_P0VGGNN, h_Delta_antiproton_P1, h_Delta_antiproton_P2, h_Delta_antiproton_P3, h_Delta_antiproton_P4, h_Delta_antiproton_P5);  // This funciton will break histogram content, need to be fixed before active it.


    Plot_Chi2(Range1_FirstIndex, Range1_LastIndex, Range2_FirstIndex, Range2_LastIndex, Range3_FirstIndex, Range3_LastIndex, Range4_FirstIndex, Range4_LastIndex, g_FitChi2_P0, g_FitChi2_P1, g_FitChi2_P2, g_FitChi2_P3, g_FitChi2_P4, g_FitChi2_P5);

    //// Make Relative Statistic Error
    for (int i = 1; i<=h_Statistic_Error_Relative_mixedPattern->GetNbinsX(); i++){ 
        h_Statistic_Error_Relative_mixedPattern->SetBinContent(i, g_StatisticalError_mixedPattern->GetY()[i-1]/g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern->GetY()[i-1]*100);
    }


    ////Save
    TFile checkrootfile( (string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/") + string("unfolded_results_MixedPattern") + suffix + string(".root")).c_str(),"RECREATE");
    // Ratio
    g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern->Write("g_Ratio_unfolded_without_AcceptanceCorrection_mixedPattern");
    g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern   ->Write("g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern");
    g_HighResult_P0             ->Write("g_HighResult_P0");
    g_HighResult_P0VGGNN        ->Write("g_HighResult_P0VGGNN");
    g_HighResult_P1             ->Write("g_HighResult_P1");
    g_HighResult_P2             ->Write("g_HighResult_P2");
    g_HighResult_P3             ->Write("g_HighResult_P3");
    g_HighResult_P4             ->Write("g_HighResult_P4");
    g_HighResult_P5             ->Write("g_HighResult_P5");
    g_Raw_Unfold_Compare_P0VGGNN->Write("g_Raw_Unfold_Compare_P0VGGNN");
    g_Raw_Unfold_Compare_P1     ->Write("g_Raw_Unfold_Compare_P1");
    g_Raw_Unfold_Compare_P2     ->Write("g_Raw_Unfold_Compare_P2"); 
    g_Raw_Unfold_Compare_P4     ->Write("g_Raw_Unfold_Compare_P4");
    // Numbers
    h_Antiproton_number_unfolded_High_mixedPattern->Write("h_Antiproton_number_unfolded_High_mixedPattern");
    h_Proton_number_unfolded_High_mixedPattern    ->Write("h_Proton_number_unfolded_High_mixedPattern");
    h_Antiproton_number_Raw_High_mixedPattern ->Write("h_Antiproton_number_Raw_High_mixedPattern");
    h_Proton_number_Raw_High_mixedPattern->Write("h_Proton_number_Raw_High_mixedPattern");

    h_Antiproton_number_unfolded_High_P0          ->Write("h_Antiproton_number_unfolded_High_P0");
    h_Antiproton_number_unfolded_High_P0VGGNN     ->Write("h_Antiproton_number_unfolded_High_P0VGGNN");
    h_Antiproton_number_unfolded_High_P1          ->Write("h_Antiproton_number_unfolded_High_P1");
    h_Antiproton_number_unfolded_High_P2          ->Write("h_Antiproton_number_unfolded_High_P2");
    h_Antiproton_number_unfolded_High_P3          ->Write("h_Antiproton_number_unfolded_High_P3");
    h_Antiproton_number_unfolded_High_P4          ->Write("h_Antiproton_number_unfolded_High_P4");
    h_Antiproton_number_unfolded_High_P5          ->Write("h_Antiproton_number_unfolded_High_P5");
 
    h_Proton_number_unfolded_High_P0     ->Write("h_Proton_number_unfolded_High_P0");        
    h_Proton_number_unfolded_High_P0VGGNN->Write("h_Proton_number_unfolded_High_P0VGGNN");
    h_Proton_number_unfolded_High_P1     ->Write("h_Proton_number_unfolded_High_P1");
    h_Proton_number_unfolded_High_P2     ->Write("h_Proton_number_unfolded_High_P2");
    h_Proton_number_unfolded_High_P3     ->Write("h_Proton_number_unfolded_High_P3");
    h_Proton_number_unfolded_High_P4     ->Write("h_Proton_number_unfolded_High_P4");
    h_Proton_number_unfolded_High_P5     ->Write("h_Proton_number_unfolded_High_P5");

    h_Antiproton_number_Raw_High_P0VGGNN ->Write("h_Antiproton_number_Raw_High_P0VGGNN");
    h_Antiproton_number_Raw_High_P0      ->Write("h_Antiproton_number_Raw_High_P0");
    h_Antiproton_number_Raw_High_P1      ->Write("h_Antiproton_number_Raw_High_P1");
    h_Antiproton_number_Raw_High_P2      ->Write("h_Antiproton_number_Raw_High_P2");
    h_Antiproton_number_Raw_High_P3      ->Write("h_Antiproton_number_Raw_High_P3");
    h_Antiproton_number_Raw_High_P4      ->Write("h_Antiproton_number_Raw_High_P4");
    h_Antiproton_number_Raw_High_P5      ->Write("h_Antiproton_number_Raw_High_P5");

    h_Proton_number_Raw_High_P0VGGNN     ->Write("h_Proton_number_Raw_High_P0VGGNN");
    h_Proton_number_Raw_High_P0          ->Write("h_Proton_number_Raw_High_P0");
    h_Proton_number_Raw_High_P1          ->Write("h_Proton_number_Raw_High_P1");
    h_Proton_number_Raw_High_P2          ->Write("h_Proton_number_Raw_High_P2");
    h_Proton_number_Raw_High_P3          ->Write("h_Proton_number_Raw_High_P3");
    h_Proton_number_Raw_High_P4          ->Write("h_Proton_number_Raw_High_P4");
    h_Proton_number_Raw_High_P5          ->Write("h_Proton_number_Raw_High_P5");

    // Acceptance
    Effective_Acceptance_Ratio_mixedPattern->Write("Effective_Acceptance_mixedPattern");
    Effective_Acceptance_Ratio_P0->Write("Effective_Acceptance_P0");
    Effective_Acceptance_Ratio_P1->Write("Effective_Acceptance_P1");
    Effective_Acceptance_Ratio_P2->Write("Effective_Acceptance_P2");
    Effective_Acceptance_Ratio_P4->Write("Effective_Acceptance_P4");
    // Error
    g_StatisticalError_mixedPattern->Write("g_StatisticalError_mixedPattern");
    h_Statistic_Error_Relative_mixedPattern->Write("h_Statistic_Error_Relative_mixedPattern");
    g_System_CC_High_mixedPattern->Write("g_System_CC_High_mixedPattern");

    checkrootfile.Close();    




}
