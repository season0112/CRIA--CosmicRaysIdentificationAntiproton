//// Note:
//
#include "ExampleAnalysisTree.hh"
#include "BayesUnfoldingWithCutoff.hh"
#include "AnalysisSettings.hh"
// Binning
#include "AntiprotonBinning.hh"
#include "BinningDefinition.hh"
// ACsoft includes
#include "AnalysisEvent.hh"
#include "ConfigHandler.hh"
#include "EventFactory.hh"
#include "FileManager.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "TreeWriter.hh"
#include "Environment.hh"
#include "ObjectManager.hh"
#include "McSpectrumScaler.hh"
#include "TemplateFitter.hh"
#include <iostream>
#include <cassert>
#include <TH2.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TPaletteAxis.h>
#include <TCollection.h>
#include <TTree.h>
#include <TText.h>
#include "TemplateFitter2D.hh"
#include <vector>
#include <TCanvas.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include "IterativeUnfolding.hh"
#include "AcceptanceUnfolding.hh"
#include "Utilities.hh"

#include "Unfolding_Bayes.hh"

using namespace std;
#define INFO_OUT_TAG "Unfolding_Bayes"
#include "debugging.hh"


int main(int argc, char* argv[]) {

    //// Parsing Auguments
    Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
    config.ReadCommandLine(argc, argv);

    config.SetProgramHelpText("Unfolding_Bayes",
                            "Bayes Unfolding method");
    config.AddHelpExample("Unfolding_Bayes", "--binningversion 525version --issversion pass7.8");

    std::string binningversion = "";
    config.GetValue("OPTIONS", "binningversion", binningversion,
                  "The binningversion is");

    std::string issversion = "";
    config.GetValue("OPTIONS", "issversion", issversion,
                  "The issversion is:");

    string pattern = "";
    config.GetValue("OPTIONS", "pattern", pattern,
                  "The choosen tracker pattern is:");

    string ifVGGNN = "";
    config.GetValue("OPTIONS", "ifVGGNN", ifVGGNN,
                  "The choosen of CC estimator is:");

    if (binningversion == "" || issversion == "" || pattern == "" || ifVGGNN == "" ) {
    WARN_OUT << "Some arguments are not given! Please check the help example! " << std::endl;
    return EXIT_FAIL_CONFIG; 
    }

    string NNsuffix;
    string NNsuffixName;
    if ( ifVGGNN == "No" ){
        NNsuffix = ""; 
        NNsuffixName = ""; }
    else if( ifVGGNN == "Yes" ){
        NNsuffix = "_VGG16NN"; 
        NNsuffixName = "VGG16NN"; }
    
    int bins = -1;
    if ( binningversion == "450version" ){
        bins = 31;  // 14.1-450
    }
    else if( binningversion == "525version" ){
        bins = 32;  // 14.1-525
    }


    //// Define binning
    std::vector<double> subrange_450( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+27, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+59);
    std::vector<double> subrange_525( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_zhili525.begin()+27, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_zhili525.begin()+60);
    std::vector<double> subrangecenter_450( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+27, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+58);           // subrangecenter_450.at(0): 14.7 (14.1-15.3);     subrangecenter_450.at(30):354.5 (259-450);      subrangecenter_450.size()=29;
    std::vector<double> subrangecenter_525( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_zhili525.begin()+27, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_zhili525.begin()+59); // subrangecenter_zhili525.at(0):14.7 (14.1-15.3); subrangecenter_zhili525.at(31):427.5 (330-525); subrangecenter_525.size():30

    std::vector<double> subrangeused;
    std::vector<double> subrangepointused;
    if ( binningversion == "450version" ){
        subrangeused.assign(subrange_450.begin(), subrange_450.end());
        subrangepointused.assign(subrangecenter_450.begin(), subrangecenter_450.end());
    }
    else if( binningversion == "525version" ){
        subrangeused.assign(subrange_525.begin(), subrange_525.end());
        subrangepointused.assign(subrangecenter_525.begin(), subrangecenter_525.end());
    }


    //// Load Raw Rersult Before Unfolding.
    string suffix = "";
    TFile *f1 = new TFile();
    if (issversion == "pass7.8"){
        suffix = "";}
    else if(issversion == "published2016"){
        suffix = "_May2015";}
    else if(issversion == "PhyRep2021"){
        suffix = "_Nov2017";}
    chdir(( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/")).c_str() );

    f1 = new TFile( (string("RawRatio_Pattern_") + string(pattern) + NNsuffix + string("_") + string(binningversion) + suffix + string(".root")).c_str() );
    TH1D *hAntiproton_raw       = (TH1D*)f1 ->Get("antiproton_number");
    TH1D *hProton_raw           = (TH1D*)f1 ->Get("proton_number");
    TH1D *hMeasuringTime        = (TH1D*)f1 ->Get("MeasuringTime");
    TH1D *hAcceptance           = (TH1D*)f1 ->Get("Acceptance");
    TH1D *hTriggerEfficiency    = (TH1D*)f1 ->Get("TriggerEfficiency");
    TH1D *hProton_published     = (TH1D*)f1 ->Get("proton_number_published");
    TH1D *hAntiproton_published = (TH1D*)f1 ->Get("antiproton_number_published");
    TTree *tstatistic_error     = (TTree*)f1->Get("tstatistic_error");
    TTree *tdelta_antiproton    = (TTree*)f1->Get("tdelta_antiproton");
    TTree *tsystem_CC = (TTree*)f1->Get("tsystem_CC");
    TGraph *g_fitchi2 = (TGraph*)f1->Get("g_fitchi2");


    //// Load 2016PRL Result
    TGraphErrors gPublishedRatioError;
    TH1D hPublished_Statistic_error            = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Statistic_error_relative   = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Systematic_error           = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Systematic_error_relative  = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_total_error                = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_total_error_relative       = TH1D("", "", 31, subrange_450.data());
    TH1D hPublished_Statistic_error_Proportion = TH1D("", "", 31, subrange_450.data());
    tie(gPublishedRatioError, hPublished_Statistic_error, hPublished_Statistic_error_relative, hPublished_Systematic_error, hPublished_Systematic_error_relative, hPublished_total_error, hPublished_total_error_relative, hPublished_Statistic_error_Proportion) = Load2016PRLResult(subrange_450);


    //// Load MeasuringTime (To calculate MaxMeasuringTime only)
    TFile *fMeasuringTime = new TFile("/hpcwork/jara0052/sichen/Measuringtime/MeasuringTime_pass7.8_06_2020_GEOMETRIC35_1.2.root");
    TH1D *fIntegratedMeasuringTimeOverCutOff = (TH1D*)fMeasuringTime->Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff"); //GetBinLowEdge[2]=1.0;
    double MaxMeasuringTime = fIntegratedMeasuringTimeOverCutOff->GetMaximum();


    //// Load Unfolding Matrices
    TFile *f2 = new TFile( (string("/hpcwork/jara0052/sichen/Unfolding_Matrices/high/") + string("Unfolding_MatricesTH2D_fill_Pattern_") + string(pattern) + string("_") + string(binningversion) + string(".root")).c_str() );
    TH2D *hMigrationMatrix = (TH2D*)f2->Get("Unfolding_Matrices");


    //// Load Effective Acceptance
    TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton;
    TGraphAsymmErrors *g_EffectiveAcceptance_Proton;
    TH1D h_EffectiveAcceptance_Antiproton;
    TH1D h_EffectiveAcceptance_Proton;
    tie(g_EffectiveAcceptance_Antiproton, g_EffectiveAcceptance_Proton, h_EffectiveAcceptance_Antiproton, h_EffectiveAcceptance_Proton) = LoadEffectiveAcceptance(binningversion, pattern, subrange_450, subrange_525);


    //// Calcualte Effective Acceptance Ratio (Manually calculated Error)
    TGraphErrors *g_EffectiveAcceptanceRatio_withManualError;
    tie(g_EffectiveAcceptanceRatio_withManualError) = CalculateAcceptanceRatio(bins, g_EffectiveAcceptance_Antiproton, g_EffectiveAcceptance_Proton);


    //// Fit in Effective_Acceptance_Ratio
    TF1  *function1 = new TF1("function1","[0]*log(log(x))+[1]",14,525);
    g_EffectiveAcceptanceRatio_withManualError->Fit(function1, "", "", 14, 525);
    TF1 *fittedfuction1 = g_EffectiveAcceptanceRatio_withManualError->GetFunction("function1");
    TGraph *gFitFunction = new TGraph(fittedfuction1);

    //// Unfolding
    // Define Unfolding Iterations according to patterns
    const int verbosity = 0;
    int bayesIterations;
    if (pattern == "0"){
        bayesIterations = 4;}
    else if (pattern == "1"){
        bayesIterations = 3;}
    else if (pattern == "2"){
        bayesIterations = 5;}
    else if (pattern == "4"){ 
        bayesIterations = 1;}
    // Perform Unfolding
    BayesUnfoldingWithCutoff unfolding_antiproton(hAntiproton_raw, hMeasuringTime, hAcceptance, hTriggerEfficiency, hMigrationMatrix, MaxMeasuringTime, bayesIterations, verbosity);
    TH1D* hAntiproton_unfolded = unfolding_antiproton.UnfoldedEventCounts();
    BayesUnfoldingWithCutoff unfolding_proton    (hProton_raw    , hMeasuringTime, hAcceptance, hTriggerEfficiency, hMigrationMatrix, MaxMeasuringTime, bayesIterations, verbosity);
    TH1D* hProton_unfolded =unfolding_proton.UnfoldedEventCounts();
    // Calculate Ratio from unfolded Pbar and Proton numbers
    TH1D *ratio_raw       = new TH1D(*hAntiproton_unfolded);
    ratio_raw->Divide(hAntiproton_raw,hProton_raw,1,1);
    TH1D *ratio_unfolded  = new TH1D(*hAntiproton_unfolded);
    ratio_unfolded->Divide(hAntiproton_unfolded,hProton_unfolded,1,1);
    TH1D *ratio_published = new TH1D(*hProton_published);
    ratio_published->Divide(hAntiproton_published,hProton_published,1,1);
    TGraph gRatio_published = TGraph(ratio_published);


    //// Correct Pbar/P Ratio (Using Parameterization of Effective Acceptance Ratio)
    TH1D ratio_unfolded_with_effective_correction = TH1D("","",bins, subrangeused.data());
    for (int q = 1; q <= bins; q++){
        ratio_unfolded_with_effective_correction.SetBinContent(q, ratio_unfolded->GetBinContent(q) *  gFitFunction->Eval( subrangepointused.at(q-1) ));
    }
    TH1D ratio_raw_with_effective_correction = TH1D("","",bins, subrangeused.data());
    for (int q = 1; q <= bins; q++){
        ratio_raw_with_effective_correction.SetBinContent(q, ratio_raw->GetBinContent(q) *  gFitFunction->Eval( subrangepointused.at(q-1) ));
    }


    //// Convert histo to graph for Pbar/Proton Ratio (with and without acceptance correction), and define error.
    Double_t Ratio_unfolded[bins];
    Double_t Ratio_unfoldederror[bins];
    Double_t Ratio_raw[bins];
    Double_t Ratio_rawerror[bins];
    Double_t Ratio_raw_with_acceptanceCorrection[bins];

    Double_t Statistic_error[bins];
    Double_t Delta_antiproton[bins];
    Double_t System_CC[bins];

    Double_t ex_NN[bins] = {0};

    Double_t statistic_error;
    Double_t delta_antiproton;
    Double_t system_CC;

    tstatistic_error->SetBranchAddress("statistic_error",&statistic_error);
    tdelta_antiproton->SetBranchAddress("delta_antiproton",&delta_antiproton);
    tsystem_CC->SetBranchAddress("system_CC",&system_CC);

    int nentries = tstatistic_error->GetEntries();
    for (int i=0; i<nentries; i++) { 
        tstatistic_error->GetEntry(i);
        tdelta_antiproton->GetEntry(i); 
        tsystem_CC->GetEntry(i);
        Statistic_error[i] = statistic_error;
        Delta_antiproton[i] = delta_antiproton;
        System_CC[i] = system_CC;

        //Ratio_unfolded[i]      = ratio_unfolded->GetBinContent(i+1);                                   // plus 1 because TH1 from 1, c++(Statistic_error) from 0.
        Ratio_unfolded[i]      = ratio_unfolded_with_effective_correction.GetBinContent(i+1);            // plus 1 because TH1 from 1, c++(Statistic_error) from 0.
        //Ratio_unfoldederror[i] = ratio_unfolded_with_effective_correction.GetBinError(i+1);            // plus 1 because TH1 from 1, c++(Statistic_error) from 0. 
        //Ratio_unfoldederror[i] = ratio_unfolded->GetBinError(i+1);                                       // plus 1 because TH1 from 1, c++(Statistic_error) from 0.
        Ratio_unfoldederror[i] = ratio_raw->GetBinError(i+1);                                            // plus 1 because TH1 from 1, c++(Statistic_error) from 0.

        Ratio_raw[i]      = ratio_raw->GetBinContent(i+1);                                               // plus 1 because TH1 from 1, c++(Statistic_error) from 0.
        Ratio_rawerror[i] = ratio_raw->GetBinError(i+1);                                                 // plus 1 because TH1 from 1, c++(Statistic_error) from 0.  
        Ratio_raw_with_acceptanceCorrection[i] = ratio_raw_with_effective_correction.GetBinContent(i+1); // plus 1 because TH1 from 1, c++(Statistic_error) from 0. 
    }

    TGraphErrors gRatio_unfolded = TGraphErrors(bins, subrangepointused.data(), Ratio_unfolded, ex_NN, Statistic_error);
    TGraphErrors g_ratio_raw_with_effective_correction = TGraphErrors(bins, subrangepointused.data(), Ratio_raw_with_acceptanceCorrection, ex_NN, Statistic_error); 


    //// Comparision between Raw and Unfolded.
    Double_t raw_unfold_com[bins]; 
    Double_t raw_unfold_com_error[bins];
    for (int j = 0; j < bins; ++j) {
        raw_unfold_com[j]       = (Ratio_raw[j] - Ratio_unfolded[j])/Ratio_unfolded[j] ;
        raw_unfold_com_error[j] = sqrt( pow(Ratio_rawerror[j]/Ratio_unfolded[j], 2) + pow(Ratio_raw[j]*Ratio_unfoldederror[j]/pow(Ratio_unfolded[j],2), 2) ) ;
        cout<< "\n" << endl;
        cout<< "raw:"            << Ratio_raw[j]            << endl;
        cout<< "raw error:"      << Ratio_rawerror[j]       << endl;
        cout<< "unfolded:"       << Ratio_unfolded[j]       << endl;
        cout<< "unfolded error:" << Ratio_unfoldederror[j]  << endl;
        cout<< "rersidual error:"<< raw_unfold_com_error[j] << endl; 
    }
    //TGraph g_Raw_Unfold_Compare = TGraph(bins, subrangepointused.data(), raw_unfold_com);
    TGraphErrors g_Raw_Unfold_Compare = TGraphErrors(bins, subrangepointused.data(), raw_unfold_com, 0, raw_unfold_com_error);

    //// Create TRoot for Statistical Error and Sysmatic Error
    //// Statistical error
    Double_t Statistic_error_relative[bins];
    for (int k = 0; k < bins; ++k) {
        Statistic_error_relative[k] = Statistic_error[k]/Ratio_unfolded[k]*100 ;
    }
    TGraph * gStatistic_error_relative = new TGraph(bins, subrangepointused.data(), Statistic_error_relative); 
    TH1D hStatistic_error_relative = TH1D("", "", bins, subrangeused.data());  
    Utilities::ConvertToHistogram ( gStatistic_error_relative, hStatistic_error_relative);

    TGraph * gStatistic_error = new TGraph(bins, subrangepointused.data(), Statistic_error);
    TH1D hStatistic_error = TH1D("", "", bins, subrangeused.data()); 
    Utilities::ConvertToHistogram ( gStatistic_error, hStatistic_error);

    //// Systetmatic error
    Double_t System_CC_relative[bins];
    for (int k = 0; k < bins; ++k) {
        System_CC_relative[k] = System_CC[k]/Ratio_unfolded[k]*100 ; 
    }
    TGraph * gSystem_CC_relative = new TGraph(bins, subrangepointused.data(), System_CC_relative);
    TH1D hSystem_CC_relative = TH1D("", "", bins, subrangeused.data());
    Utilities::ConvertToHistogram ( gSystem_CC_relative, hSystem_CC_relative);

    TGraph * gSystem_CC = new TGraph(bins, subrangepointused.data(), System_CC);
    TH1D hSystem_CC = TH1D("", "", bins, subrangeused.data());
    Utilities::ConvertToHistogram ( gSystem_CC, hSystem_CC);

    //// Antiproton number uncertainty
    TGraph * gDelta_antiproton = new TGraph(bins, subrangepointused.data(), Delta_antiproton);
    TH1D hDelta_antiproton = TH1D("", "", bins, subrangeused.data()); 
    Utilities::ConvertToHistogram ( gDelta_antiproton, hDelta_antiproton);

     
    //// Plots
    Plot_EffectiveAcceptance_Parametrilised   (g_EffectiveAcceptanceRatio_withManualError, NNsuffix);        
    Plot_UnfoledRatio                         (gRatio_unfolded, g_ratio_raw_with_effective_correction, gRatio_published, issversion, pattern, NNsuffix);
    Plot_Raw_Unfolded_Comparison              (g_Raw_Unfold_Compare, issversion, pattern, NNsuffixName);
    Plot_FitChi2                              (g_fitchi2, issversion, pattern, NNsuffix);
    Plot_MM                                   (hMigrationMatrix, pattern, NNsuffix); // Also Normalize Histogram
    Plot_Statistic_Error                      (hStatistic_error_relative, issversion, pattern, NNsuffix);
    Plot_Statistic_Error_Absolute             (hStatistic_error, issversion, pattern, NNsuffix);
    Plot_Statistic_Error_Compare              (hStatistic_error, hPublished_Statistic_error, issversion, pattern, NNsuffix);
    Plot_Error_Compare                        (hStatistic_error, hPublished_Statistic_error, hPublished_Systematic_error, issversion, pattern, NNsuffix);
    Plot_Error_compare_with_total             (hStatistic_error, hPublished_total_error, hPublished_Statistic_error, hPublished_Systematic_error, issversion, pattern, NNsuffix);
    Plot_Statistic_Error_Compare_Relative     (hStatistic_error_relative, hPublished_Statistic_error_relative, issversion, pattern, NNsuffix);
    Plot_Error_Compare_Relative               (hStatistic_error_relative, hPublished_Statistic_error_relative, hPublished_Systematic_error_relative, issversion, pattern, NNsuffix);
    Plot_Error_Compare_Relative_With_Total    (hStatistic_error_relative, hPublished_Statistic_error_relative, hPublished_Systematic_error_relative, hPublished_total_error_relative, issversion, pattern, NNsuffix);
    Plot_Statistic_Error_Proportion_Published (hPublished_Statistic_error_Proportion, issversion, pattern, NNsuffix);


    //// Save Result.
    TFile *f = new TFile( (string("unfolded_results_Pattern_") + string(pattern) + NNsuffix + suffix + string(".root")).c_str(), "RECREATE");
    // Chi2
    g_fitchi2->Write("g_fitchi2");
    // Numbers
    hAntiproton_unfolded->Write("hAntiproton_unfolded");
    hProton_unfolded    ->Write("hProton_unfolded");
    hAntiproton_raw     ->Write("hAntiproton_raw");
    hProton_raw         ->Write("hProton_raw");
    //Number error
    hDelta_antiproton.Write("Delta_antiproton");
    // Ratio
    ratio_published->Write("ratio_published");
    ratio_raw      ->Write("ratio_raw");                                                        // Raw Ratio, No Effective Acceptance Correction, in TH1D. 
    ratio_raw_with_effective_correction  .Write("ratio_raw_with_effective_correction");         // Raw Ratio, With Effective Acceptance Correction, in TH1D.
    g_ratio_raw_with_effective_correction.Write("g_ratio_raw_with_effective_correction");       // Raw Ratio, With Effective Acceptance Correction, in TGraphErrors.
    ratio_unfolded->Write("ratio_unfolded");                                                    // Unfolded Ratio, No Effective Acceptance Correction, in TH1D.
    ratio_unfolded_with_effective_correction.Write("ratio_unfolded_with_effective_correction"); // Unfolded Ratio, with Effective Acceptance Correction, in TH1D.
    gRatio_unfolded.Write("My_ratio");                                                          // Unfolded Ratio, with Effective Acceptance Correction, in TGraphErrors. 
    g_Raw_Unfold_Compare.Write("g_Raw_Unfold_Compare");
    gPublishedRatioError.Write("g_PRL2016_Ratio");
    // Statistical Error
    gStatistic_error->Write("g_StatisticError");
    hStatistic_error.Write("h_Statistic_error");
    gStatistic_error_relative->Write("g_Statistic_error_relative");
    hStatistic_error_relative.Write("h_Statistic_Error_Relative");
    hPublished_Statistic_error.Write("Published_Statistic_error");
    // Systematic Error
    gSystem_CC->Write("gSystem_CC");
    hSystem_CC.Write("hSystem_CC");
    gSystem_CC_relative->Write("gSystem_CC_relative");
    hSystem_CC_relative.Write("hSystem_CC_relative");
    hPublished_Systematic_error.Write("Published_Systematic_error");
    //Effective Acceptance
    g_EffectiveAcceptance_Antiproton->Write("g_EffectiveAcceptance_Antiproton");   // Effective Acceptance Antiproton (in TGraphAsymmErrors)
    g_EffectiveAcceptance_Proton    ->Write("g_EffectiveAcceptance_Proton");       // Effective Acceptance Proton     (in TGraphAsymmErrors)
    h_EffectiveAcceptance_Antiproton.Write("h_EffectiveAcceptance_Antiproton");    // Effective Acceptance Antiproton (in TH1D)
    h_EffectiveAcceptance_Proton    .Write("h_EffectiveAcceptance_Proton");        // Effective Acceptance Proton     (in TH1D)
    g_EffectiveAcceptanceRatio_withManualError->Write("g_EffectiveAcceptanceRatio_withManualError"); // Effective Acceptance Ratio. (Error is calculated manually.)
    gFitFunction->Write("Fitfunction");  //Parameterization of Effective Acceptance Ratio

    f->Close();

}


