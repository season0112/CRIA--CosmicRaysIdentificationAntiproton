#include "AntiprotonIntermediateEnergyTree.hh"
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

#include "main.hh"

using namespace std;

#define INFO_OUT_TAG "Unfolding_Bayes_Intermediate"
#include "debugging.hh"

int main(int argc, char* argv[]) {

    Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
    config.ReadCommandLine(argc, argv);

    config.SetProgramHelpText("Unfolding_Bayes_Intermediate",
                            "Bayes Unfolding method");

    config.AddHelpExample("Unfolding_Bayes_Intermediate", "");

    std::string issversion = "";
    config.GetValue("OPTIONS", "issversion", issversion,
                  "The Issversion is");

    chdir( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/") ).c_str());
    string lowpath = ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/") ).c_str(); 


    //// Binning Defination
    std::vector<double> subrangepointused( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+10, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+30); // subrangepointused.at(0): 3.13 (2.97-3.29); subrangepointused.at(19):17.3 (16.6-18.0); subrangepointused.size()=20;
    std::vector<double> subrange_intermediate( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+10, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+31);
    std::vector<double> RigidityBinPoint_Published( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.end());


    //// Load files and Root objects.
    // Open root file for fit result and uncertainty fit result
    TFile *f_ResultOriginal    = new TFile();
    TFile *f_ResultUncertainty = new TFile();
    if (issversion == "pass7.8"){
        f_ResultOriginal    = new TFile("Time_Averaged_ratio/binmerge1/intermediate_0124_free_pass7.8binmerge1.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/intermediate_0124_free_pass7.8binmerge1_uncertainty.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/SysErr_Shape_pass7.8.root");
        f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/SysErr_SigEff_pass7.8.root");
    }
    else if (issversion == "2016paper"){
        f_ResultOriginal    = new TFile("Time_Averaged_ratio/binmerge1/intermediate_0124_free_2016paperbinmerge1.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/intermediate_0124_free_2016paperbinmerge1_uncertainty.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/SysErr_Shape_2016paper.root");
        f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/SysErr_SigEff_2016paper.root");
    }
    else if (issversion == "PhysicsReport"){
        f_ResultOriginal    = new TFile("Time_Averaged_ratio/binmerge1/intermediate_0124_free_PhysicsReportbinmerge1.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/intermediate_0124_free_PhysicsReportbinmerge1_uncertainty.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/SysErr_Shape_PhysicsReport.root");
        f_ResultUncertainty = new TFile("Time_Averaged_ratio/binmerge1/SysErr_SigEff_PhysicsReport.root");
    }

    // Load intermediate result (before unfolding)
    // Temporary object to host
    TGraph *g_chi2dof                  = (TGraph*)f_ResultOriginal->Get("g_chi2dof_TRDeff_0.90");
    TGraph *g_antiproton               = (TGraph*)f_ResultOriginal->Get("g_antiproton_TRDeff_0.90");
    TGraph *g_proton                   = (TGraph*)f_ResultOriginal->Get("g_proton_TRDeff_0.90");
    TGraphErrors *Effective_Acceptance = (TGraphErrors*)f_ResultOriginal->Get("Effective_Acceptance_TRDeff_0.90");
    TGraph *g_error                    = (TGraph*)f_ResultOriginal->Get("g_error_TRDeff_0.90");
    Effective_Acceptance->RemovePoint(61);
    Effective_Acceptance->RemovePoint(0);

    // Load intermediate result with uncertainty
    //TGraph *g_SysError     = (TGraph*)f_ResultUncertainty->Get("g_SysErrorShape");
    TGraph *g_SysError     = (TGraph*)f_ResultUncertainty->Get("g_SysError_SigEff");

    // Fit setting according to mean of ratio distribution
    std::vector<double> *OfficialTRDEff_All;
    f_ResultUncertainty->GetObject( string("OfficialTRDEff").c_str(), OfficialTRDEff_All);

    for (int i=0; i<g_antiproton->GetN(); i++){
        TGraph *g_chi2dof_temporary    = (TGraph*)f_ResultOriginal->Get( (string("g_chi2dof_TRDeff_")    + to_string_with_precision(OfficialTRDEff_All->at(i), 2)).c_str() );
        TGraph *g_antiproton_temporary = (TGraph*)f_ResultOriginal->Get( (string("g_antiproton_TRDeff_") + to_string_with_precision(OfficialTRDEff_All->at(i), 2)).c_str() );
        TGraph *g_proton_temporary     = (TGraph*)f_ResultOriginal->Get( (string("g_proton_TRDeff_")     + to_string_with_precision(OfficialTRDEff_All->at(i), 2)).c_str() );
        TGraph *g_error_temporary      = (TGraph*)f_ResultOriginal->Get( (string("g_error_TRDeff_")      + to_string_with_precision(OfficialTRDEff_All->at(i), 2)).c_str() );

        g_chi2dof    ->SetPoint(i, g_chi2dof   ->GetX()[i], g_chi2dof_temporary->GetY()[i]); 
        g_antiproton ->SetPoint(i, g_antiproton->GetX()[i], g_antiproton_temporary->GetY()[i]);
        g_proton     ->SetPoint(i, g_proton    ->GetX()[i], g_proton_temporary->GetY()[i]);
        g_error      ->SetPoint(i, g_error     ->GetX()[i], g_error_temporary->GetY()[i]);
    }


    // Get Effective Acceptance for Only Intermediate Range 
    TGraphErrors *Effective_Acceptance_IntermediateOnly = new TGraphErrors(20);
    for (int k = 0; k < 20; k++){  // Effective_Acceptance->GetX()[9]=3.13(2.97-3.29GV),  Effective_Acceptance->GetX()[28]=17.3(16.6-18.0GV).
        Effective_Acceptance_IntermediateOnly->SetPoint( k, Effective_Acceptance->GetX()[k+9], Effective_Acceptance->GetY()[k+9] );
        Effective_Acceptance_IntermediateOnly->SetPointError( k, 0.0, Effective_Acceptance->GetErrorY(k+9));
    }

    // Get Effective Acceptance for High Range from Intermediate Range
    TGraphErrors *Effective_Acceptance_Intermediate_inHigh = new TGraphErrors(*Effective_Acceptance);
    for (int i=0; i<27; i++){
        Effective_Acceptance_Intermediate_inHigh->RemovePoint(0);
    }

    // Parametrization of ratio of effective acceptance
    TF1  *function1 = new TF1("function1","[0]*log(log(x))+[1]",0,10);
    Effective_Acceptance->Fit(function1, "", "", 2.97, 18.0 );
    TF1 *fittedfuction1 = Effective_Acceptance->GetFunction("function1");
    TGraph *gFitFunction = new TGraph(fittedfuction1);

    Effective_Acceptance_IntermediateOnly->Fit(function1, "", "", 2.97, 18.0 );
    //TF1 *fittedfuction_IntermediateOnly = Effective_Acceptance_IntermediateOnly->GetFunction("function1");
    //TGraph *gFitFunction_IntermediateOnly = new TGraph(fittedfuction_IntermediateOnly);    
    
    TF1  *function2 = new TF1("function2","[0]*log(log(x))+[1]*log(x)+[2]",14.0,525);
    Effective_Acceptance_Intermediate_inHigh->Fit(function2, "", "", 14.0, 525 );

    // Convert Graph to Histogram
    TH1D *hantiproton_number_raw = new TH1D("", "", 20, subrange_intermediate.data());  
    TH1D *hproton_number_raw     = new TH1D("", "", 20, subrange_intermediate.data());  
    Utilities::ConvertToHistogram ( g_antiproton, *hantiproton_number_raw);
    Utilities::ConvertToHistogram ( g_proton, *hproton_number_raw);


    //// Load MeasuringTime
    TFile *f2 = new TFile("/hpcwork/jara0052/sichen/Measuringtime/MeasuringTime_pass7.8_06_2020_GEOMETRIC35_1.2.root");
    TH1D *fIntegratedMeasuringTimeOverCutOff = (TH1D*)f2->Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff"); 
    TH1D *hMeasuringTime = new TH1D("", "", 20, subrange_intermediate.data()); 
    double MaxMeasuringTime = fIntegratedMeasuringTimeOverCutOff->GetMaximum();


    //// Load EffectiveAcceptance (not EffectiveAcceptance ratio, because for each unfolding, we need EffectiveAcceptance of proton and antiproton respectively.)
    TFile *f3 = new TFile( (lowpath + string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root")).c_str() );
    TFile *f4 = new TFile( (lowpath + string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root")).c_str() );
    TGraphAsymmErrors *Acceptance_antiproton_all = (TGraphAsymmErrors*)f3->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_proton_all     = (TGraphAsymmErrors*)f4->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TH1D *hAcceptance_antiproton                 = new TH1D("", "", 20, subrange_intermediate.data());
    TH1D *hAcceptance_proton                     = new TH1D("", "", 20, subrange_intermediate.data());


    //// Load TriggerEfficiency //FIX ME: Load proton TriggerEff
    TFile *f5 = new TFile( (lowpath + string("/TriggerEff_B1042_antipr.pl1.1800_7.6_all.root")).c_str() );
    TFile *f6 = new TFile( (lowpath + string("/TriggerEff_B1042_antipr.pl1.1800_7.6_all.root")).c_str() );
    TH1F *Trig_antiproton_noprescaling_all = (TH1F*)f5->Get("TriggerEff_noprescaling");
    TH1F *Trig_proton_noprescaling_all     = (TH1F*)f6->Get("TriggerEff_noprescaling");
    TH1D *Trig_antiproton_noprescaling     = new TH1D("", "", 20, subrange_intermediate.data());
    TH1D *Trig_proton_noprescaling         = new TH1D("", "", 20, subrange_intermediate.data());


    //// Load Unfolding_Matrices
    TFile *f7 = new TFile("/hpcwork/jara0052/sichen/Unfolding_Matrices/intermediate/Unfolding_MatricesTH2D_fill.root");
    TH2D *Unfolding_Matrices = (TH2D*)f7->Get("Unfolding_Matrices");


    //// Fill histogram in intermediate range.
    for (int k = 1; k <= 20; ++k){
        hMeasuringTime->SetBinContent(k, fIntegratedMeasuringTimeOverCutOff->GetBinContent(k+10)); //11-30:2.97-18.0
        hAcceptance_antiproton->SetBinContent(k, Acceptance_antiproton_all->GetY()[k+9]); //10-29:2.97-18.0
        hAcceptance_proton->SetBinContent(k, Acceptance_proton_all->GetY()[k+9]);  //10-29:2.97-18.0
        Trig_antiproton_noprescaling->SetBinContent(k, Trig_antiproton_noprescaling_all->GetBinContent(k+10)); //11-30:2.97-18.0
        Trig_proton_noprescaling->SetBinContent(k,  Trig_proton_noprescaling_all->GetBinContent(k+10)); //11-30:2.97-18.0
    }
    //double MaxMeasuringTime = hMeasuringTime->GetBinContent(hMeasuringTime->GetNbinsX());
   

    //// Perform Unfolding
    //const int bayesIterations = 1;
    const int bayesIterations = 6;
    const int verbosity = 0;
    BayesUnfoldingWithCutoff unfolding_antiproton(hantiproton_number_raw, hMeasuringTime, hAcceptance_antiproton, Trig_antiproton_noprescaling, Unfolding_Matrices, MaxMeasuringTime, bayesIterations, verbosity);
    TH1D* hAntiproton_unfolded = unfolding_antiproton.UnfoldedEventCounts();
    BayesUnfoldingWithCutoff unfolding_proton(hproton_number_raw, hMeasuringTime, hAcceptance_proton, Trig_proton_noprescaling, Unfolding_Matrices, MaxMeasuringTime, bayesIterations, verbosity);
    TH1D* hProton_unfolded =unfolding_proton.UnfoldedEventCounts();


    //// Calculate Ratio and Save result in TGraphErrors and TH1D
    TGraphErrors gPublishedRatio = TGraphErrors(57, RigidityBinPoint_Published.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.data(), 0, AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.data());

    TH1D *ratio_raw      = new TH1D("", "", 20, subrange_intermediate.data());
    TH1D *ratio_unfolded = new TH1D("", "", 20, subrange_intermediate.data());
    ratio_raw     ->Divide(hantiproton_number_raw, hproton_number_raw, 1, 1);
    ratio_unfolded->Divide(hAntiproton_unfolded  , hProton_unfolded,   1, 1);

    TH1D *ratio_unfolded_with_effective_correction = new TH1D("", "", 20, subrange_intermediate.data());
    TH1D *ratio_raw_with_effective_correction      = new TH1D("", "", 20, subrange_intermediate.data());
    for (int q = 1; q <= 20; q++){
        ratio_unfolded_with_effective_correction->SetBinContent(q, ratio_unfolded->GetBinContent(q) *  gFitFunction->Eval( subrangepointused.at(q-1) ));
        ratio_raw_with_effective_correction     ->SetBinContent(q, ratio_raw->GetBinContent(q) *  gFitFunction->Eval( subrangepointused.at(q-1) ));
    }

    TGraphErrors *gRatio_unfolded = new TGraphErrors (ratio_unfolded_with_effective_correction);
    TGraphErrors *gRatio_Raw      = new TGraphErrors (ratio_raw_with_effective_correction);


    //// Reset points error with Statistical Error only, Calculate Relative Statistical Error
    TGraphErrors *g_StatisticalRelError = new TGraphErrors(ratio_unfolded_with_effective_correction);
    TH1D *h_StatisticalRelError         = new TH1D("", "", 20, subrange_intermediate.data());
    for (int q = 1; q <= 20; q++){
        gRatio_unfolded->SetPointError( q-1, 0, g_error->GetY()[q-1] );
        gRatio_Raw     ->SetPointError( q-1, 0, g_error->GetY()[q-1] );
        g_StatisticalRelError->SetPoint     (q-1, g_error->GetX()[q-1], g_error->GetY()[q-1]/gRatio_unfolded->GetY()[q-1]*100); // in %.
        h_StatisticalRelError->SetBinContent(q,   g_error->GetY()[q-1]/gRatio_unfolded->GetY()[q-1]*100);                       // in %.
    }


    //// Calculate SystematicError
    TGraphErrors *g_SystematicError_TRD = new TGraphErrors(ratio_unfolded_with_effective_correction);
    TH1D *h_SystematicError_TRD         = new TH1D("", "", 20, subrange_intermediate.data());
    for (int q = 1; q <= 20; q++){
        g_SystematicError_TRD->SetPoint(     q-1, g_error->GetX()[q-1], g_SysError->GetY()[q-1]); 
        g_SystematicError_TRD->SetPointError(q-1, 0, 0);
        h_SystematicError_TRD->SetBinContent(  q, g_SystematicError_TRD->GetY()[q-1]);
    }


    //// Plot unfolded result 
    Plot_Unfolded_Ratio(gRatio_unfolded, gPublishedRatio, issversion);
    Plot_EffectiveAcceptance(Effective_Acceptance);
    Plot_Error(g_error, g_SystematicError_TRD, issversion); 


    //// Save in Root file
    /*
    string ResultFileName;
    if (bayesIterations == 1) {
        ResultFileName = string("unfolded_results_") + issversion;
    }
    else{
        ResultFileName = string("unfolded_results_UnfoldingIteration_") + to_string(bayesIterations) + issversion;
    }
    */
    string ResultFileName = string("unfolded_results_UnfoldingIteration_") + to_string(bayesIterations) + issversion;

    TFile *f = new TFile( (string("Time_Averaged_ratio/unfolding/") + ResultFileName + string(".root")).c_str(),"RECREATE");
    // Chi2
    g_chi2dof                               ->Write("g_chi2dof_intermediate");
    // Error
    g_error                                 ->Write("g_StatisticError");
    g_SystematicError_TRD                   ->Write("g_SystematicError_TRD");
    h_SystematicError_TRD                   ->Write("h_SystematicError_TRD");
    g_StatisticalRelError                   ->Write("g_StatisticalRelError");
    h_StatisticalRelError                   ->Write("h_StatisticalRelError");
    // Number
    hantiproton_number_raw                  ->Write("hantiproton_number_raw");
    hproton_number_raw                      ->Write("hproton_number_raw");
    hAntiproton_unfolded                    ->Write("hAntiproton_number_unfolded");
    hProton_unfolded                        ->Write("hProton_number_unfolded");
    // Acceptance
    Effective_Acceptance                    ->Write("Effective_Acceptance");
    Effective_Acceptance_IntermediateOnly   ->Write("Effective_Acceptance_IntermediateOnly");
    Effective_Acceptance_Intermediate_inHigh->Write("Effective_Acceptance_Intermediate_inHigh");
    // Ratio
    gRatio_unfolded                         ->Write("gRatio_unfolded");
    gRatio_Raw                              ->Write("gRatio_Raw");
    gPublishedRatio.Write("gPublishedRatio");
    f->Close();

}
