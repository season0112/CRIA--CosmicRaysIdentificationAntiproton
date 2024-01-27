#include "AntiprotonLowEnergyTree.hh"
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

#define INFO_OUT_TAG "Unfolding_Bayes_Low"
#include "debugging.hh"

int main(int argc, char* argv[]) {

    Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
    config.ReadCommandLine(argc, argv);

    config.SetProgramHelpText("Unfolding_Bayes_Low",
                              "Bayes Unfolding method");

    config.AddHelpExample("Unfolding_Bayes_Low", "");

    std::string issversion = "";
    config.GetValue("OPTIONS", "issversion", issversion,
                    "The Issversion is");


    chdir( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/") ).c_str());
    int startindex=1;
    int endindex=18;


    //// Define binnings
    std::vector<double> subrange_low( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+startindex, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+endindex); // AntiprotonBinning_450.begin()+1: 1.0-1.16GV;
    std::vector<double> subrangepointused( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+startindex, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+endindex-1); // AntiprotonBinningCenter_450.begin()+1: 1.08GV;    AntiprotonBinningCenter_450.begin()+16: 5.635GV; subrangepointused.size()=16;
    std::vector<double> RigidityBinPoint_Published( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.end());
    std::vector<double> publishedcenter ( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.end()); // since first point (0.8-1.0) don't have result, therefore it should be removed.
    std::vector<double> PhysicsReportCenter ( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_525.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_525.end()); // since first point (0.8-1.0) don't have result, therefore it should be removed.

    //// Load files and Root objects.

    // Load low result without unfolding 
    TFile *f_ResultOriginal = new TFile();
    TFile *f_ResultUncertainty = new TFile();
    if (issversion == "pass7.8"){
        f_ResultOriginal    = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/Ratio_pass7.8.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/SysErr_Shape_pass7.8.root");
        f_ResultUncertainty = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/SysErr_SigEff_pass7.8.root");
    }
    else if (issversion == "2016paper"){
        f_ResultOriginal    = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/Ratio_2016paper.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/SysErr_Shape_2016paper.root");
        f_ResultUncertainty = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/SysErr_SigEff_2016paper.root");
    }
    else if (issversion == "PhysicsReport"){
        f_ResultOriginal    = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/Ratio_PhysicsReport.root");
        //f_ResultUncertainty = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/SysErr_Shape_PhysicsReport.root");
        f_ResultUncertainty = new TFile("Time_Averaged_ratio_Low/binmerge1/plots/SysErr_SigEff_PhysicsReport.root");
    }

    // Temporary object to host    
    TGraph *g_antiproton_number = (TGraph*)f_ResultOriginal      ->Get("g_antiproton_number_TRDeff_0.94_TOFeff_0.95");
    TGraph *g_proton_number     = (TGraph*)f_ResultOriginal      ->Get("g_proton_number_TRDeff_0.94_TOFeff_0.95");
    TGraphErrors *g_error       = (TGraphErrors*)f_ResultOriginal->Get("g_error_TRDeff_0.94_TOFeff_0.95");
    TGraphErrors *g_chi2_tof    = (TGraphErrors*)f_ResultOriginal->Get("chi2_tof_TRDeff_0.94_TOFeff_0.95");
    /*
    TGraph *g_antiproton_number = (TGraph*)f_ResultOriginal      ->Get("g_antiproton_number_TRDeff_0.84_TOFeff_0.72");
    TGraph *g_proton_number     = (TGraph*)f_ResultOriginal      ->Get("g_proton_number_TRDeff_0.84_TOFeff_0.72");
    TGraphErrors *g_error       = (TGraphErrors*)f_ResultOriginal->Get("g_error_TRDeff_0.84_TOFeff_0.72");
    TGraphErrors *g_chi2_tof    = (TGraphErrors*)f_ResultOriginal->Get("chi2_tof_TRDeff_0.84_TOFeff_0.72");
    */

    /* 
    // Fit setting for first four points
    for (int i=0; i<2; i++){
        TGraph *g_antiproton_number_temporary = (TGraph*)f_ResultOriginal      ->Get( (string("g_antiproton_number_TRDeff_0.84") + string("_TOFeff_0.72") ).c_str() );
        TGraph *g_proton_number_temporary     = (TGraph*)f_ResultOriginal      ->Get( (string("g_proton_number_TRDeff_0.84")     + string("_TOFeff_0.72") ).c_str() );
        //TGraphErrors *g_error_temporary       = (TGraphErrors*)f_ResultOriginal->Get( (string("g_error_TRDeff_0.84")             + string("_TOFeff_0.72") ).c_str() );
        //TGraphErrors *g_chi2_tof_temporary    = (TGraphErrors*)f_ResultOriginal->Get( (string("chi2_tof_TRDeff_0.84")            + string("_TOFeff_0.72") ).c_str() );

        g_antiproton_number->SetPoint(i, g_antiproton_number->GetX()[i], g_antiproton_number_temporary->GetY()[i]);
        g_proton_number    ->SetPoint(i, g_proton_number->GetX()[i]    , g_proton_number_temporary->GetY()[i]);
        //g_error            ->SetPoint(i, g_error->GetX()[i]            , g_error_temporary->GetY()[i]);
        //g_chi2_tof         ->SetPoint(i, g_chi2_tof->GetX()[i]         , g_chi2_tof_temporary->GetY()[i]);
    }
    */

    /*    
    // Fit setting according to mean of ratio distribution
    std::vector<double> *OfficialTRDEff_All;
    std::vector<double> *OfficialTOFEff_All;
    f_ResultUncertainty->GetObject( string("OfficialTRDEff").c_str(), OfficialTRDEff_All);
    f_ResultUncertainty->GetObject( string("OfficialTOFEff").c_str(), OfficialTOFEff_All);

    for (int i=0; i<g_antiproton_number->GetN(); i++){
        TGraph *g_antiproton_number_temporary = (TGraph*)f_ResultOriginal->Get( (string("g_antiproton_number_TRDeff_") + to_string_with_precision(OfficialTRDEff_All->at(i), 2) + string("_TOFeff_") + to_string_with_precision(OfficialTOFEff_All->at(i), 2)).c_str() );
        TGraph *g_proton_number_temporary     = (TGraph*)f_ResultOriginal->Get( (string("g_proton_number_TRDeff_") + to_string_with_precision(OfficialTRDEff_All->at(i), 2) + string("_TOFeff_") + to_string_with_precision(OfficialTOFEff_All->at(i), 2)).c_str() );
        TGraphErrors *g_error_temporary       = (TGraphErrors*)f_ResultOriginal->Get( (string("g_error_TRDeff_") + to_string_with_precision(OfficialTRDEff_All->at(i), 2) + string("_TOFeff_") + to_string_with_precision(OfficialTOFEff_All->at(i), 2)).c_str() );
        TGraphErrors *g_chi2_tof_temporary    = (TGraphErrors*)f_ResultOriginal->Get( (string("chi2_tof_TRDeff_") + to_string_with_precision(OfficialTRDEff_All->at(i), 2) + string("_TOFeff_") + to_string_with_precision(OfficialTOFEff_All->at(i), 2)).c_str() );

        g_antiproton_number->SetPoint(i, g_antiproton_number->GetX()[i], g_antiproton_number_temporary->GetY()[i]);
        g_proton_number    ->SetPoint(i, g_proton_number->GetX()[i]    , g_proton_number_temporary->GetY()[i]); 
        g_error            ->SetPoint(i, g_error->GetX()[i]            , g_error_temporary->GetY()[i]);
        g_chi2_tof         ->SetPoint(i, g_chi2_tof->GetX()[i]         , g_chi2_tof_temporary->GetY()[i]);
    }
    */


    // Load low result with uncertainty
    //TGraph *g_SysError     = (TGraph*)f_ResultUncertainty->Get("g_SysErrorShape");
    TGraph *g_SysError     = (TGraph*)f_ResultUncertainty->Get("g_SysError_SigEff");


    // Convert graph to histogram
    TH1D *hantiproton = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    TH1D *hproton     = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    Utilities::ConvertToHistogram ( g_antiproton_number, *hantiproton);
    Utilities::ConvertToHistogram ( g_proton_number, *hproton);


    //// Load Published Result and Error
    TGraphErrors gPublishedRatio = TGraphErrors(57, RigidityBinPoint_Published.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.data(), 0, AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.data());    
    // PRL Result
    TGraph *PublishedError                           = new TGraph(AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.size(), publishedcenter.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.data() );
    TGraph *PublishedRatioStatisticErrorPRL          = new TGraph(AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticErrorPRL.size(), publishedcenter.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticErrorPRL.data() );
    TGraph *PublishedRatioSystematicErrorPRL         = new TGraph(AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicErrorPRL.size(), publishedcenter.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicErrorPRL.data() );
    TGraph *PublishedRatioStatisticRelativeErrorPRL  = new TGraph(AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.size(), publishedcenter.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.data() );
    TGraph *PublishedRatioSystematicRelativeErrorPRL = new TGraph(AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.size(), publishedcenter.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.data() );
    TGraph *PublishedRatioRelativeErrorPRL           = new TGraph(AntiprotonNewBinning::AntiprotonResults::PublishedRatioRelativeErrorPRL.size(), publishedcenter.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioRelativeErrorPRL.data() );
    // PhysicsReport Result
    TGraph *PhysicsReportRatioError                   = new TGraph(AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioError.size(), PhysicsReportCenter.data(), AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioError.data() );
    TGraph *PhysicsReportRatioStatisticError          = new TGraph(AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioStatisticError.size(), PhysicsReportCenter.data(), AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioStatisticError.data() );
    TGraph *PhysicsReportRatioSystematicError         = new TGraph(AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioSystematicError.size(), PhysicsReportCenter.data(), AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioSystematicError.data() );
    TGraph *PhysicsReportRatioRelativeError           = new TGraph(AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioRelativeError.size(), PhysicsReportCenter.data(), AntiprotonNewBinning::AntiprotonResults::PhysicsReportRatioRelativeError.data() );
    TGraph *PhysicsReportStatisticRatioRelativeError  = new TGraph(AntiprotonNewBinning::AntiprotonResults::PhysicsReportStatisticRatioRelativeError.size(), PhysicsReportCenter.data(), AntiprotonNewBinning::AntiprotonResults::PhysicsReportStatisticRatioRelativeError.data() );
    TGraph *PhysicsReportSystematicRatioRelativeError = new TGraph(AntiprotonNewBinning::AntiprotonResults::PhysicsReportSystematicRatioRelativeError.size(), PhysicsReportCenter.data(), AntiprotonNewBinning::AntiprotonResults::PhysicsReportSystematicRatioRelativeError.data() );


    //// Load MeasuringTime
    TFile *f2 = new TFile("/hpcwork/jara0052/sichen/Measuringtime/MeasuringTime_pass7.8_06_2020_GEOMETRIC35_1.2.root");
    TH1D *fIntegratedMeasuringTimeOverCutOff = (TH1D*)f2->Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff"); //GetBinLowEdge[2]=1.0;
    TH1D *hMeasuringTime = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    double MaxMeasuringTime = fIntegratedMeasuringTimeOverCutOff->GetMaximum();


    //// Load EffectiveAcceptance (not EffectiveAcceptance ratio, because for each unfolding, we need EffectiveAcceptance of proton and antiproton respectively.)
    TFile *f3 = new TFile("EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root");
    TFile *f4 = new TFile("EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root");
    TGraphAsymmErrors *Acceptance_antiproton_all = (TGraphAsymmErrors*)f3->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_proton_all     = (TGraphAsymmErrors*)f4->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TH1D *hAcceptance_antiproton                 = new TH1D("", "", subrangepointused.size(), subrange_low.data()); //GetX()[1]=1.08(1.0-1.16)
    TH1D *hAcceptance_proton                     = new TH1D("", "", subrangepointused.size(), subrange_low.data()); //GetX()[1]=1.08(1.0-1.16)


    //// Load TriggerEfficiency ////FIX ME: Load proton TriggerEff
    TFile *f5 = new TFile("TriggerEff_B1042_antipr.pl1.1800_7.6_all.root");
    TFile *f6 = new TFile("TriggerEff_B1042_antipr.pl1.1800_7.6_all.root");
    TH1F *Trig_antiproton_noprescaling_all = (TH1F*)f5->Get("TriggerEff_noprescaling"); //GetBinLowEdge(2)=1.0
    TH1F *Trig_proton_noprescaling_all     = (TH1F*)f6->Get("TriggerEff_noprescaling"); //GetBinLowEdge(2)=1.0
    TH1D *Trig_antiproton_noprescaling     = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    TH1D *Trig_proton_noprescaling         = new TH1D("", "", subrangepointused.size(), subrange_low.data());


    //// Fill histogram in intermediate range for MeasuringTime, EffectiveAcceptance, TriggerEfficiency
    for (long unsigned int k = 1; k <= subrangepointused.size(); ++k){
        hAcceptance_antiproton      ->SetBinContent(k, Acceptance_antiproton_all->GetY()[k+startindex-1]); 
        hAcceptance_antiproton      ->SetBinError(k  , Acceptance_antiproton_all->GetErrorY(k+startindex-1));
        hAcceptance_proton          ->SetBinContent(k, Acceptance_proton_all->GetY()[k+startindex-1]);
        hAcceptance_proton          ->SetBinError(k  , Acceptance_proton_all->GetErrorY(k+startindex-1));
        hMeasuringTime              ->SetBinContent(k, fIntegratedMeasuringTimeOverCutOff->GetBinContent(k+startindex)); 
        Trig_antiproton_noprescaling->SetBinContent(k, Trig_antiproton_noprescaling_all->GetBinContent(k+startindex)); 
        Trig_proton_noprescaling    ->SetBinContent(k, Trig_proton_noprescaling_all->GetBinContent(k+startindex)); 
    }


    //// Load Unfolding_Matrices
    TFile *f7 = new TFile("/hpcwork/jara0052/sichen/Unfolding_Matrices/low/Unfolding_MatricesTH2D_fill.root"); //## reconstructed x, true y.
    TH2D *Unfolding_Matrices = (TH2D*)f7->Get("Unfolding_Matrices"); //Unfolding_Matrices->GetX(Y)axis()->GetBinLowEdge(1)=1.00; GetBinContent(0,0)=0, GetBinContent(1,1)!=0;


    //// Calculate Effective_Acceptance ratio, Reset Effective_Acceptance error, Parametrization of ratio of effective acceptance
    TH1D *hEffective_Acceptance;
    TGraphErrors *Effective_Acceptance;
    TGraphErrors *g_Effective_Acceptance_all;
    TF1 *fittedfuction1;
    TGraph *gFitFunction;
    tie(hEffective_Acceptance, Effective_Acceptance, g_Effective_Acceptance_all, fittedfuction1, gFitFunction) = CalculateEffectiveAcceptanceRatio_Parametrization(hAcceptance_proton, hAcceptance_antiproton, Acceptance_antiproton_all, Acceptance_proton_all, subrangepointused);


    //// Perform Unfolding
    const int bayesIterations = 2;
    const int verbosity = 0;
    BayesUnfoldingWithCutoff unfolding_antiproton(hantiproton, hMeasuringTime, hAcceptance_antiproton, Trig_antiproton_noprescaling, Unfolding_Matrices, MaxMeasuringTime, bayesIterations, verbosity);
    TH1D* hAntiproton_unfolded = unfolding_antiproton.UnfoldedEventCounts();
    BayesUnfoldingWithCutoff unfolding_proton    (hproton    , hMeasuringTime, hAcceptance_proton    , Trig_proton_noprescaling    , Unfolding_Matrices, MaxMeasuringTime, bayesIterations, verbosity);
    TH1D* hProton_unfolded =unfolding_proton.UnfoldedEventCounts();


    //// Calculate Ratio, Reset Error as Statistical Error, Fill StatisticalRelError
    TGraphErrors *gRatio_unfolded; 
    TGraphErrors *gRatio_Raw;
    TGraphErrors *g_StatisticalRelError;
    TH1D *h_StatisticalRelError;
    tie(gRatio_unfolded, gRatio_Raw, g_StatisticalRelError, h_StatisticalRelError) = CalculateRatio_And_ResetError(subrangepointused, subrange_low, hantiproton, hproton, hAntiproton_unfolded, hProton_unfolded, gFitFunction, g_error);


    /*
    //// Calculate SystematicError
    TGraphErrors *g_SystematicError = new TGraphErrors(h_StatisticalRelError);
    TH1D *h_SystematicError         = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    for (int q = 1; q <= subrangepointused.size(); q++){
        g_SystematicError->SetPoint(     q-1, g_error->GetX()[q-1], abs(g_antiproton_number_uncertainty->GetY()[q-1] - g_antiproton_number->GetY()[q-1])/g_proton_number->GetY()[q-1]);
        g_SystematicError->SetPointError(q-1, 0, 0);
        h_SystematicError->SetBinContent(q  , abs(g_antiproton_number_uncertainty->GetY()[q-1] - g_antiproton_number->GetY()[q-1])/g_proton_number->GetY()[q-1]);
    }    
    */

    //// Calculate SystematicError
    TGraphErrors *g_SystematicError = new TGraphErrors(h_StatisticalRelError);
    TH1D *h_SystematicError         = new TH1D("", "", subrangepointused.size(), subrange_low.data());
    for (int q = 1; q <= subrangepointused.size(); q++){
        g_SystematicError->SetPoint(     q-1, g_error->GetX()[q-1], g_SysError->GetY()[q-1]/100000);
        g_SystematicError->SetPointError(q-1, 0, 0);
        h_SystematicError->SetBinContent(q  , g_SystematicError->GetY()[q-1]);
    }


    //// Plot 
    Plot_unfolded_ratio(gRatio_unfolded, gPublishedRatio, issversion);
    Plot_Effective_Acceptance(Effective_Acceptance);


    //// Save in Root file
    TFile *f = new TFile( (string("Time_Averaged_ratio_Low/plots/unfolded_results") + issversion + string(".root") ).c_str(),"RECREATE");

    // Chi2
    g_chi2_tof                               ->Write("g_chi2_tof");
    // Error
    g_SystematicError                        ->Write("g_SystematicError");
    h_SystematicError                        ->Write("h_SystematicError");
    g_error                                  ->Write("g_StatisticalError");
    g_StatisticalRelError                    ->Write("g_StatisticalRelError");
    h_StatisticalRelError                    ->Write("h_StatisticalRelError");
    PublishedError                           ->Write("g_PublishedError");
    PublishedRatioStatisticErrorPRL          ->Write("g_PublishedRatioStatisticErrorPRL");
    PublishedRatioSystematicErrorPRL         ->Write("g_PublishedRatioSystematicErrorPRL");
    PublishedRatioStatisticRelativeErrorPRL  ->Write("g_PublishedRatioStatisticRelativeErrorPRL");
    PublishedRatioSystematicRelativeErrorPRL ->Write("g_PublishedRatioSystematicRelativeErrorPRL");
    PublishedRatioRelativeErrorPRL           ->Write("g_PublishedRatioRelativeErrorPRL");
    PhysicsReportRatioError                  ->Write("g_PhysicsReportRatioError");
    PhysicsReportRatioStatisticError         ->Write("g_PhysicsReportRatioStatisticError");
    PhysicsReportRatioSystematicError        ->Write("g_PhysicsReportRatioSystematicError");
    PhysicsReportRatioRelativeError          ->Write("g_PhysicsReportRatioRelativeError");
    PhysicsReportStatisticRatioRelativeError ->Write("g_PhysicsReportStatisticRatioRelativeError");
    PhysicsReportSystematicRatioRelativeError->Write("g_PhysicsReportSystematicRatioRelativeError");
    // Number
    g_antiproton_number                      ->Write("g_antiproton_number_raw");
    g_proton_number                          ->Write("g_proton_number_raw");
    hAntiproton_unfolded                     ->Write("hAntiproton_number_unfolded");
    hProton_unfolded                         ->Write("hProton_number_unfolded");
    // Acceptance
    hAcceptance_antiproton                   ->Write("hAcceptance_antiproton");
    hAcceptance_proton                       ->Write("hAcceptance_proton");
    hEffective_Acceptance                    ->Write("hEffective_Acceptance");
    Effective_Acceptance                     ->Write("g_Effective_Acceptance");
    g_Effective_Acceptance_all               ->Write("g_Effective_Acceptance_all");
    fittedfuction1                           ->Write("fittedfuction1");
    gFitFunction                             ->Write("gFitFunction");
    // Ratio
    gRatio_unfolded                          ->Write("gRatio_unfolded");
    gRatio_Raw                               ->Write("gRatio_Raw");
    gPublishedRatio                          .Write("gPublishedRatio");

    f->Close();

}

