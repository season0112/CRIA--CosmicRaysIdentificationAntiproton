#include "AntiprotonLowEnergyTree.hh"
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
#include <iostream>
#include <string>
#include <cassert>
#include <TH2.h>
#include <TCut.h>
#include <TFile.h>
#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <TCanvas.h>
#include <TF1.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include "Utilities.hh"
#include <TChain.h>

#include "AntiprotonBinning.hh"
#include "TemplateFitterforLowEnergy2D.hh"

#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>
#include <regex> 
#include <iomanip>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "AntiprotonAnalysisTools.hh"

#include "CommonHeader.hh"

using namespace std;

#define INFO_OUT_TAG "Antiproton_LowTF1D"
#include "debugging.hh"


int main(int argc, char* argv[]) {

Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
config.ReadCommandLine(argc, argv);

config.SetProgramHelpText("Antiproton_LowTF1D",
                            "Template fit for antiproton signal determination in Low Rigidity Range.");
config.AddHelpExample("Antiproton_LowTF1D", "--issversion pass7.8 --rigidity_start 4 --rigidity_end 6 --binmerge 1 --Testmode Yes --ParametrilizedMode No --GenerateNumbers 1");

std::string issversion = "";
config.GetValue("OPTIONS", "issversion", issversion,
              "The Issversion is");

std::string rigidity_start = "";
config.GetValue("OPTIONS", "rigidity_start", rigidity_start,
              "The rigidity_start_index is");

std::string rigidity_end = "";
config.GetValue("OPTIONS", "rigidity_end", rigidity_end,
              "The rigidity_end_index is");

std::string binmerge = "";
config.GetValue("OPTIONS", "binmerge", binmerge,
              "how many bins used for each fit");

std::string TestMode = "";
config.GetValue("OPTIONS", "Testmode", TestMode,
              "If you want to use TestMode to run program quickly");

std::string ParametrilizedMode = "";
config.GetValue("OPTIONS", "ParametrilizedMode", ParametrilizedMode,
              "Parametrilized Template or not");

double GenerateNumbers = -999;
config.GetValue("OPTIONS", "GenerateNumbers", GenerateNumbers,
              "GenerateNumbers is");

std::string FullRange = "";
config.GetValue("OPTIONS", "FullRange", FullRange,
              "If you want to use FullRange of the fit range");

if (issversion == "") {
  WARN_OUT << "No issversion is given! Please give a issversion." << std::endl;
  return EXIT_FAIL_CONFIG;
}

if (rigidity_start == "") {
  WARN_OUT << "No rigidity_start_index is given! Please give a rigidity_start_index." << std::endl;
  return EXIT_FAIL_CONFIG;
}

if (rigidity_end == "") {
  WARN_OUT << "No rigidity_end_index is given! Please give a rigidity_end_index." << std::endl;
  return EXIT_FAIL_CONFIG;
}


chdir( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/")).c_str());
string lowpath = (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/")).c_str();

std::string suffix;
if (ParametrilizedMode == "No"){
    suffix = "";}
else if (ParametrilizedMode == "Yes"){
    suffix = "_uncertainty";}
double RescaleFactor;
if (TestMode == "Yes"){
    RescaleFactor = 20;}
else{
    RescaleFactor = 1;}


//// Prepare: binning defined, result vector defined
vector<double> binning     (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());
vector<double> bincenter   (AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins());
vector<double> subbinedge  (binning.begin()+stoi(rigidity_start)  , binning.begin()+stoi(rigidity_end)+1+(stoi(binmerge)-1));
vector<double> subbincenter(bincenter.begin()+stoi(rigidity_start), bincenter.begin()+stoi(rigidity_end)+(stoi(binmerge)-1));

vector<double> subbincenter_binmerge;
vector<double> subbinedge_binmerge;

for(long unsigned int p=stoi(binmerge)-1; p<subbincenter.size(); p=p+stoi(binmerge)) {  // for binmerge2: from p=1 (1.61:1.51-1.71) because 1.33-1.51 is not included. for binmerge1: from 1.33-1.51
    cout<< p <<endl;
    subbincenter_binmerge.push_back( (subbincenter.at(p)+subbincenter.at(p+stoi(binmerge)-1))/2 );}

for(long unsigned int p=stoi(binmerge)-1; p<subbinedge.size(); p=p+stoi(binmerge)) {
    cout<< p <<endl;
    subbinedge_binmerge.push_back( subbinedge.at(p) );}

vector<double> result_tof, ResultError_tof, v_ratio_result_tof, v_chi2_tof, v_ratio_result_tof_with_effective, v_error, v_error_relative, v_antiproton_number, v_proton_number;


//// Load
// Load Published result
TGraphErrors *g_ratio_published;
TH1D h_Published_Statistic_Error_Relative;
TH1D h_Published_Systematic_Error_Relative;
tie(g_ratio_published, h_Published_Statistic_Error_Relative, h_Published_Systematic_Error_Relative) = LoadPublishedResult(bincenter, rigidity_start, rigidity_end, subbincenter, subbinedge);

// Load Effective acceptance
TGraphErrors *Effective_Acceptance;
int EffectiveAcceptanceRange = 17; // Number17: 0.8-5.9GV 
tie(Effective_Acceptance) = LoadEffectiveAcceptance(EffectiveAcceptanceRange);


//// Set TRD Efficiency for Systematic Error Uncertainty Calculation
vector<string> TRDEffAll;
//double TRDEffStart = 0.70;
//double TRDEffEnd   = 0.96;
//double TRDEffStep  = 0.02;

double TRDEffStart = 0.60;
//double TRDEffStart = 0.92;
double TRDEffEnd   = 1.00;
double TRDEffStep  = 0.02;
for (double i=TRDEffStart; i<=TRDEffEnd; i=i+TRDEffStep){
    TRDEffAll.push_back( to_string_with_precision(i,2) );
}

//// Set TOF Beta Efficiency for Systematic Error Uncertainty Calculation
vector<string> TOFEffAll;
//double TOFEffStart = 0.70;
//double TOFEffEnd   = 0.96;
//double TOFEffStep  = 0.02;

double TOFEffStart = 0.61;
//double TOFEffStart = 0.93;
double TOFEffEnd   = 0.99;
double TOFEffStep  = 0.02;
for (double i=TOFEffStart; i<=TOFEffEnd; i=i+TOFEffStep){
    TOFEffAll.push_back( to_string_with_precision(i,2) );
}


//// Open ROOT File to Save result
std::string ResultRootName;
if (FullRange == "Yes"){
    ResultRootName = "Ratio_fullrange_";;
    TRDEffAll.clear();
    TOFEffAll.clear();
    TRDEffAll.push_back( to_string_with_precision(1.00,2) );
    TOFEffAll.push_back( to_string_with_precision(1.00,2) );
}
else if (FullRange == "No"){
    ResultRootName = "Ratio";
}
TFile ResultROOT((string("Time_Averaged_ratio_Low/binmerge") + binmerge + string("/plots/") + ResultRootName + string("_") + issversion + suffix + string(".root")).c_str(),"RECREATE");


//// Loop in TRD Efficiency
for (auto &trdeff: TRDEffAll){  // Trd Efficiency loop start

    //// Loop in TOF Efficiency
    for (auto &tofeff: TOFEffAll){  // Tof Efficiency loop start

        //// Loop in generate number template loop
        for (int numberindex = 0; numberindex < GenerateNumbers; numberindex=numberindex+1){ // generate number loop start

            //// Loop in Individual Rigidity Loop
            for (int index = stoi(rigidity_start)+(stoi(binmerge)-1) ; index<stoi(rigidity_end)+(stoi(binmerge)-1); index=index+stoi(binmerge)) {   // rigidity loop start

                //cout<< "\n" <<endl;
                //cout<< string("Index now is:") + to_string(index) <<endl;


                //// Load Histograms   
                TH2F *TofTRD_data_pass7_positive_template  = new TH2F();
                TH2F *TofTRD_data_pass7_positive_data      = new TH2F();
                TH2F *TofTRD_data_pass7_positive_data_test = new TH2F();
                TH2F *TofTRD_data_pass7_negative           = new TH2F();
                TH2F *TofTRD_template_ElectronData         = new TH2F();
                TH2F *TofTRD_template_PionData             = new TH2F();
                double ProtonNumber;
                if (FullRange == "No"){
                    tie(TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data, TofTRD_data_pass7_positive_data_test, TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData, ProtonNumber) = LoadTemplates(binmerge, binning, index, issversion, ParametrilizedMode, numberindex, trdeff, tofeff, FullRange);
                }
                else if (FullRange == "Yes"){
                    double RECITOFBETALOW    = -0.3;
                    double RECITOFBETAHIGH   = 0.2;
                    double RECITOFBETANUMBER = 40;
                    double TrdLOW            = 0.70;
                    double TrdHIGH           = 1.6;
                    double TrdBINNUMBER      = 20;
                    tie(TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data, TofTRD_data_pass7_positive_data_test, TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData, ProtonNumber) = LoadTemplates(binmerge, binning, index, issversion, ParametrilizedMode, numberindex, trdeff, tofeff, FullRange, RECITOFBETALOW, RECITOFBETAHIGH, RECITOFBETANUMBER, TrdLOW, TrdHIGH, TrdBINNUMBER);
                }
                TofTRD_data_pass7_positive_template->SetTitle("Antiprotons");
                TofTRD_template_ElectronData       ->SetTitle("Electrons");
                TofTRD_template_PionData           ->SetTitle("Secondaries"); 

                /* No smooth: othervise the chi2 get worse.
                if (FullRange == "Yes"){
                    TofTRD_template_ElectronData->Smooth(1);
                    TofTRD_template_PionData->Smooth(1);
                }
                */
 
                //// 2D Template Fit
                //cout<< "Begin Template Fit now !"<<endl;
                MYUtilities::TemplateFitter2D templateFitter_toftrd_2D(-1);  //printlevel 0:standard, -1:quiet
                // Add Signal
                templateFitter_toftrd_2D.AddTemplateHistogram(TofTRD_data_pass7_positive_template);
                // Add Background
                templateFitter_toftrd_2D.AddTemplateHistogram(TofTRD_template_ElectronData);
                templateFitter_toftrd_2D.AddTemplateHistogram(TofTRD_template_PionData);
                // Add Data
                templateFitter_toftrd_2D.SetDataHistogram(TofTRD_data_pass7_negative);
                // Perfer Fit
                templateFitter_toftrd_2D.Fit(1,0,0);           //2D Template Fit do not support fit ranges. First index:0->chisqr_fit, 1->likelihood_fit
                

                //// Fill Template Fit Result into Vector (Proton numbers dependes on different cases)
                // Get Fit Result
                assert(result_tof.empty());
                assert(ResultError_tof.empty());
                result_tof     .assign(templateFitter_toftrd_2D.GetAbsoluteResult().begin()     , templateFitter_toftrd_2D.GetAbsoluteResult().end());
                ResultError_tof.assign(templateFitter_toftrd_2D.GetAbsoluteResultError().begin(), templateFitter_toftrd_2D.GetAbsoluteResultError().end());
                double Chi2_tof    = templateFitter_toftrd_2D.Chi2();
                int NDF_tof        = templateFitter_toftrd_2D.NDF();
                double CHI2dof_tof = Chi2_tof/NDF_tof;    

                // Effective Acceptance Correction (handle with binmerge=1 or 2)
                double effective_all=0, effective_mean=0;
                for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                    effective_all = effective_all + Effective_Acceptance->GetY()[index+binindex];
                }
                effective_mean = effective_all/stod(binmerge);
                
                // Fill into vector
                v_ratio_result_tof.push_back(result_tof[0] / (TofTRD_data_pass7_positive_data->GetEntries()*RescaleFactor) *100000); // proton from total number after cut
                //v_ratio_result_tof.push_back(result_tof[0]/result_positive[0]*100000);                                             //proton from template fit
                //v_ratio_result_tof.push_back(result_tof[0]/data_pass7_positive->GetEntries()*100000);                              // tof fit, proton from total number before cut
                v_antiproton_number              .push_back(result_tof[0]);
                v_proton_number                  .push_back(TofTRD_data_pass7_positive_data->GetEntries() * RescaleFactor);
                v_ratio_result_tof_with_effective.push_back(v_ratio_result_tof.back() * effective_mean);

                v_chi2_tof                       .push_back(CHI2dof_tof);
                v_error                          .push_back(ResultError_tof[0] / (TofTRD_data_pass7_positive_data->GetEntries()*RescaleFactor) * 100000);
                v_error_relative                 .push_back( (ResultError_tof[0] / (TofTRD_data_pass7_positive_data->GetEntries()*RescaleFactor) ) / (result_tof.at(0) / (TofTRD_data_pass7_positive_data->GetEntries()*RescaleFactor) * effective_mean ) *100);
               
                //// Plot for only original template fit
                if ( (ParametrilizedMode == "No" and trdeff == "0.94" and tofeff == "0.95") or FullRange == "Yes"){
                    
                    // Print fit result
                    PrintFitResult(result_tof, TofTRD_data_pass7_positive_data, RescaleFactor, TofTRD_data_pass7_positive_data_test, CHI2dof_tof, Chi2_tof, NDF_tof);

                    
                    // Plot Template and data
                    string AnalysisMode = "TimeAveraged";
                    int TimeIndex       = 99999;
                    Plot_TemplatePion_Data      (binmerge, binning, index, issversion, TofTRD_template_PionData           , ParametrilizedMode, trdeff, tofeff, AnalysisMode, TimeIndex);
                    Plot_TemplateElectron_Data  (binmerge, binning, index, issversion, TofTRD_template_ElectronData       , ParametrilizedMode, trdeff, tofeff, AnalysisMode, TimeIndex);
                    Plot_ISSNegative            (binmerge, binning, index, issversion, TofTRD_data_pass7_negative         , ParametrilizedMode, trdeff, tofeff, AnalysisMode, TimeIndex);
                    Plot_TemplateAntiproton_Data(binmerge, binning, index, issversion, TofTRD_data_pass7_positive_template, ParametrilizedMode, trdeff, tofeff, AnalysisMode, TimeIndex);   
                     

                    // Plot Fit Result and Projections 
                    if (templateFitter_toftrd_2D.HasBeenFit() == true){
                        templateFitter_toftrd_2D.CreateResultDrawing                ("Fit_Result_tof",1000,500)->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/FitResult_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + suffix + std::string(".pdf")).c_str());  

                        double OneOverTOFCut_LowEdge_InTRDProjection  = -0.1; //-0.03
                        double OneOverTOFCut_HighEdge_InTRDProjection = 0.3;
                        int BinRemovedNumber_FromLeft = 0;
                        int RebinNumber_X = 1;
                        templateFitter_toftrd_2D.CreateResultDrawingXprojection     ("TrdLikelihood_X",1000,800, OneOverTOFCut_LowEdge_InTRDProjection, OneOverTOFCut_HighEdge_InTRDProjection, BinRemovedNumber_FromLeft, RebinNumber_X)->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TrdLikelihood_X_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + suffix + std::string(".pdf")).c_str());   // ("TrdLikelihood_X",1000,800, RECITOFBETALOW, RECITOFBETAHIGH)->SaveAs(  //RECITOFBETALOW=-0.03: One bin to show for 2.2-2.4GV Generam meeting 06.2021

                        templateFitter_toftrd_2D.CreateResultDrawingXprojection_LogY("TrdLikelihood_X_LogY",1000,800, OneOverTOFCut_LowEdge_InTRDProjection, OneOverTOFCut_HighEdge_InTRDProjection, BinRemovedNumber_FromLeft, RebinNumber_X)->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TrdLikelihood_X_tof_LogY_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + suffix + std::string(".pdf")).c_str());  // ("TrdLikelihood_X_LogY",1000,800, RECITOFBETALOW, RECITOFBETAHIGH)->SaveAs(  // RECITOFBETALOW=-0.03: One bin to show for 2.2-2.4GV Generam meeting 06.2021

                        double TRDCut_LowEdge_InTOFProjection  = 0.90;   
                        double TRDCut_HighEdge_InTOFProjection = 2.0;
                        int BinRemovedNumber_FromRight = 0;
                        int RebinNumber_Y = 2;
                        templateFitter_toftrd_2D.CreateResultDrawingYprojection     ("1/TofBeta_Y",1000,800, TRDCut_LowEdge_InTOFProjection, TRDCut_HighEdge_InTOFProjection, BinRemovedNumber_FromRight, RebinNumber_Y)->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/1_TofBeta_Y_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + suffix + std::string(".pdf")).c_str());   // ("1/TofBeta_Y",1000,800, TrdLOW, TrdHIGH)->SaveAs( // TrdLOW = 0.92: One bin to show for 2.2-2.4GV Generam meeting 06.2021

                        templateFitter_toftrd_2D.CreateResultDrawingYprojection_LogY("1/TofBeta_Y_LogY",1000,800, TRDCut_LowEdge_InTOFProjection, TRDCut_HighEdge_InTOFProjection, BinRemovedNumber_FromRight, RebinNumber_Y)->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/1_TofBeta_Y_tof_LogY_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + suffix + std::string(".pdf")).c_str()); // ("1/TofBeta_Y_LogY",1000,800, TrdLOW, TrdHIGH)->SaveAs // TrdLOW = 0.92: One bin to show for 2.2-2.4GV Generam meeting 06.2021
                    }
                    
                }


                //// Clear Fit Result Vectors
                result_tof.clear();
                ResultError_tof.clear(); 


            } // rigidity loop end


            //// Make TGraph and TH1D to hold ratio result
            TGraph *g_antiproton_number              = new TGraph( v_antiproton_number.size()     , subbincenter_binmerge.data(), v_antiproton_number.data());
            TGraph *g_proton_number                  = new TGraph( v_proton_number.size()         , subbincenter_binmerge.data(), v_proton_number.data() );
            TGraph *g_ratio_tof                      = new TGraph( v_ratio_result_tof.size()      , subbincenter_binmerge.data(), v_ratio_result_tof.data());
            TGraphErrors *g_ratio_tof_with_effective = new TGraphErrors( v_ratio_result_tof.size(), subbincenter_binmerge.data(), v_ratio_result_tof_with_effective.data(), 0, v_error.data());
            TGraph *g_chi2_tof                       = new TGraph( v_chi2_tof.size()              , subbincenter_binmerge.data(), v_chi2_tof.data());
            TGraph *g_error                          = new TGraph( v_error.size()                 , subbincenter_binmerge.data(), v_error.data());
            TGraph *g_error_relative                 = new TGraph( v_error_relative.size()        , subbincenter_binmerge.data(), v_error_relative.data());
            TH1D h_error_relative = TH1D("", "", v_error_relative.size(), subbinedge_binmerge.data());
            Utilities::ConvertToHistogram(g_error_relative,h_error_relative);

            //// Plot Ratio only for original template fit
            //if (ParametrilizedMode == "No" and trdeff == "0.90" and tofeff == "0.90"){
            if (ParametrilizedMode == "No" and FullRange == "No"){
                Plot_Low_Ratio_Tof                         (g_ratio_published, g_ratio_tof               , binmerge  , issversion, ParametrilizedMode, trdeff, tofeff);
                Plot_Low_Ratio_Tof_WithAcceptanceCorrection(g_ratio_published, g_ratio_tof_with_effective, binmerge  , issversion, ParametrilizedMode, trdeff, tofeff);
                Plot_Chi2                                  (g_chi2_tof       , binmerge                  , issversion, ParametrilizedMode, trdeff, tofeff);
                Plot_error_compare                         (h_error_relative , h_Published_Statistic_Error_Relative, h_Published_Systematic_Error_Relative, binmerge, issversion, ParametrilizedMode, trdeff, tofeff);
                Plot_PbarNumbers                           (g_antiproton_number, binmerge, issversion, ParametrilizedMode, trdeff, tofeff);
            }


            //// Save Fit Result in ROOT File
            ResultROOT.cd();
            if (ParametrilizedMode == "No"){
                Effective_Acceptance      ->Write( (std::string("Effective_Acceptance_")     + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_antiproton_number       ->Write( (std::string("g_antiproton_number_")      + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_proton_number           ->Write( (std::string("g_proton_number_")          + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_ratio_tof               ->Write( (std::string("ratio_tof_")                + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_ratio_tof_with_effective->Write( (std::string("ratio_tof_with_effective_") + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_chi2_tof                ->Write( (std::string("chi2_tof_")                 + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_error                   ->Write( (std::string("g_error_")                  + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_error_relative          ->Write( (std::string("g_error_relative_")         + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                g_ratio_published         ->Write( (std::string("published_")                + std::string("TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
            }
            else if (ParametrilizedMode == "Yes"){
                g_antiproton_number       ->Write( (std::string("g_antiproton_number")       + to_string(numberindex)).c_str() );
                g_ratio_tof_with_effective->Write( (std::string("ratio_tof_with_effective")  + to_string(numberindex)).c_str() );
                g_chi2_tof                ->Write( (std::string("chi2_tof")                  + to_string(numberindex)).c_str() );
            }
            

            //// Clear Ratio Vectors
            v_ratio_result_tof.clear();
            v_ratio_result_tof_with_effective.clear();
            v_chi2_tof.clear();
            v_antiproton_number.clear();
            v_proton_number.clear();
            v_error.clear();
            v_error_relative.clear();


        } // generate number loop end

    } // Tof Efficiency loop end

} // Trd Efficiency loop end

ResultROOT.Close();


return EXIT_SUCCESS;
}
