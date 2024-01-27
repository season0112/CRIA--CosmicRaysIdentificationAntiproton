#include "AntiprotonIntermediateEnergyTree.hh"
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
#include <TFile.h>
#include <TMath.h>
#include <TCut.h>
#include <TChain.h>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include "Utilities.hh"
#include <cstdlib>

#include "TemplateFitterforIntermediateEnergy.hh"
#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>

#include "AntiprotonAnalysisTools.hh"
#include "AntiprotonBinning.hh"

#include "main.hh"
#include "CutDefinition.hh"

using namespace std;

#define INFO_OUT_TAG "Antiproton_IntermediateTF1D"
#include "debugging.hh"


int main(int argc, char* argv[]) {

Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
config.ReadCommandLine(argc, argv);

std::string issversion = "";
config.GetValue("OPTIONS", "issversion", issversion,
              "The Issversion is");

std::string rigidity_start = "10";
config.GetValue("OPTIONS", "rigidity_start", rigidity_start,
              "The rigidity_start_index is");

std::string rigidity_end = "28";
config.GetValue("OPTIONS", "rigidity_end", rigidity_end,
              "The rigidity_end_index is");

std::string trackerpattern = "";
config.GetValue("OPTIONS", "trackerpattern", trackerpattern,
              "The trackerpattern is");

std::string richcut = "";
config.GetValue("OPTIONS", "richcut", richcut,
              "The richcut is");

std::string signalefficiency = "";
config.GetValue("OPTIONS", "signalefficiency", signalefficiency,
              "whose the signalefficiency then calculated RichBeta cut value is used.");

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
              "If you want to use FullRange of the TRDLikelihood");

if (issversion == "") {
  WARN_OUT << "No issversion is given! Please give a issversion." << std::endl;
  return EXIT_FAIL_CONFIG;
}

if (trackerpattern == "") {
  WARN_OUT << "No trackerpattern is given! Please give a trackerpattern." << std::endl;
  return EXIT_FAIL_CONFIG;
}


string lowpath = (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/")).c_str();
chdir( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/")).c_str());

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
vector<double> subbincenter(bincenter.begin()+stoi(rigidity_start), bincenter.begin()+stoi(rigidity_end));
vector<double> subbinedge  (binning.begin()+stoi(rigidity_start)  , binning.begin()+stoi(rigidity_end)+1);
std::vector<double> RigidityBinPoint_Published( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.end());
vector<double> subbincenter_binmerge;
assert(subbincenter_binmerge.empty());
for(long unsigned int p=0; p<subbincenter.size(); p=p+stoi(binmerge)) {
    subbincenter_binmerge.push_back(subbincenter.at(p));}
vector<double> subbinedge_binmerge;
assert(subbinedge_binmerge.empty());
for(long unsigned int p=0; p<subbinedge.size(); p=p+stoi(binmerge)) {
    subbinedge_binmerge.push_back(subbinedge.at(p));}

vector<double> v_pattern0_percentage, v_pattern1_percentage, v_pattern2_percentage, v_pattern4_percentage;
vector<double> result, ResultError, v_ratio_result, v_ratio_result_with_effective_acceptance, v_chi2dof, v_antiproton, v_proton, v_error, v_error_relative;
vector<double> v_StatisticError_relative(AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.begin()+stoi(rigidity_start), AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.begin()+stoi(rigidity_end));
vector<double> v_SystematicError_relative(AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.begin()+stoi(rigidity_start), AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.begin()+stoi(rigidity_end));
TGraphErrors g_PRLResult = TGraphErrors(57, RigidityBinPoint_Published.data(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.data(), 0, AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.data());


//// Load
// Load Published Result Error
TH1D h_Published_Statistic_Error_Relative;
TH1D h_Published_Systematic_Error_Relative;
tie(h_Published_Statistic_Error_Relative, h_Published_Systematic_Error_Relative) = LoadPublishedResultError(subbincenter, v_StatisticError_relative, subbinedge, v_SystematicError_relative);
// Load Effective acceptance
TGraphErrors *Effective_Acceptance;
tie(Effective_Acceptance) = LoadEffectiveAcceptance(lowpath);


//// Set TRD Efficiency for Systematic Error Uncertainty Calculation
vector<string> TRDEffAll;
double TRDEffStart = 0.70;
double TRDEffEnd   = 0.96;
double TRDEffStep  = 0.02;
/*
double TRDEffStart = 0.60;
double TRDEffEnd   = 1.00;
double TRDEffStep  = 0.01;
*/
for (double i=TRDEffStart; i<=TRDEffEnd; i=i+TRDEffStep){
    TRDEffAll.push_back( to_string_with_precision(i,2) );
}


//// Open ROOT File to Save result
std::string ResultRootName;
if (FullRange == "Yes"){
    ResultRootName = "/intermediate_fullrange_";;
    TRDEffAll.clear();
    TRDEffAll.push_back( to_string_with_precision(1.00,2) );
}
else if (FullRange == "No"){
    ResultRootName = "/intermediate_";
}
TFile ResultROOT(( string("Time_Averaged_ratio/binmerge") + binmerge + ResultRootName + trackerpattern + "_" + richcut + "_" + issversion + "binmerge" + binmerge + suffix + string(".root") ).c_str(),"RECREATE");


//// Loop in TRD Efficiency
for (auto &trdeff: TRDEffAll){  // Trd Efficiency loop start

    //// Loop in generate number template loop
    for (int numberindex = 0; numberindex < GenerateNumbers; numberindex=numberindex+1){ // generate number loop start
        cout<< "numberindex is: " << numberindex << endl;

        //// Loop in Rigidity
        for (int index = stoi(rigidity_start); index < stoi(rigidity_end); index = index + stoi(binmerge)){   // rigidity loop start
            cout<< "\n" << endl;
            cout<< "Current Index is " << index << endl;
            cout<< "Current Rigidity Bining is " << binning[index] << endl;

            //// Load Histograms
            TH1F *data_pass7_positive    = new TH1F();
            TH1F *template_electron_Data = new TH1F();
            TH1F *template_pion_Data     = new TH1F();
            TH1F *data_pass7_negative    = new TH1F();
            double ProtonNumber;
            tie(data_pass7_positive, template_electron_Data, template_pion_Data, data_pass7_negative, ProtonNumber) = LoadTemplates(lowpath, binmerge, binning, index, issversion, ParametrilizedMode, RescaleFactor, numberindex, trdeff, FullRange);


            //// Plot Tracker Pattern Percentage 
            //Plot_Pattern_Percentage(richcut, issversion, binmerge, v_pattern0_percentage, v_pattern1_percentage, v_pattern2_percentage, v_pattern4_percentage, subbincenter, fpass7_negative, RichBetaCut, pattern0_percentage, pattern1_percentage, pattern2_percentage, pattern4_percentage, AllCut);

            //// 1D Template Fit
            MYUtilities::TemplateFitter templateFitter(-1); ////printlevel 0:standard, -1:quiet
            templateFitter.AddTemplateHistogram(data_pass7_positive);
            templateFitter.AddTemplateHistogram(template_electron_Data);
            //templateFitter.AddTemplateHistogram(template_pion_Data);
            templateFitter.SetDataHistogram(data_pass7_negative);

            double fitrange = -0.6805105;
            if (FullRange == "Yes"){
                //templateFitter.Fit(1, -2.0, 0);}
                templateFitter.Fit(1, -2.0, fitrange);} //1D Template Fit support fit ranges. First index:0->chisqr_fit, 1->likelihood_fit
            else{
                templateFitter.Fit(1, -2.0, 0);}
            /* 
            if (index == 27 || index== 28 || index== 29){
                templateFitter.AddTemplateHistogram(data_pass7_positive);
                templateFitter.AddTemplateHistogram(template_electron_Data);
                templateFitter.SetDataHistogram(data_pass7_negative);
                templateFitter.Fit(1, -2.0,0);
            }
            else{
                templateFitter.AddTemplateHistogram(data_pass7_positive);
                templateFitter.AddTemplateHistogram(template_electron_Data);
                templateFitter.AddTemplateHistogram(template_pion_Data);
                //templateFitter.AddTemplateHistogram(template_electron);
                //templateFitter.AddTemplateHistogram(template_pion);

                templateFitter.SetDataHistogram(data_pass7_negative);
                //templateFitter.SetStartValue(0, 24);
                //templateFitter.SetStartValue(1, 27);
                templateFitter.Fit(1, -2.0, 0);
            }
            */


            //// Fill Template Fit Result into Vector (Proton numbers dependes on different cases)
            // Get Fit Result
            assert(result.empty());
            assert(ResultError.empty());
            result     .assign(templateFitter.GetAbsoluteResult().begin()     , templateFitter.GetAbsoluteResult().end());
            ResultError.assign(templateFitter.GetAbsoluteResultError().begin(), templateFitter.GetAbsoluteResultError().end());
            double Chi2    = templateFitter.Chi2(); 
            int NDF        = templateFitter.NDF();
            double CHI2dof = Chi2/NDF;
            std::cout<< Chi2 << std::endl;
            std::cout<< NDF << std::endl;

            // Effective Acceptance Correction (handle with binmerge=1 or 2)
            double effective_all=0, effective_mean=0;
            for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                effective_all = effective_all+ Effective_Acceptance->GetY()[index+binindex];
            }
            effective_mean = effective_all/stod(binmerge);

            // Fill into vector
            v_ratio_result.push_back(result.at(0)/ProtonNumber);
            v_ratio_result_with_effective_acceptance.push_back(result.at(0)/ProtonNumber*effective_mean );
            v_chi2dof.push_back(CHI2dof);
            v_antiproton.push_back(result.at(0));
            v_proton.push_back(ProtonNumber);
            v_error.push_back(ResultError[0]/ProtonNumber);
            v_error_relative.push_back( (ResultError[0]/ProtonNumber) / (result.at(0)/ProtonNumber*effective_mean) *100);


            //// Plot for only original template fit
            if (ParametrilizedMode == "No" && (trdeff == "0.90" || trdeff == "1.00") ){
                // Print Result 
                PrintFitResult(result, ResultError, CHI2dof, ProtonNumber);

                // Plot Templates and Data
                Plot_ISSPositiveData       (data_pass7_positive   , trackerpattern, richcut, issversion, binning, index, binmerge, trdeff);
                Plot_ElectronDataTemplate  (template_electron_Data, trackerpattern, richcut, issversion, binning, index, binmerge, trdeff);
                Plot_PionDataTemplate      (template_pion_Data    , trackerpattern, richcut, issversion, binning, index, binmerge, trdeff);
                Plot_ISSNegativeData       (data_pass7_negative   , trackerpattern, richcut, issversion, binning, index, binmerge, trdeff);

                // Plot Fit Result and Projections
                templateFitter.CreateResultDrawing("Fit_Result",800,600)          ->SaveAs((string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/FitResult_")      + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + suffix + string(".pdf")).c_str());
                templateFitter.CreateResultDrawing_LogY("Fit_Result_LogY",800,600, fitrange)->SaveAs((string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/FitResult_LogY_") + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + suffix + string(".pdf")).c_str());
                templateFitter.CreateContourPlot(0,1,"ContourPlot",800,600)       ->SaveAs((string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/ContourPlot_")    + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + suffix + string(".pdf")).c_str());

            }


            //// Clear Fit Result Vectors
            result.clear();
            ResultError.clear();     


        }   // rigidity loop end


        //// Make TGraph and TH1D to hold ratio result
        cout<< v_ratio_result_with_effective_acceptance.at(0) << endl;
        TGraph *g_ratio                      = new TGraph(v_ratio_result.size()      , subbincenter_binmerge.data(), v_ratio_result.data());
        TGraphErrors *g_ratio_with_effective = new TGraphErrors(v_ratio_result.size(), subbincenter_binmerge.data(), v_ratio_result_with_effective_acceptance.data(), 0, v_error.data());
        TGraph *g_chi2dof                    = new TGraph(v_chi2dof.size()           , subbincenter_binmerge.data(), v_chi2dof.data());
        TGraph *g_antiproton                 = new TGraph(v_antiproton.size()        , subbincenter_binmerge.data(), v_antiproton.data());
        TGraph *g_proton                     = new TGraph(v_proton.size()            , subbincenter_binmerge.data(), v_proton.data());
        TGraph *g_error                      = new TGraph(v_error.size()             , subbincenter_binmerge.data(), v_error.data());
        TGraph *g_error_relative             = new TGraph(v_error_relative.size()    , subbincenter_binmerge.data(), v_error_relative.data());
        TH1D h_error_relative = TH1D("", "", v_error_relative.size(), subbinedge_binmerge.data());
        Utilities::ConvertToHistogram(g_error_relative,h_error_relative);


        //// Plot Ratio only for original template fit
        //if (ParametrilizedMode == "No" and trdeff == "0.90"){
        if (ParametrilizedMode == "No" and FullRange == "No"){
            Plot_ratio(g_PRLResult, g_ratio, trackerpattern, richcut, issversion, binmerge, trdeff);
            Plot_ratio_with_effective_acceptance(g_PRLResult, g_ratio_with_effective, issversion, trackerpattern, richcut, binmerge, trdeff);
            Plot_Error_Compare(h_error_relative, h_Published_Statistic_Error_Relative, h_Published_Systematic_Error_Relative, trackerpattern, richcut, issversion, binmerge, trdeff);
            Plot_Chi2(g_chi2dof, trackerpattern, richcut, issversion, binmerge, trdeff);
        }    


        //// Save Fit Result in ROOT File
        ResultROOT.cd();
        if (ParametrilizedMode == "No" and FullRange == "No"){ 
            Effective_Acceptance                ->Write( (std::string("Effective_Acceptance")                  + std::string("_TRDeff_") + trdeff).c_str() );
            g_chi2dof                           ->Write( (std::string("g_chi2dof")                             + std::string("_TRDeff_") + trdeff).c_str() );
            g_antiproton                        ->Write( (std::string("g_antiproton")                          + std::string("_TRDeff_") + trdeff).c_str() );
            g_proton                            ->Write( (std::string("g_proton")                              + std::string("_TRDeff_") + trdeff).c_str() );
            g_ratio                             ->Write( (std::string("g_ratio")                               + std::string("_TRDeff_") + trdeff).c_str() );
            g_ratio_with_effective              ->Write( (std::string("g_ratio_with_effective_acceptance")     + std::string("_TRDeff_") + trdeff).c_str() );
            g_PRLResult                          .Write( (std::string("publisehd_ratio")                       + std::string("_TRDeff_") + trdeff).c_str() );    
            h_Published_Statistic_Error_Relative .Write( (std::string("h_Published_Statistic_Error_Relative")  + std::string("_TRDeff_") + trdeff).c_str() );
            h_Published_Systematic_Error_Relative.Write( (std::string("h_Published_Systematic_Error_Relative") + std::string("_TRDeff_") + trdeff).c_str() );
            g_error                             ->Write( (std::string("g_error")                               + std::string("_TRDeff_") + trdeff).c_str() );
            g_error_relative                    ->Write( (std::string("g_error_relative")                      + std::string("_TRDeff_") + trdeff).c_str() );
            h_error_relative                     .Write( (std::string("h_error_relative")                      + std::string("_TRDeff_") + trdeff).c_str() );
        }
        else if (ParametrilizedMode == "Yes"){
            cout<< "save result:" << g_ratio_with_effective->GetY()[0] << endl;
            g_ratio_with_effective           ->Write( (std::string("g_ratio_with_effective_acceptance") + to_string(numberindex)).c_str() );
        }


        //// Clear Ratio Vectors
        v_ratio_result.clear();
        v_ratio_result_with_effective_acceptance.clear();
        v_chi2dof.clear();
        v_antiproton.clear();
        v_proton.clear();
        v_error.clear();
        v_error_relative.clear();


    } // generate number loop end

} // Trd Efficiency loop end 

ResultROOT.Close();

return EXIT_SUCCESS;
}
