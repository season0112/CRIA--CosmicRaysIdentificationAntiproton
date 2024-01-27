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

#include "CutDefinition.hh"
#include "main.hh"

using namespace std;

#define INFO_OUT_TAG "Antiproton_IntermediateTF1D_TimeStamp_3and6_BartalRotations"
#include "debugging.hh"


int main(int argc, char* argv[]) {

Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
config.ReadCommandLine(argc, argv);

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

std::string timemode = "";
config.GetValue("OPTIONS", "timemode", timemode,
              "Time mode you choose.");

std::string binmerge = "";
config.GetValue("OPTIONS", "binmerge", binmerge,
              "how many bins used for each fit");

std::string ParametrilizedMode = "";
config.GetValue("OPTIONS", "ParametrilizedMode", ParametrilizedMode,
              "Parametrilized Template or not");

double GenerateNumbers = -999;
config.GetValue("OPTIONS", "GenerateNumbers", GenerateNumbers,
              "GenerateNumbers is");

std::string SysErrEffStudyMode = "";
config.GetValue("OPTIONS", "SysErrEffStudyMode", SysErrEffStudyMode,
              "SysErrEffStudyMode for SysErr from Signal efficiency");

std::string StartTimeIndex = "";
config.GetValue("OPTIONS", "StartTimeIndex", StartTimeIndex,
              "If you need to start from a certain TimeIndex");

if (trackerpattern == "") {
  WARN_OUT << "No trackerpattern is given! Please give a trackerpattern." << std::endl;
  return EXIT_FAIL_CONFIG;
}

if (richcut == "") {
  WARN_OUT << "No richcut is given! Please give a richcut." << std::endl;
  return EXIT_FAIL_CONFIG;
}

std::string FullRange = "";
config.GetValue("OPTIONS", "FullRange", FullRange,
              "If you want to use FullRange of the TRDLikelihood");

string lowpath = (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/")).c_str();
chdir( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/")).c_str());


//// Define binnings
vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());
vector<double> bincenter (AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins());
vector<double> subbincenter(bincenter.begin()+stoi(rigidity_start), bincenter.begin()+stoi(rigidity_end));
vector<double> subbinedge(binning.begin()+stoi(rigidity_start), binning.begin()+stoi(rigidity_end)+1);
vector<double> subbincenter_binmerge;
assert(subbincenter_binmerge.empty());
for(long unsigned int p=0; p<subbincenter.size(); p=p+stoi(binmerge)) {
    subbincenter_binmerge.push_back(subbincenter.at(p));}

vector<double> v_ratio_result, v_chi2, v_error, v_antiproton_number, v_proton_number;


//// Load TimeStape for Setings
string bartals = "";
string timeindex = "";
int maxindex;
tie(bartals, timeindex, maxindex) = LoadTimeStapeSetings(timemode);


//// Set TRD Efficiency for Systematic Error Uncertainty Calculation
// For Systematic Error Uncertainty Calculation
vector<string> TRDEffAll;
if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes"){
    /*
    double TRDEffStart = 0.70;
    double TRDEffEnd   = 0.96;
    double TRDEffStep  = 0.02;
    */
    double TRDEffStart = 0.60;
    double TRDEffEnd   = 1.00;
    double TRDEffStep  = 0.01;

    for (double i=TRDEffStart; i<=TRDEffEnd; i=i+TRDEffStep){
        TRDEffAll.push_back( to_string_with_precision(i,2) );
    }
}
// For representative point
if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
    //TRDEffAll.push_back( to_string_with_precision(0.94, 2) );
    //TRDEffAll.push_back( to_string_with_precision(0.98, 2) );
    TRDEffAll.push_back( to_string_with_precision(0.99, 2) );
}

if (FullRange == "Yes"){
    TRDEffAll.clear();
    TRDEffAll.push_back( to_string_with_precision(1.00,2) );
}



//// Create txt files to save result.
if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No" and FullRange == "No"){
    for (int RigidityIndex=stoi(rigidity_start); RigidityIndex < stoi(rigidity_end); RigidityIndex=RigidityIndex+stoi(binmerge)){  
        CreateTxtFiles(timemode, binning, RigidityIndex, binmerge, bartals, trackerpattern, richcut);
    }
}


//// Loop Time bins
for (int i = stoi(StartTimeIndex); i < maxindex; i = i + stoi(bartals)){   //i from 0: first bartel rotation included, i from 1:first bartel rotation removed.  // begin the timestamp loop
    std::cout << "Now Timebins is " << i << std::endl;

    //// Open ROOT File to Save result
    TFile *ResultROOT_SysErr_Shape;
    TFile *ResultROOT_SysErr_Eff;
    if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes" and FullRange == "No"){
        ResultROOT_SysErr_Shape = new TFile ((string("SysError_Shape/") + string("Ratio_") + timemode + string("_TimeIndex") + to_string(i) + string(".root")).c_str(),"RECREATE");
        ResultROOT_SysErr_Eff   = new TFile ((string("SysError_Eff/")   + string("Ratio_") + timemode + string("_TimeIndex") + to_string(i) + string(".root")).c_str(),"RECREATE");
    }

    //// Loop in TRD Efficiency
    for (auto &trdeff: TRDEffAll){  // Trd Efficiency loop start

        //// Loop in generate number template loop
        for (int numberindex = 0; numberindex < GenerateNumbers; numberindex=numberindex+1){ // generate number loop start    

            //// Rigidity Loop
            for ( int rigidityindex = stoi(rigidity_start); rigidityindex< stoi(rigidity_end); rigidityindex=rigidityindex+stoi(binmerge) ){  // Rigidity loop start 

                if (timemode == "3BartalRotation" and i == 123){
                    //// fill vectors forr fit result                
                    v_antiproton_number.push_back(0);
                    v_proton_number    .push_back(0);
                    v_ratio_result     .push_back(0); 
                    v_chi2             .push_back(0);
                    v_error            .push_back(0);

                    vector<double> result       = {0, 0, 0};
                    double         ProtonNumber = 0;
                    vector<double> ResultError  = {0, 0, 0};
                    double         CHI2dof      = 0;

                    //// Save result in txt files.
                    if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
                        SaveResult(bartals, binning, rigidityindex, binmerge, result, ProtonNumber, ResultError, CHI2dof, trackerpattern, richcut);
                    }                
                }

                else{
                    //// Load Histograms
                    TH1F *data_pass7_negative;
                    TH1F *template_electron_Data;
                    TH1F *template_pion_Data;
                    TH1F *data_pass7_positive;
                    tie(data_pass7_negative, template_electron_Data, template_pion_Data, data_pass7_positive) = LoadTemplates_TimeDependent(binmerge, binning, rigidityindex, timemode, i, ParametrilizedMode, numberindex, SysErrEffStudyMode, trdeff, FullRange);
                    data_pass7_positive->SetTitle("Antiproton");
                    template_electron_Data->SetTitle("Electron");

                    //// Template Fit
                    MYUtilities::TemplateFitter templateFitter(-1); // output: 0 standard, -1 quiet
                    // Add templates
                    //templateFitter.AddTemplateHistogram(template_antiproton);
                    templateFitter.AddTemplateHistogram(data_pass7_positive);
                    //templateFitter.AddTemplateHistogram(template_electron);
                    //templateFitter.AddTemplateHistogram(template_pion);
                    templateFitter.AddTemplateHistogram(template_electron_Data);
                    //templateFitter.AddTemplateHistogram(template_pion_Data);
                    // Add data
                    templateFitter.SetDataHistogram(data_pass7_negative);

                    // Fit
                    // double fitrange = -0.6805105;
                    double fitrange = 999;
                    templateFitter.Fit(1);
                    /*
                    if (FullRange == "Yes"){
                        templateFitter.Fit(1, -2.0, fitrange);
                    }
                    else{
                        templateFitter.Fit(1);  // Don't set fit range here, user default.  //chi-square (i==0) or likelihood fit (i==1)(default).
                    }
                    */

                    //// Plots
                    if (ParametrilizedMode == "No" && SysErrEffStudyMode == "No" && (trdeff == "0.90" || trdeff == "1.00")){
                        TCanvas * c1 = templateFitter.CreateResultDrawing("Fit_Result",800,600);
                        c1->SaveAs(string("results/fit_plot/FitResult_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + "_" + doubleToString(binning[rigidityindex+stoi(binmerge)]) + string("_index_") + i + std::string("_TRDEff_") + trdeff + string(".pdf"));
                        TCanvas * c2 = templateFitter.CreateResultDrawing_LogY("Fit_Result",800,600, fitrange);
                        c2->SaveAs(string("results/fit_plot/FitResult_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + "_" + doubleToString(binning[rigidityindex+stoi(binmerge)]) + string("_index_") + i + std::string("_TRDEff_") + trdeff + string("_LogY.pdf"));
                        //TCanvas * c5 = templateFitter.CreateContourPlot(0,1,"ContourPlot",800,600);
                        //c5->SaveAs("ContourPlot.pdf");
                    }

                    //// Template Fit result
                    vector<double> result, ResultError;
                    assert(result.empty());
                    assert(ResultError.empty());
                    result.assign(templateFitter.GetAbsoluteResult().begin(), templateFitter.GetAbsoluteResult().end());
                    ResultError.assign(templateFitter.GetAbsoluteResultError().begin(), templateFitter.GetAbsoluteResultError().end());
                    double Chi2    = templateFitter.Chi2(); 
                    int NDF        = templateFitter.NDF();
                    double CHI2dof = Chi2/NDF;

                    //// fill vectors forr fit result
                    v_antiproton_number.push_back(result[0]);
                    double ProtonNumber = data_pass7_positive->GetEntries();
                    v_proton_number.push_back(ProtonNumber);
                    v_ratio_result.push_back(result[0]/ProtonNumber*100000); // proton from total number before cut
                    v_chi2.push_back(CHI2dof);
                    v_error.push_back(ResultError[0]/ProtonNumber*100000);


                    //// Save result in txt files.
                    if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No" and FullRange == "No"){
                        SaveResult(bartals, binning, rigidityindex, binmerge, result, ProtonNumber, ResultError, CHI2dof, trackerpattern, richcut);
                    }


                    //// Clear
                    result.clear();
                    ResultError.clear();            
                }

            }  // rigidity loop end.


            // Make TGraph and TH1D to hold ratio result
            TGraph *g_ratio                      = new TGraph( v_ratio_result.size()      , subbincenter_binmerge.data(), v_ratio_result.data());


            //// Save Fit Result in ROOT File
            if (ParametrilizedMode == "Yes" and SysErrEffStudyMode == "No" and FullRange == "No"){
                ResultROOT_SysErr_Shape->cd();
                g_ratio->Write( (std::string("g_ratio") + to_string(numberindex)).c_str() );
            }

            if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes" and FullRange == "No"){
                ResultROOT_SysErr_Eff->cd();
                g_ratio->Write( (std::string("g_ratio") + std::string("_TRDeff_") + trdeff).c_str() );
            }

            //// Clear
            v_ratio_result.clear();
            v_chi2.clear();
            v_error.clear();
            v_antiproton_number.clear();
            v_proton_number.clear();
            

        } // generate number loop end

    } // Trd Efficiency loop end

    if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes" and FullRange == "No"){
        ResultROOT_SysErr_Shape->Close();
        ResultROOT_SysErr_Eff  ->Close();
    }
  
} // Time index loop end


    return EXIT_SUCCESS;
}


