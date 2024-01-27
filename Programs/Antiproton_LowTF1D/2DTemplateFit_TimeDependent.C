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

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include "AntiprotonAnalysisTools.hh"

#include "CommonHeader.hh"
#include "TimeDependentHeader.hh"

using namespace std;

#define INFO_OUT_TAG "Antiproton_LowTF1D_TimeStamp_3and6_BartalRotations"
#include "debugging.hh"


int main(int argc, char* argv[]) {

Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
config.ReadCommandLine(argc, argv);

config.SetProgramHelpText("Antiproton_LowTF1D",
                            "Template fit for antiproton signal determination in Low Rigidity Range.");

config.AddHelpExample("Antiproton_LowTF1D", "--timemode 6months --rigidity_start 4 --rigidity_end 6 --binmerge 2 --ParametrilizedMode No --GenerateNumbers 1500");

std::string rigidity_start = "";
config.GetValue("OPTIONS", "rigidity_start", rigidity_start,
              "The rigidity_start_index is");

std::string rigidity_end = "";
config.GetValue("OPTIONS", "rigidity_end", rigidity_end,
              "The rigidity_end_index is");

std::string binmerge = "";
config.GetValue("OPTIONS", "binmerge", binmerge,
              "how many bins used for each fit");

std::string timemode = "";
config.GetValue("OPTIONS", "timemode", timemode,
              "Time mode you choose.");

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

std::string EndTimeIndex = "";
config.GetValue("OPTIONS", "EndTimeIndex",EndTimeIndex,
              "If you need to end from a certain TimeIndex");

if (rigidity_start == "") {
  WARN_OUT << "No rigidity_start_index is given! Please give a rigidity_start_index." << std::endl;
  return EXIT_FAIL_CONFIG;
}

if (rigidity_end == "") {
  WARN_OUT << "No rigidity_end_index is given! Please give a rigidity_end_index." << std::endl;
  return EXIT_FAIL_CONFIG;
}

chdir( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/")).c_str());


//// Define binnings
vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());
vector<double> bincenter (AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins());
vector<double> subbincenter(bincenter.begin()+stoi(rigidity_start), bincenter.begin()+stoi(rigidity_end)+2);
vector<double> subbincenter_binmerge;
for(long unsigned int BinIndex=stoi(binmerge)-1; BinIndex<subbincenter.size()-(stoi(binmerge)-1)+1 ; BinIndex=BinIndex+stoi(binmerge)) {  // for binmerge2: from p=1 (1.61:1.51-1.71) because 1.33-1.51 is not included. for binmerge1: from 1.33-1.51
    subbincenter_binmerge.push_back( (subbincenter.at(BinIndex)+subbincenter.at(BinIndex+stoi(binmerge)-1))/2 );
    cout<< "bins: " << (subbincenter.at(BinIndex)+subbincenter.at(BinIndex+stoi(binmerge)-1))/2 << endl; // check here!
}
vector<double> v_ratio_result_tof, v_chi2_tof, v_error, v_antiproton_number, v_proton_number;


//// Load TimeStape for Setings
string bartals   = "";
string timeindex = "";
int maxindex;
tie(bartals, timeindex, maxindex) = LoadTimeStapeSetings(timemode);
//// Rest Time loop end if EndTimeIndex has value
if (EndTimeIndex != ""){
    maxindex = stoi(EndTimeIndex);
}

//// Create txt files to save result.
if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
    for (int RigidityIndex=stoi(rigidity_start); RigidityIndex < stoi(rigidity_end); RigidityIndex=RigidityIndex+stoi(binmerge)){  
        CreateTxtFiles(timemode, binning, RigidityIndex, binmerge, bartals);
    }
}


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
    double TRDEffStep  = 0.02;

    for (double i=TRDEffStart; i<=TRDEffEnd; i=i+TRDEffStep){
        TRDEffAll.push_back( to_string_with_precision(i,2) );
    }
}
// For representative point
if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
    //TRDEffAll.push_back( to_string_with_precision(0.94,2) );
    TRDEffAll.push_back( to_string_with_precision(0.98,2) );
}


//// Set TOF Beta Efficiency for Systematic Error Uncertainty Calculation
// For Systematic Error Uncertainty Calculation
vector<string> TOFEffAll;
if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes"){
    /*
    double TOFEffStart = 0.70;
    double TOFEffEnd   = 0.96;
    double TOFEffStep  = 0.02;
    */
    double TOFEffStart = 0.61;
    double TOFEffEnd   = 0.99;
    double TOFEffStep  = 0.02;

    for (double i=TOFEffStart; i<=TOFEffEnd; i=i+TOFEffStep){
        TOFEffAll.push_back( to_string_with_precision(i,2) );
    }
}
// For representative point
if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
    //TOFEffAll.push_back( to_string_with_precision(0.95,2) );
    TOFEffAll.push_back( to_string_with_precision(0.97,2) );
}


//// Loop Time bins
for ( int i = stoi(StartTimeIndex); i < maxindex; i = i+stoi(bartals) ){ //i from 0:first bartel rotation included, i from 1:first bartel rotation removed.  // Time index loop start
    std::cout << "Now Timebins is " << i << std::endl;

    //// Open ROOT File to Save result
    TFile *ResultROOT_SysErr_Shape;
    TFile *ResultROOT_SysErr_Eff;
    if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes"){    
        ResultROOT_SysErr_Shape = new TFile ((string("SysError_Shape/") + string("Ratio_") + timemode + string("_TimeIndex") + to_string(i) + string(".root")).c_str(),"RECREATE");
        ResultROOT_SysErr_Eff   = new TFile ((string("SysError_Eff/")   + string("Ratio_") + timemode + string("_TimeIndex") + to_string(i) + string(".root")).c_str(),"RECREATE");
    }


    //// Loop in TRD Efficiency
    for (auto &trdeff: TRDEffAll){  // Trd Efficiency loop start

        //// Loop in TOF Efficiency
        for (auto &tofeff: TOFEffAll){  // Tof Efficiency loop start 

            //// Loop in generate number template loop
            for (int numberindex = 0; numberindex < GenerateNumbers; numberindex=numberindex+1){ // generate number loop start

                //// Rigidity Loop
                for (int index=stoi(rigidity_start); index < stoi(rigidity_end); index=index+stoi(binmerge)){  // Rigidity loop start

                    if (timemode == "3BartalRotation" and i == 123){
                        //// fill vectors forr fit result
                        v_antiproton_number.push_back(0);
                        v_proton_number    .push_back(0);
                        v_ratio_result_tof .push_back(0);
                        v_chi2_tof         .push_back(0);
                        v_error            .push_back(0);

                        // Save Result in original mode (no parametrizalised templates)
                        vector<double> result_tof       = {0, 0, 0};
                        double         ProtonNumber     = 0;
                        vector<double> ResultError_tof  = {0, 0, 0};
                        double         CHI2dof_tof      = 0;

                        if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
                            SaveResult(bartals, binning, index, binmerge, result_tof, ProtonNumber, ResultError_tof, CHI2dof_tof);
                        }                        


                    }

                    else{
                        //// Load Histograms
                        TH2F *TofTRD_data_pass7_negative          = new TH2F();
                        TH2F *TofTRD_template_ElectronData        = new TH2F();
                        TH2F *TofTRD_template_PionData            = new TH2F();
                        TH2F *TofTRD_data_pass7_positive_template = new TH2F();        
                        TH2F *TofTRD_data_pass7_positive_data     = new TH2F();
                        tie(TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData, TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data) = LoadTemplates_TimeDependent(binmerge, binning, index, i, timemode, ParametrilizedMode, numberindex, SysErrEffStudyMode, trdeff, tofeff);


                        // Template Fit
                        MYUtilities::TemplateFitter2D templateFitter_tof(-1); //printlevel 0:standard, -1:quiet
                        // Signal
                        templateFitter_tof.AddTemplateHistogram(TofTRD_data_pass7_positive_template);
                        // Background
                        templateFitter_tof.AddTemplateHistogram(TofTRD_template_ElectronData);
                        templateFitter_tof.AddTemplateHistogram(TofTRD_template_PionData);
                        // Data
                        templateFitter_tof.SetDataHistogram(TofTRD_data_pass7_negative);
                        // Fit
                        templateFitter_tof.Fit(1,0,0); //2D Template Fit do not support fit ranges.

                        // Template Fit result  
                        vector<double> result_tof, ResultError_tof;
                        assert(result_tof.empty());
                        assert(ResultError_tof.empty());
                        result_tof.assign(templateFitter_tof.GetAbsoluteResult().begin(), templateFitter_tof.GetAbsoluteResult().end());
                        ResultError_tof.assign(templateFitter_tof.GetAbsoluteResultError().begin(), templateFitter_tof.GetAbsoluteResultError().end());
                        double Chi2_tof    = templateFitter_tof.Chi2();
                        int NDF_tof        = templateFitter_tof.NDF();
                        double CHI2dof_tof = Chi2_tof/NDF_tof;

                        // fill vectors forr fit result
                        v_antiproton_number.push_back(result_tof[0]);
                        double ProtonNumber = TofTRD_data_pass7_positive_data->GetEntries();
                        v_proton_number.push_back(ProtonNumber);
                        v_ratio_result_tof.push_back(result_tof[0]/ProtonNumber*100000); // proton from total number before cut
                        v_chi2_tof.push_back(CHI2dof_tof);
                        v_error.push_back(ResultError_tof[0]/ProtonNumber*100000);

                        // Plot for only original mode
                        if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
                            string issversion   = "";
                            string AnalysisMode = "TimeDependent";
                            // Plot Template and data
                            Plot_TemplatePion_Data      (binmerge, binning, index, issversion, TofTRD_template_PionData           , ParametrilizedMode, trdeff, tofeff, AnalysisMode, i);
                            Plot_TemplateElectron_Data  (binmerge, binning, index, issversion, TofTRD_template_ElectronData       , ParametrilizedMode, trdeff, tofeff, AnalysisMode, i);
                            Plot_ISSNegative            (binmerge, binning, index, issversion, TofTRD_data_pass7_negative         , ParametrilizedMode, trdeff, tofeff, AnalysisMode, i);
                            Plot_TemplateAntiproton_Data(binmerge, binning, index, issversion, TofTRD_data_pass7_positive_template, ParametrilizedMode, trdeff, tofeff, AnalysisMode, i);
                            if (templateFitter_tof.HasBeenFit() == true){                            
                                // Plot Fit Result and Projections
                                templateFitter_tof.CreateResultDrawing("Fit_Result_tof",1000,500)->SaveAs( (std::string("results/fitplot/FitResult_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::to_string(i) + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string(".pdf")).c_str()); 
                                templateFitter_tof.CreateResultDrawingXprojection("TrdLikelihood_X",1000,800)->SaveAs( (std::string("results/fitplot/TrdLikelihood_X_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::to_string(i) + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string(".pdf")).c_str());
                                templateFitter_tof.CreateResultDrawingYprojection("1/TofBeta_Y",1000,800)->SaveAs( (std::string("results/fitplot/1_TofBeta_Y_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::to_string(i) + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string(".pdf")).c_str());
                            }
                            
                        }
                        

                        // Save Result in original mode (no parametrizalised templates)
                        if (ParametrilizedMode == "No" and SysErrEffStudyMode == "No"){
                            SaveResult(bartals, binning, index, binmerge, result_tof, ProtonNumber, ResultError_tof, CHI2dof_tof);
                        }

                        // Clear
                        result_tof.clear();
                        ResultError_tof.clear();
          
                    }
 
                }  // rigidity loop end.


                // Make TGraph and TH1D to hold ratio result
                TGraph *g_ratio_tof   = new TGraph( v_ratio_result_tof.size(), subbincenter_binmerge.data(), v_ratio_result_tof.data());

                //// Save Fit Result in ROOT File
                if (ParametrilizedMode == "Yes" and SysErrEffStudyMode == "No"){
                    ResultROOT_SysErr_Shape->cd();
                    g_ratio_tof->Write( (std::string("g_ratio_tof") + to_string(numberindex)).c_str() );
                }

                if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes"){
                    ResultROOT_SysErr_Eff->cd();
                    g_ratio_tof->Write( (std::string("g_ratio_tof") + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() ); 
                }

                //// Clear
                v_ratio_result_tof.clear();
                v_chi2_tof.clear();
                v_error.clear();
                v_antiproton_number.clear();
                v_proton_number.clear();

            }// generate number loop end

        } // Tof Efficiency loop end

    } // Trd Efficiency loop end

    if (ParametrilizedMode == "No" and SysErrEffStudyMode == "Yes"){
        ResultROOT_SysErr_Shape->Close();
        ResultROOT_SysErr_Eff  ->Close();
    }

} // Time index loop end



return EXIT_SUCCESS;
}


