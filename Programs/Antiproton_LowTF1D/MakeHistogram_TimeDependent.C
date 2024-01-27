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

#define INFO_OUT_TAG "MakeHistogram_TimeDependent"
#include "debugging.hh"


int main(int argc, char* argv[]) {

Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
config.ReadCommandLine(argc, argv);

config.SetProgramHelpText("MakeHistogram_TimeDependent",
                            "Make time dependent Templates in Low Rigidity Range for template fit");

config.AddHelpExample("Antiproton_LowMakeTimeDependentTemplate", "--timemode 6months --rigidity_start 4 --rigidity_end 6 --binmerge 2");

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


string namesuffix;
if (timemode == "3BartalRotation"){
    namesuffix="_3BartalRotation";}
else if (timemode == "6BartalRotation"){
    namesuffix="_6BartalRotation";}
else if (timemode == "6months"){
    namesuffix="_6months";}


//// Define binnings
vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());


//// Set TRD Efficiency for Systematic Error Uncertainty Calculation
vector<string> TRDEffAll;
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


//// Set TOF Beta Efficiency for Systematic Error Uncertainty Calculation
vector<string> TOFEffAll;
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


//// Load TimeStape for Setings
string bartals = "";
string timeindex = "";
int maxindex;
tie(bartals, timeindex, maxindex) = LoadTimeStapeSetings(timemode);
//// Rest Time loop end if EndTimeIndex has value
if (EndTimeIndex != ""){
    maxindex = stoi(EndTimeIndex);
}


//// Loop in Rigidity
for (int RigidityIndex=stoi(rigidity_start); RigidityIndex < stoi(rigidity_end); RigidityIndex = RigidityIndex+stoi(binmerge)){  // Rigidity Loop Start

    //// Loop Time bins 
    for (int TimeIndex=stoi(StartTimeIndex); TimeIndex<maxindex; TimeIndex = TimeIndex + stoi(bartals)) { //TimeIndex from 0:first bartel rotation included, i from 1:first bartel rotation removed.    // Time mode index start

        std::cout << "Now Timebins is " << TimeIndex << std::endl;

        //// Create root file to save Template (Each Rigidity bin has a root file)
        std::string left  = doubleToString(binning[RigidityIndex]);
        std::string right = doubleToString(binning[RigidityIndex + stoi(binmerge)]);
        if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
            left = to_string_with_precision(binning[RigidityIndex], 1);}
        if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
            right = to_string_with_precision(binning[RigidityIndex + stoi(binmerge)], 1);}
        TFile HistoRootFile( ((string("rootfiles") + namesuffix + string("_template/TimeDependentTemplatesAndData_") + left + string("_") + right + string("_") + std::to_string(TimeIndex)) + string(".root")).c_str(), "RECREATE");


        //// Loop in TRD Efficiency
        for (auto &trdeff: TRDEffAll){  // Trd Efficiency loop start

            //// Load Trdlikelihood Cut
            vector<double> v_TRDcut;
            tie(v_TRDcut) = LoadTRDLogLikelihoodCut(trdeff, binmerge);

            //// Loop in TOFBeta
            for (auto &tofeff: TOFEffAll){  // Tof Efficiency loop start                

                //// Load TOFBeta Cut
                vector<double> v_TOFcut;
                tie(v_TOFcut) = LoadTOFBetaCut(tofeff, binmerge);            


                // Fill histograms
                // 2D templates
                //double RECITOFBETALOW = -0.7;
                double RECITOFBETALOW   = v_TOFcut.at( (RigidityIndex - stoi(rigidity_start))/stoi(binmerge) );
                double RECITOFBETAHIGH  = 0.3;
                int RECITOFBETANUMBER   = 40;
                double TrdLOW           = v_TRDcut.at( (RigidityIndex - stoi(rigidity_start))/stoi(binmerge) );
                double TrdHIGH          = 1.7;
                int TrdBINNUMBER        = 20;

                // Cuts definations (Most of them are in TimeDependentHeader.hh)
                std::string trackerpattern = "0";
                TCut PatternCut = (string("Pattern==") + trackerpattern).c_str();
                TCut TrdLikelihoodCut = (std::string("TrdLogLikelihoodRatioElectronProtonTracker>") + doubleToString(TrdLOW) + std::string("|| TrdLogLikelihoodRatioElectronProtonTracker==-1.5") ).c_str();

                // Templates and data cuts
                //TCut ElectronCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut ;
                //TCut PionCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && RichBetaCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
                TCut AntiprotonCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;

                // v10.0
                //TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
                //TCut PositiveCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
                // v11.0
                TCut TrdNumberOfHitsCut_P = "TrdNumberOfHits<40";
                TCut NegativeCut   = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut;
                TCut PositiveCut   = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TRDVTracksSizeCut && TrdNumberOfHitsCut_P;

                // Pion and electron template from Data  (Mass_electron=0.000511)
                // Set: ele14,pion145
                //TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002 && TrdSegmentsXZNumber==1 && TrdSegmentsYZNumber==1 && TrdNumberOfHits<35" ;
                //TCut PionTemplateDataCut = (std::string("RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>") + std::to_string(TrdLOW) + std::string("&& TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>50")).c_str();
                //// Set: ele24,pion245  (official for v10.0)
                //TCut ElectronTemplateDataCut = "EcalBDT_EnergyD>0.5 && TrdSegmentsXZNumber==1 && TrdSegmentsYZNumber==1 && TrdNumberOfHits<35";
                //TCut PionTemplateDataCut = (std::string("EcalBDT_EnergyD>-1 && EcalBDT_EnergyD<0.5 && TrdLogLikelihoodRatioElectronProtonTracker>") + std::to_string(TrdLOW) + std::string("&& TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>50")).c_str();
                //// Set: v11.0, same as time averaged
                TCut ElectronTemplateDataCut = "RichIsNaF==0 && RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2)>-0.002 && TrdSegmentsXZNumber==1 && TrdSegmentsYZNumber==1 && TrdNumberOfHits<35";
                TCut PionTemplateDataCut = (std::string("RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>") + std::to_string(TrdLOW) + std::string("&& TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>50")).c_str();        


                // Load data, apply cut and draw histograms
                // ISS negative data
                TH2F *TofTRD_data_pass7_negative;
                TH2F *TofTRD_template_ElectronData;
                TH2F *TofTRD_template_PionData;
                tie(TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData) = LoadISSNegativeData_TimeDependent(binmerge, binning, RigidityIndex, NegativeCut, ElectronTemplateDataCut, PionTemplateDataCut, TrdBINNUMBER, TrdLOW, TrdHIGH, RECITOFBETANUMBER, RECITOFBETALOW, RECITOFBETAHIGH, TimeIndex, bartals, maxindex, timemode);
                // ISS positive data
                TH2F *TofTRD_data_pass7_positive_template;
                TH2F *TofTRD_data_pass7_positive_data;
                tie(TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data) = LoadISSPositiveData_TimeDependent(binmerge, binning, RigidityIndex, AntiprotonCut, PositiveCut, TrdBINNUMBER, TrdLOW, TrdHIGH, RECITOFBETANUMBER, RECITOFBETALOW, RECITOFBETAHIGH, TimeIndex, bartals, maxindex, timemode);


                //// Save Time Dependent Template
                HistoRootFile.cd();
                TofTRD_data_pass7_negative          ->Write( (std::string("TofTRD_data_pass7_negative")          + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                TofTRD_template_ElectronData        ->Write( (std::string("TofTRD_template_ElectronData")        + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                TofTRD_template_PionData            ->Write( (std::string("TofTRD_template_PionData")            + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                TofTRD_data_pass7_positive_template ->Write( (std::string("TofTRD_data_pass7_positive_template") + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );
                TofTRD_data_pass7_positive_data     ->Write( (std::string("TofTRD_data_pass7_positive_data")     + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff).c_str() );

            }  // Tof Efficiency loop start

        }  // Trd Efficiency loop start

        HistoRootFile.Close();

    } // Time mode index end

}  // rigidity loop end.

  return EXIT_SUCCESS;
}


