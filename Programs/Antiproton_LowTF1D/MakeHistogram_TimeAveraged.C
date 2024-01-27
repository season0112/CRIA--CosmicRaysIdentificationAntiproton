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

#define INFO_OUT_TAG "MakeHistogram_TimeAveraged"
#include "debugging.hh"


int main(int argc, char* argv[]) {

Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
config.ReadCommandLine(argc, argv);

config.SetProgramHelpText("Antiproton_LowMakeTemplate",
                            "Make Templates in Low Rigidity Range for template fit.");

config.AddHelpExample("Antiproton_LowMakeTemplate", "--issversion pass7.8 --rigidity_start 2 --rigidity_end 18 --binmerge 1  --Testmode Yes");

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

std::string FullRange = "";
config.GetValue("OPTIONS", "FullRange", FullRange,
              "If you want to use FullRange of the Fit Range");

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


//// Define Binnings
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


//// Load root files and apply cuts

//// Loop in Rigidity
for (int index = stoi(rigidity_start)+(stoi(binmerge)-1) ; index<stoi(rigidity_end)+(stoi(binmerge)-1); index=index+stoi(binmerge)) {   // rigidity loop start

    //// Create root file to save Template (Each Rigidity bin has a root file)
    std::string RootFileName;
    std::string RootFileOpenState;
    if (FullRange == "Yes"){
        RootFileName = "averaged_ratio_fullRange_";
        TRDEffAll.clear();
        TOFEffAll.clear();
        TRDEffAll.push_back( to_string_with_precision(1.00,2) );
        TOFEffAll.push_back( to_string_with_precision(1.00,2) );
        RootFileOpenState = "RECREATE";
    }
    else if (FullRange == "No"){
        RootFileName = "averaged_ratio_";;
        RootFileOpenState = "RECREATE";
    }

    TFile HistoRootFile((string("Time_Averaged_ratio_Low/binmerge") + binmerge + string("/tof/") + RootFileName + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+stoi(binmerge)]) + string("_") +issversion + string(".root")).c_str(), (RootFileOpenState).c_str() );

    cout<< string("Index now is:") + to_string(index) <<endl;

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


            //// standard fit range
            double RECITOFBETALOW;
            double RECITOFBETAHIGH;
            int RECITOFBETANUMBER;
            double TrdLOW;
            double TrdHIGH;
            int TrdBINNUMBER;

            if (FullRange == "Yes"){
                /*
                //Onebintoshow (XX-XXGV)
                RECITOFBETALOW  = -0.3;
                RECITOFBETAHIGH = 0.2;
                RECITOFBETANUMBER  = 20;
                TrdLOW    = (v_TRDcut.at(index-1));
                TrdHIGH   = 1.7;
                TrdBINNUMBER = 15;            
                */
                RECITOFBETALOW    = -0.3; // no name version:-0.3
                RECITOFBETAHIGH   = 0.2;
                RECITOFBETANUMBER = 40;
                TrdLOW            = 0.70; // no name version:0.88  
                TrdHIGH           = 1.6;
                TrdBINNUMBER      = 20;
            }
            else{
                //double RECITOFBETALOW = -0.7;
                RECITOFBETALOW    = v_TOFcut.at( (index - stoi(rigidity_start))/stoi(binmerge) );
                RECITOFBETAHIGH   = 0.3;
                RECITOFBETANUMBER = 40;
                TrdLOW            = v_TRDcut.at( (index - stoi(rigidity_start))/stoi(binmerge) );
                TrdHIGH           = 1.7;
                TrdBINNUMBER      = 20;
            }



            //// Allcut definations.     TRD_electron_proton: all, TRD_He_Proton: Antiproton Template, positivedata,
            // Tracker
            std::string trackerpattern = "0";
            TCut PatternCut                            = (string("Pattern==") + trackerpattern).c_str();
            TCut ExtrapolatedPhotoElectronsFromTracker ="NPhotoElectrons-ExtrapolatedRichExpectedPhotoElectronsProton<20 && NPhotoElectrons>-999 && ExtrapolatedRichTileIndex == TileIndex || NPhotoElectrons==-999 || ExtrapolatedRichTileIndex != TileIndex";
            TCut ProtonCCMVABDCut                      = "ProtonCCMVABDT > 0.9";
            // Trigger
            TCut TriggerPhysics                        = "TriggerFlags != 1 && TriggerFlags != 64 && TriggerFlags != 65";
            // TRD
            //TCut TrdLikelihoodCut                    = (string("TrdLogLikelihoodRatioElectronProtonTracker>") + doubleToString(TrdLOW) + string("|| TrdLogLikelihoodRatioElectronProtonTracker==-1.5") ).c_str(); // Antiproton:1to1.2, Electron:0.5
            TCut TrdLikelihoodCut                      = (string("TrdLogLikelihoodRatioElectronProtonTracker>") + doubleToString(TrdLOW)).c_str(); // Antiproton:1to1.2, Electron:0.5
            TCut TrdLikelihoodHeProtonCut              = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.1";
            // ECAL
            TCut EcalBDT_EnergyDCut                    = "EcalBDT_EnergyD < -0.9";
            // RICH  
            //TCut RichBetaCut                         = (string("RichBeta<") + to_string(v_richbetacut.at((index-stoi(rigidity_start)-(stoi(binmerge)-1))/stoi(binmerge) )) ).c_str();
            TCut RichBetaCut                           = "RichBeta==0 || RichIsNaF==1";
            //TCut RichBetaCut                         = "RichIsNaF==1";
            TCut RichNaFCut                            = "RichIsNaF==1";
            TCut RichAglCut                            = "RichIsNaF==0";
            TCut RichPhotoElectron                     = "NPhotoElectrons-NExpectedPhotoElectrons<10  || NPhotoElectrons==-999";
            TCut RichChargeCut                         = "RichCharge<2 || RichCharge==0";
            // New cuts added:
            TCut TrdSegmentsXZNumberCut                = "TrdSegmentsXZNumber==1";
            TCut TrdSegmentsYZNumberCut                = "TrdSegmentsYZNumber==1";
            TCut TRDVTracksSizeCut                     = "TRDVTracksSize==1";
            TCut TrdNumberOfHitsCut                    = "TrdNumberOfHits<35";
            TCut ACCHitsCut                            = "ACCHits==0";

            vector<double> v_TRDHit = {55, 48, 45, 43, 44, 43, 42, 42, 42, 40, 41, 41, 39, 40, 40, 40}; 
            double TRDHit = v_TRDHit.at(index-1); 
            //TCut TrdNumberOfHitsCut_P                = (string("TrdNumberOfHits<") + doubleToString(TRDHit)).c_str();  // To match for PhyReport2021 (with NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut)
            TCut TrdNumberOfHitsCut_P                  = "TrdNumberOfHits<40";

            // Templates and data cuts
            TCut AntiprotonCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
            //TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
            //TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && "TrdSegmentsXZNumber==1 || TrdSegmentsXZNumber==2" && "TrdSegmentsYZNumber==1 || TrdSegmentsYZNumber==2" && "TRDVTracksSize==1" && "TrdNumberOfHits<40" && "ACCHits==0 || ACCHits==1"; //set5
            
            // new
            TCut NegativeCut      = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut;    
            TCut PositiveCut      = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TRDVTracksSizeCut && TrdNumberOfHitsCut_P;
            TCut PositiveCut_test = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TRDVTracksSizeCut;
            // old
            //TCut NegativeCut      = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
            //TCut PositiveCut      = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut; // (originali)
            //TCut PositiveCut_test = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;

            //// Pion and electron template from Data  (Mass_electron=0.000511)
            //// Set1:
            //TCut PionTemplateDataCut     = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>0.72";
            //TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-1)<0.002 && TrdLogLikelihoodRatioElectronProtonTracker<0.55";
            //TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-1)<0.002";
            //TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002 && TrdSegmentsXZNumber==1 && TrdSegmentsYZNumber==1 && TrdNumberOfHits<35" ;
            //TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker<0.55";
            //TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002 && EcalBDT_EnergyD>0.9";
            //// Set2:
            //TCut ElectronTemplateDataCut = "EcalBDT_EnergyD>-999 && EcalBDT_EnergyD>0.5";
            //TCut PionTemplateDataCut     = "EcalBDT_EnergyD>-999 && EcalBDT_EnergyD<-0.5";
            //TCut PionTemplateDataCut;
            //// Set: ele14,pion145
            //TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002 && TrdSegmentsXZNumber==1 && TrdSegmentsYZNumber==1 && TrdNumberOfHits<35" ;
            TCut ElectronTemplateDataCut = "RichIsNaF==0 && RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2)>-0.002 && TrdSegmentsXZNumber==1 && TrdSegmentsYZNumber==1 && TrdNumberOfHits<35" ;
            TCut PionTemplateDataCut = (std::string("RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>") + std::to_string(TrdLOW) + std::string("&& TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>50")).c_str();
            //// Set: ele24,pion245
            //TCut ElectronTemplateDataCut = "EcalBDT_EnergyD>0.5 && TrdSegmentsXZNumber==1 && TrdSegmentsYZNumber==1 && TrdNumberOfHits<35"; 
            //TCut PionTemplateDataCut = (std::string("EcalBDT_EnergyD>-1 && EcalBDT_EnergyD<0.5 && TrdLogLikelihoodRatioElectronProtonTracker>") + std::to_string(TrdLOW) + std::string("&& TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>50")).c_str();     
            
            TCut TrdCut = ( string("TrdLogLikelihoodRatioElectronProtonTracker>") + doubleToString(TrdLOW) + string("&&") + string("TrdLogLikelihoodRatioElectronProtonTracker<") + doubleToString(TrdHIGH) ).c_str();
            TCut TOFBETACut =  (string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ) ") + string(">") + to_string(RECITOFBETALOW) + string("&&") + string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ) ") + string("<") + to_string(RECITOFBETAHIGH) ).c_str();


            //// Load data, Apply cut and Make histograms
            // ISS positive data
            TH2F *TofTRD_data_pass7_positive_template;
            TH2F *TofTRD_data_pass7_positive_data;
            TH2F *TofTRD_data_pass7_positive_data_test;
            tie(TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data, TofTRD_data_pass7_positive_data_test) = LoadISSPositiveData(binmerge, binning, index, TOFBETACut, TrdCut, PositiveCut, PositiveCut_test, issversion, TestMode, TrdBINNUMBER, TrdLOW, TrdHIGH, RECITOFBETANUMBER, RECITOFBETALOW, RECITOFBETAHIGH);
            // ISS negaitve data (NegativeCut, ElectronTemplateDataCut, PionTemplateDataCut)
            TH2F *TofTRD_data_pass7_negative;
            TH2F *TofTRD_template_ElectronData;
            TH2F *TofTRD_template_PionData;
            tie(TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData) = LoadISSNegativeData(binmerge, binning, index, TOFBETACut, TrdCut, NegativeCut, ElectronTemplateDataCut, PionTemplateDataCut, issversion, TestMode, TrdBINNUMBER, TrdLOW, TrdHIGH, RECITOFBETANUMBER, RECITOFBETALOW, RECITOFBETAHIGH);

            //// Save Template in Root File
            HistoRootFile.cd();

            if (FullRange == "No"){            
                TofTRD_data_pass7_positive_template  ->Write( (std::string("TofTRD_data_pass7_positive_template")  + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff ).c_str() );
                TofTRD_template_ElectronData         ->Write( (std::string("TofTRD_template_ElectronData")         + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff ).c_str() );
                TofTRD_template_PionData             ->Write( (std::string("TofTRD_template_PionData")             + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff ).c_str() );
                TofTRD_data_pass7_negative           ->Write( (std::string("TofTRD_data_pass7_negative")           + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff ).c_str() );
                TofTRD_data_pass7_positive_data      ->Write( (std::string("TofTRD_data_pass7_positive_data")      + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff ).c_str() );
                TofTRD_data_pass7_positive_data_test ->Write( (std::string("TofTRD_data_pass7_positive_data_test") + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff ).c_str() );
            }
            else if (FullRange == "Yes"){
                TofTRD_data_pass7_positive_template  ->Write( (std::string("TofTRD_data_pass7_positive_template")  + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string("_RECITOFBETALOW_") + to_string_with_precision(RECITOFBETALOW,2) + std::string("_RECITOFBETAHIGH_") + to_string_with_precision(RECITOFBETAHIGH,2) + std::string("_RECITOFBETANUMBER_") + to_string_with_precision(RECITOFBETANUMBER,2) + std::string("_TrdLOW_") + to_string_with_precision(TrdLOW,2) + std::string("_TrdHIGH_") + to_string_with_precision(TrdHIGH,2) + std::string("_TrdBINNUMBER_") + to_string_with_precision(TrdBINNUMBER,2)      ).c_str() );
                TofTRD_template_ElectronData         ->Write( (std::string("TofTRD_template_ElectronData")         + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string("_RECITOFBETALOW_") + to_string_with_precision(RECITOFBETALOW,2) + std::string("_RECITOFBETAHIGH_") + to_string_with_precision(RECITOFBETAHIGH,2) + std::string("_RECITOFBETANUMBER_") + to_string_with_precision(RECITOFBETANUMBER,2) + std::string("_TrdLOW_") + to_string_with_precision(TrdLOW,2) + std::string("_TrdHIGH_") + to_string_with_precision(TrdHIGH,2) + std::string("_TrdBINNUMBER_") + to_string_with_precision(TrdBINNUMBER,2) ).c_str() );
                TofTRD_template_PionData             ->Write( (std::string("TofTRD_template_PionData")             + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string("_RECITOFBETALOW_") + to_string_with_precision(RECITOFBETALOW,2) + std::string("_RECITOFBETAHIGH_") + to_string_with_precision(RECITOFBETAHIGH,2) + std::string("_RECITOFBETANUMBER_") + to_string_with_precision(RECITOFBETANUMBER,2) + std::string("_TrdLOW_") + to_string_with_precision(TrdLOW,2) + std::string("_TrdHIGH_") + to_string_with_precision(TrdHIGH,2) + std::string("_TrdBINNUMBER_") + to_string_with_precision(TrdBINNUMBER,2)).c_str() );
                TofTRD_data_pass7_negative           ->Write( (std::string("TofTRD_data_pass7_negative")           + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string("_RECITOFBETALOW_") + to_string_with_precision(RECITOFBETALOW,2) + std::string("_RECITOFBETAHIGH_") + to_string_with_precision(RECITOFBETAHIGH,2) + std::string("_RECITOFBETANUMBER_") + to_string_with_precision(RECITOFBETANUMBER,2) + std::string("_TrdLOW_") + to_string_with_precision(TrdLOW,2) + std::string("_TrdHIGH_") + to_string_with_precision(TrdHIGH,2) + std::string("_TrdBINNUMBER_") + to_string_with_precision(TrdBINNUMBER,2)).c_str() );
                TofTRD_data_pass7_positive_data      ->Write( (std::string("TofTRD_data_pass7_positive_data")      + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string("_RECITOFBETALOW_") + to_string_with_precision(RECITOFBETALOW,2) + std::string("_RECITOFBETAHIGH_") + to_string_with_precision(RECITOFBETAHIGH,2) + std::string("_RECITOFBETANUMBER_") + to_string_with_precision(RECITOFBETANUMBER,2) + std::string("_TrdLOW_") + to_string_with_precision(TrdLOW,2) + std::string("_TrdHIGH_") + to_string_with_precision(TrdHIGH,2) + std::string("_TrdBINNUMBER_") + to_string_with_precision(TrdBINNUMBER,2)).c_str() );
                TofTRD_data_pass7_positive_data_test ->Write( (std::string("TofTRD_data_pass7_positive_data_test") + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string("_RECITOFBETALOW_") + to_string_with_precision(RECITOFBETALOW,2) + std::string("_RECITOFBETAHIGH_") + to_string_with_precision(RECITOFBETAHIGH,2) + std::string("_RECITOFBETANUMBER_") + to_string_with_precision(RECITOFBETANUMBER,2) + std::string("_TrdLOW_") + to_string_with_precision(TrdLOW,2) + std::string("_TrdHIGH_") + to_string_with_precision(TrdHIGH,2) + std::string("_TrdBINNUMBER_") + to_string_with_precision(TrdBINNUMBER,2)).c_str() );
            }

        } // Tof Efficiency loop end

    }  // Trd Efficiency loop end

    HistoRootFile.Close();
        
}  // rigidity loop end.


  return EXIT_SUCCESS;
}
