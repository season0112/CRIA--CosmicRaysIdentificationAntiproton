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

#define INFO_OUT_TAG "Antiproton_IntermediateMakeTemplate"
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

std::string binmerge = "";
config.GetValue("OPTIONS", "binmerge", binmerge,
              "how many bins used for each fit");

std::string signalefficiency = "";
config.GetValue("OPTIONS", "signalefficiency", signalefficiency,
              "The rich beta efficiency is");

std::string TestMode = "";
config.GetValue("OPTIONS", "Testmode", TestMode,
              "If you want to use TestMode to run program quickly");

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
string AnylsisMode = "TimeAveraged";


//// Define Binnings
vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());
vector<double> bincenter (AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins());
//vector<double> subbincenter(bincenter.begin()+stoi(rigidity_start), bincenter.begin()+stoi(rigidity_end));
//vector<double> v_pattern0_percentage, v_pattern1_percentage, v_pattern2_percentage, v_pattern4_percentage;


//// Set TRD Efficiency for Systematic Error Uncertainty Calculation
vector<string> TRDEffAll;
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


//// Load RichbetaCut
vector<double> v_richbetacut;
tie(v_richbetacut) = LoadRichbetaCut(lowpath, signalefficiency, binmerge);



//// Load root files and apply cuts

//// Loop in Rigidity
for (int index = stoi(rigidity_start); index < stoi(rigidity_end); index = index + stoi(binmerge)){   // rigidity loop


    //// Create root file to save Template (Each Rigidity bin has a root file)
    std::string left  = doubleToString(binning[index]);
    std::string right = doubleToString(binning[index+stoi(binmerge)]); 
    if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
        left = to_string_with_precision(binning[index],1);}
    if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
        right = to_string_with_precision(binning[index+stoi(binmerge)], 1);}

    std::string RootFileName;
    if (FullRange == "Yes"){
        RootFileName = "averaged_ratio_fullRange_";
        TRDEffAll.clear();
        TRDEffAll.push_back( to_string_with_precision(1.00,2) );
    }
    else if (FullRange == "No"){
        RootFileName = "averaged_ratio_";;
    } 

    TFile HistoRootFile((string("Time_Averaged_ratio/binmerge") + binmerge + string("/RigidityRootFiles/") + RootFileName + left + string("_") + right + string("_") + issversion + string(".root")).c_str(), "RECREATE");
    cout<< "\n" <<endl;
    cout<< "Index now is: "  << index          << endl;
    cout<< "Binning now is " << left << "_" << right << endl;


    //// Loop in TRD Efficiency
    for (auto &trdeff: TRDEffAll){  // Trd Efficiency loop start
        cout<< "TRD Efficiency is " << trdeff <<endl;

        //// Load Trdlikelihood Cut
        vector<double> v_TRDcut;
        tie(v_TRDcut) = LoadTRDLogLikelihoodCut(trdeff, binmerge, lowpath);

        //// standard fit range
        double TrdLOW    = -2.0;
        double TrdHIGH;
        if (FullRange == "Yes"){
            TrdHIGH   = 0.0;}
        else{
            TrdHIGH   = ( (-1.0) * v_TRDcut.at( (index - stoi(rigidity_start))/stoi(binmerge) ));}

        int TrdBINNUMBER = 50;

        //// Cut Defination
        TCut RichBetaCut = (string("RichBeta<") + to_string( v_richbetacut.at( (index - stoi(rigidity_start))/stoi(binmerge))) ).c_str();
        TCut PatternCut;
        tie(PatternCut) = patterncut(trackerpattern);
        TCut pattern0_percentage =  RichBetaCut && Pattern0;
        TCut pattern1_percentage =  RichBetaCut && Pattern1;
        TCut pattern2_percentage =  RichBetaCut && Pattern2;
        TCut pattern4_percentage =  RichBetaCut && Pattern4;
        // Templates selectron
        //TCut NegativeCut = RichBetaCut && PatternCut && TrdLikelihoodHeProtonCut_proton && EcalBDT_EnergyDCut_proton && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
        //TCut PositiveCut = RichBetaCut && PatternCut && TrdLikelihoodHeProtonCut_proton && EcalBDT_EnergyDCut_proton && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
        // old:
        //TCut NegativeCut = RichBetaCut && PatternCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
        //TCut PositiveCut = RichBetaCut && PatternCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut && TRDVTracksSizeCut && TrdNumberOfHitsCut;
        // new:
        TCut NegativeCut = RichBetaCut && PatternCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut;
        TCut PositiveCut = RichBetaCut && PatternCut && TrdSegmentsXZNumberCut && TrdSegmentsYZNumberCut;
        // Pion and electron template from Data  (Mass_electron=0.000511) (Other tries see CutDefinition.hh) 
        TCut ElectronTemplateDataCut = PatternCut && "TrdLogLikelihoodRatioElectronProtonTracker>-1.5" && "EcalBDT_EnergyD>0.5" && "TrdSegmentsXZNumber==1" && "TrdSegmentsYZNumber==1" && "TrdNumberOfHits<35";
        TCut PionTemplateDataCut = PatternCut && "TrdLogLikelihoodRatioElectronProtonTracker>-1.5" && "EcalBDT_EnergyD>-1" && "EcalBDT_EnergyD<0.5" && "TrdSegmentsXZNumber>1" && "TrdSegmentsYZNumber>1" && "TrdNumberOfHits>40";


        //// Load Root File and make Histograms
        // ISS positive data
        TH1F *data_pass7_positive;
        TChain *fpass7_positive;
        tie(data_pass7_positive, fpass7_positive) = LoadISSPositiveData(lowpath, binmerge, binning, index, PositiveCut, issversion, TestMode, TrdLOW, TrdHIGH, TrdBINNUMBER);
        // ISS negative data
        TH1F *data_pass7_negative;
        TH1F *template_electron_Data;
        TH1F *template_pion_Data;    
        TChain *fpass7_negative;
        tie(data_pass7_negative, template_electron_Data, template_pion_Data, fpass7_negative) = LoadISSNegativeData(lowpath, binmerge, binning, index, NegativeCut, ElectronTemplateDataCut, PionTemplateDataCut, issversion, TestMode, TrdLOW, TrdHIGH, TrdBINNUMBER);


        //// Plot Tracker Pattern Percentage 
        //Plot_Pattern_Percentage(richcut, issversion, binmerge, v_pattern0_percentage, v_pattern1_percentage, v_pattern2_percentage, v_pattern4_percentage, subbincenter, fpass7_negative, RichBetaCut, pattern0_percentage, pattern1_percentage, pattern2_percentage, pattern4_percentage, AllCut);


        //// Save Template in root file (Each Rigidity bin has a root file)
        HistoRootFile.cd();
        data_pass7_positive   ->Write( (std::string("data_pass7_positive") + std::string("_TRDeff_")    + trdeff).c_str() );
        template_electron_Data->Write( (std::string("template_electron_Data") + std::string("_TRDeff_") + trdeff).c_str() );
        template_pion_Data    ->Write( (std::string("template_pion_Data") + std::string("_TRDeff_")     + trdeff).c_str() );
        data_pass7_negative   ->Write( (std::string("data_pass7_negative") + std::string("_TRDeff_")    + trdeff).c_str() );


        //// Clear 
        fpass7_negative->Reset();
        delete fpass7_negative;
        fpass7_positive->Reset();
        delete fpass7_positive;

    }  // Trd Efficiency loop end

    HistoRootFile.Close();

}   // rigidity loop end



    return EXIT_SUCCESS;
}
