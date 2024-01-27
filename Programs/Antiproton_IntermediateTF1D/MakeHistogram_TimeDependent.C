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

#define INFO_OUT_TAG "Antiproton_IntermediateMakeTemplate_TimeDependent"
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

std::string binmerge = "";
config.GetValue("OPTIONS", "binmerge", binmerge,
              "how many bins used for each fit");

std::string signalefficiency = "";
config.GetValue("OPTIONS", "signalefficiency", signalefficiency,
              "The rich beta efficiency is");

std::string timemode = "";
config.GetValue("OPTIONS", "timemode", timemode,
              "Time mode you choose.");

std::string StartTimeIndex = "";
config.GetValue("OPTIONS", "StartTimeIndex", StartTimeIndex,
              "If you need to start from a certain TimeIndex");

std::string EndTimeIndex = "";
config.GetValue("OPTIONS", "EndTimeIndex",EndTimeIndex,
              "If you need to end from a certain TimeIndex");

std::string FullRange = "";
config.GetValue("OPTIONS", "FullRange", FullRange,
              "If you want to use FullRange of the TRDLikelihood");

string lowpath = (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/")).c_str();
chdir( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/")).c_str());
string AnylsisMode = "TimeAveraged";


string namesuffix;
if (timemode == "3BartalRotation"){
    namesuffix="_3BartalRotation";}
else if (timemode == "6BartalRotation"){
    namesuffix="_6BartalRotation";}
else if (timemode == "6months"){
    namesuffix="_6months";}


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
double TRDEffStep  = 0.01;

for (double i=TRDEffStart; i<=TRDEffEnd; i=i+TRDEffStep){
    TRDEffAll.push_back( to_string_with_precision(i,2) );
}


//// Load RichbetaCut
vector<double> v_richbetacut;
tie(v_richbetacut) = LoadRichbetaCut(lowpath, signalefficiency, binmerge);


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
for (int index = stoi(rigidity_start); index < stoi(rigidity_end); index = index + stoi(binmerge)){   // rigidity loop

    //// Loop Time bins
    for (int TimeIndex=stoi(StartTimeIndex); TimeIndex<maxindex; TimeIndex = TimeIndex + stoi(bartals)) { //TimeIndex from 0:first bartel rotation included, i from 1:first bartel rotation removed.    // Time mode index start
    
        std::cout << "Now Timebins is " << TimeIndex << std::endl;

        //// Create root file to save Template (Each Rigidity bin has a root file)
        std::string left  = doubleToString(binning[index]);
        std::string right = doubleToString(binning[index + stoi(binmerge)]);
        if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
            left = to_string_with_precision(binning[index], 1);}
        if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
            right = to_string_with_precision(binning[index + stoi(binmerge)], 1);}

        std::string RootFileName;
        if (FullRange == "Yes"){
            RootFileName = "TimeDependentTemplatesAndData_fullRange_";
            TRDEffAll.clear();
            TRDEffAll.push_back( to_string_with_precision(1.00,2) );
        }
        else if (FullRange == "No"){
            RootFileName = "TimeDependentTemplatesAndData_";;
        }

        TFile HistoRootFile( ((string("rootfiles") + namesuffix + string("_template/") + RootFileName + left + string("_") + right + string("_") + std::to_string(TimeIndex)) + string(".root")).c_str(), "RECREATE");
        std::cout << "left " << left << ", right: " << right << std::endl;

        //// Loop in TRD Efficiency
        for (auto &trdeff: TRDEffAll){  // Trd Efficiency loop start

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
            tie(data_pass7_positive) = LoadISSPositiveData_TimeDependent(lowpath, binmerge, binning, index, PositiveCut, timemode, TimeIndex, bartals, maxindex, TrdLOW, TrdHIGH, TrdBINNUMBER);  
            // ISS negative data
            TH1F *data_pass7_negative;
            TH1F *template_electron_Data;
            TH1F *template_pion_Data;    
            tie(data_pass7_negative, template_electron_Data, template_pion_Data) = LoadISSNegativeData_TimeDependent(lowpath, binmerge, binning, index, NegativeCut, ElectronTemplateDataCut, PionTemplateDataCut, timemode, TimeIndex, bartals, maxindex, TrdLOW, TrdHIGH, TrdBINNUMBER);

            //// Save Template in root file (Each Rigidity bin has a root file)
            HistoRootFile.cd();
            data_pass7_positive   ->Write( (std::string("data_pass7_positive")    + std::string("_TRDeff_") + trdeff).c_str() );
            template_electron_Data->Write( (std::string("template_electron_Data") + std::string("_TRDeff_") + trdeff).c_str() );
            template_pion_Data    ->Write( (std::string("template_pion_Data")     + std::string("_TRDeff_") + trdeff).c_str() );
            data_pass7_negative   ->Write( (std::string("data_pass7_negative")    + std::string("_TRDeff_") + trdeff).c_str() );


        }  // Trd Efficiency loop end

    HistoRootFile.Close();

    }  // // Time mode index end

}   // rigidity loop end



    return EXIT_SUCCESS;
}
