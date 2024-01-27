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
#include <vector>
#include <TCanvas.h>
#include <TF1.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include "Utilities.hh"

#include "AntiprotonBinning.hh"
#include "TemplateFitterforLowEnergy.hh"
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

using namespace std;

#define INFO_OUT_TAG "Antiproton_LowTF1D"
#include "debugging.hh"

int main(int argc, char* argv[]) {

Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
config.ReadCommandLine(argc, argv);

config.SetProgramHelpText("Antiproton_LowTF1D",
                            "Template fit for antiproton signal determination in Low Rigidity Range.");

config.AddHelpExample("Antiproton_LowTF1D", "--issversion pass7.8 --rigidity_start 4 --rigidity_end 6");

std::string issversion = "";
config.GetValue("OPTIONS", "issversion", issversion,
              "The Issversion is");

std::string rigidity_start = "";
config.GetValue("OPTIONS", "rigidity_start", rigidity_start,
              "The rigidity_start_index is");

std::string rigidity_end = "";
config.GetValue("OPTIONS", "rigidity_end", rigidity_end,
              "The rigidity_end_index is");

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

////////////////////////////////////////////////////
vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());
vector<double> bincenter (AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins());



vector<double> v_ratio_result, v_chi2, v_ratio_result_with_effective;
vector<double> v_ratio_result_tof, v_chi2_tof, v_ratio_result_tof_with_effective;
assert(v_ratio_result.empty());
assert(v_chi2.empty());
assert(v_ratio_result_tof.empty());
assert(v_chi2_tof.empty());
assert(v_ratio_result_with_effective.empty());
assert(v_ratio_result_tof_with_effective.empty());


TFile *f_eff_antiproton = new TFile("EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root");
TFile *f_eff_proton = new TFile("EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root");
TGraphAsymmErrors *Acceptance_antiproton = (TGraphAsymmErrors*)f_eff_antiproton->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
TGraphAsymmErrors *Acceptance_proton = (TGraphAsymmErrors*)f_eff_proton->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
int number;
number = Acceptance_antiproton->GetN();
Double_t xa[number],Aa[number],ea[number],xp[number],Ap[number],ep[number];
TGraphErrors *Effective_Acceptance = new TGraphErrors(number);
for (int p = 0; p < number; p++) {       //from 0, 0 is 0.8-1.0, last is from 822 to 1130. so for B1042MC there two are 0 (generated momentum from 1 to 800).
Acceptance_antiproton->GetPoint(p,xa[p],Aa[p]);
Acceptance_proton->GetPoint(p,xp[p],Ap[p]);
ea[p] = Acceptance_antiproton->GetErrorY(p);
ep[p] = Acceptance_proton->GetErrorY(p);
Effective_Acceptance->SetPoint(p,xa[p],Ap[p]/Aa[p]);
Effective_Acceptance->SetPointError(p, 0.0, sqrt(pow(ep[p],2)/pow(Aa[p],2)+pow(Ap[p],2)/pow(Aa[p],4)*pow(ea[p],2)));
}




for (int index=stoi(rigidity_start); index<stoi(rigidity_end); index++)
{
    string negativename;
    string positivename;
    string antiprotonname = string("B1042_antipr.pl1.1800_7.6_all_Tree_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    string electronname = string("B1091_el.pl1.0_25200_7.6_all_Tree_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    string pionname = string("B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    string deuteronname = string("B1128_d.pl1ph.021000_7.6_all_Tree_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");

    if (issversion == "pass7.8"){
        negativename = string("B1130_pass7_7.7_all_Tree_negative_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
        positivename = string("B1130_pass7_7.7_all_Tree_positive_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    }
    else if (issversion == "2016paper"){
        negativename = string("B1130_pass7_7.7_all_Tree_negative_May2015_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
        positivename = string("B1130_pass7_7.7_all_Tree_positive_May2015_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    }
    TFile *fantiproton = new TFile(antiprotonname.c_str());
    TFile *felectron = new TFile(electronname.c_str());
    TFile *fpass7_negative = new TFile(negativename.c_str());
    TFile *fpass7_positive = new TFile(positivename.c_str());
    TFile *fpion = new TFile(pionname.c_str());
    TFile *fdeuteron = new TFile(deuteronname.c_str());


    TTree *tantiproton = (TTree*)fantiproton->Get("AntiprotonIntermediateEnergyTree");
    TTree *telectron = (TTree*)felectron->Get("AntiprotonIntermediateEnergyTree");
    TTree *tpass7_negative = (TTree*)fpass7_negative->Get("AntiprotonIntermediateEnergyTree");
    TTree *tpass7_positive = (TTree*)fpass7_positive->Get("AntiprotonIntermediateEnergyTree");
    TTree *tpion = (TTree*)fpion->Get("AntiprotonIntermediateEnergyTree");
    TTree *tdeuteron = (TTree*)fdeuteron->Get("AntiprotonIntermediateEnergyTree");

/*
    // Mass ^2 template 
    TH1F *template_antiproton = new TH1F("Mass^2_p","",50,-1.3,2.5);
    TH1F *template_antiproton_parametrized = new TH1F("Mass^2_p_parametrilized","",50,-1.3,2.5);
    TH1F *template_electron = new TH1F("Mass^2_e","",50,-1.3,2.5);
    TH1F *data_pass7_negative = new TH1F("Mass^2_pass7_negative","",50,-1.3,2.5);
    TH1F *data_pass7_negative_fit_parametrized = new TH1F("Mass^2_pass7_negative_parametrilized","",50,-1.3,2.5);
    TH1F *data_pass7_positive = new TH1F("Mass^2_pass7_positive","",50,-1.3,2.5);
    TH1F *data_pass7_positive_fit_parametrized = new TH1F("Mass^2_pass7_positive_parametrilized","",50,-1.3,2.5);
    TH1F *template_pion = new TH1F("Mass^2_pion","",50,-1.3,2.5); //pion mass:0.139GeV, mass^2 = 0.01932.
    //TH1F *template_deuteron = new TH1F("Mass^2_deuteron","",50,-1.3,2.5);  //deuteron mass: 1.875GeV, mass^2 = 3.515.
    //TH1F *template_deuteron_parametrized = new TH1F("Mass^2_deuteron_parametrilized","",50,-1.3,2.5);
    TH1F *template_deuteron = new TH1F("Mass^2_deuteron","",50,-1.3,2.5);  //deuteron mass: 1.875GeV, mass^2 = 3.515.
    TH1F *template_deuteron_parametrized = new TH1F("Mass^2_deuteron_parametrilized","",50,-1.3,2.5);
*/

    // Tof_Beta template
    TH1F *Tof_template_electron = new TH1F("Tof_Beta_e","",30,0.6,1.5);
    TH1F *Tof_template_antiproton = new TH1F("Tof_Beta_p","",30,0.6,1.5);
    TH1F *Tof_template_antiproton_parametrized = new TH1F("Tof_Beta_p_parametrilized","",30,0.6,1.5);
    TH1F *Tof_data_pass7_negative = new TH1F("Tof_Beta_pass7_negative","",30,0.6,1.5);
    TH1F *Tof_data_pass7_negative_fit_parametrized = new TH1F("Tof_data_pass7_negative_fit_parametrized","",30,0.6,1.5);
    TH1F *Tof_data_pass7_positive = new TH1F("Tof_data_pass7_positive","",30,0.6,1.5);
    TH1F *Tof_data_pass7_positive_fit_parametrized = new TH1F("Tof_data_pass7_positive_fit_parametrized","",30,0.6,1.5);
    TH1F *Tof_template_pion = new TH1F("Tof_Beta_pion","",30,0.6,1.5);



    Float_t TrdLogLikelihoodRatioElectronProtonTracker_p, TrdLogLikelihoodRatioElectronProtonTracker_e, TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative, TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive, TrdLogLikelihoodRatioElectronProtonTracker_pion, TrdLogLikelihoodRatioElectronProtonTracker_d;
    tantiproton->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_p);
    telectron->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_e);
    tpass7_negative->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative);
    tpass7_positive->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive);
    tpion->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_pion);
    tdeuteron->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_d);

    Float_t TrdLogLikelihoodRatioProtonHeliumTracker_p, TrdLogLikelihoodRatioProtonHeliumTracker_e, TrdLogLikelihoodRatioProtonHeliumTracker_pass7_negative, TrdLogLikelihoodRatioProtonHeliumTracker_pass7_positive, TrdLogLikelihoodRatioProtonHeliumTracker_pion, TrdLogLikelihoodRatioProtonHeliumTracker_d;
    tantiproton->SetBranchAddress("TrdLogLikelihoodRatioProtonHeliumTracker",&TrdLogLikelihoodRatioProtonHeliumTracker_p);
    telectron->SetBranchAddress("TrdLogLikelihoodRatioProtonHeliumTracker",&TrdLogLikelihoodRatioProtonHeliumTracker_e);
    tpass7_negative->SetBranchAddress("TrdLogLikelihoodRatioProtonHeliumTracker",&TrdLogLikelihoodRatioProtonHeliumTracker_pass7_negative);
    tpass7_positive->SetBranchAddress("TrdLogLikelihoodRatioProtonHeliumTracker",&TrdLogLikelihoodRatioProtonHeliumTracker_pass7_positive);
    tpion->SetBranchAddress("TrdLogLikelihoodRatioProtonHeliumTracker",&TrdLogLikelihoodRatioProtonHeliumTracker_pion);
    tdeuteron->SetBranchAddress("TrdLogLikelihoodRatioProtonHeliumTracker",&TrdLogLikelihoodRatioProtonHeliumTracker_d);

    Float_t Rigidity_p, Rigidity_e, Rigidity_pass7_negative, Rigidity_pass7_positive, Rigidity_pion, Rigidity_d;
    tantiproton->SetBranchAddress("Rigidity",&Rigidity_p);
    telectron->SetBranchAddress("Rigidity",&Rigidity_e);
    tpass7_negative->SetBranchAddress("Rigidity",&Rigidity_pass7_negative);
    tpass7_positive->SetBranchAddress("Rigidity",&Rigidity_pass7_positive);
    tpion->SetBranchAddress("Rigidity",&Rigidity_pion);
    tdeuteron->SetBranchAddress("Rigidity",&Rigidity_d);

    Float_t RichBeta_p, RichBeta_e, RichBeta_pass7_negative, RichBeta_pass7_positive, RichBeta_pion, RichBeta_d;
    tantiproton->SetBranchAddress("RichBeta",&RichBeta_p);
    telectron->SetBranchAddress("RichBeta",&RichBeta_e);
    tpass7_negative->SetBranchAddress("RichBeta",&RichBeta_pass7_negative);
    tpass7_positive->SetBranchAddress("RichBeta",&RichBeta_pass7_positive);
    tpion->SetBranchAddress("RichBeta",&RichBeta_pion);
    tdeuteron->SetBranchAddress("RichBeta",&RichBeta_d);

    Float_t TofBeta_p, TofBeta_e, TofBeta_pass7_negative, TofBeta_pass7_positive, TofBeta_pion, TofBeta_d;
    tantiproton->SetBranchAddress("TofBeta",&TofBeta_p);
    telectron->SetBranchAddress("TofBeta",&TofBeta_e);
    tpass7_negative->SetBranchAddress("TofBeta",&TofBeta_pass7_negative);
    tpass7_positive->SetBranchAddress("TofBeta",&TofBeta_pass7_positive);
    tpion->SetBranchAddress("TofBeta",&TofBeta_pion);
    tdeuteron->SetBranchAddress("TofBeta",&TofBeta_d);

    Double_t Weight_p, Weight_e, Weight_pass7_negative, Weight_pass7_positive, Weight_pion, Weight_d;
    tantiproton->SetBranchAddress("Weight",&Weight_p);
    telectron->SetBranchAddress("Weight",&Weight_e);
    tpass7_negative->SetBranchAddress("Weight",&Weight_pass7_negative);
    tpass7_positive->SetBranchAddress("Weight",&Weight_pass7_positive);
    tpion->SetBranchAddress("Weight",&Weight_pion);
    tdeuteron->SetBranchAddress("Weight",&Weight_d);

    Float_t ProtonCCMVABDT_p, ProtonCCMVABDT_e, ProtonCCMVABDT_pass7_negative, ProtonCCMVABDT_pass7_positive, ProtonCCMVABDT_pion, ProtonCCMVABDT_d;
    tantiproton->SetBranchAddress("ProtonCCMVABDT",&ProtonCCMVABDT_p);
    telectron->SetBranchAddress("ProtonCCMVABDT",&ProtonCCMVABDT_e);
    tpass7_negative->SetBranchAddress("ProtonCCMVABDT",&ProtonCCMVABDT_pass7_negative);
    tpass7_positive->SetBranchAddress("ProtonCCMVABDT",&ProtonCCMVABDT_pass7_positive);
    tpion->SetBranchAddress("ProtonCCMVABDT",&ProtonCCMVABDT_pion);
    tdeuteron->SetBranchAddress("ProtonCCMVABDT",&ProtonCCMVABDT_d);

    Float_t EcalBDT_EnergyD_p, EcalBDT_EnergyD_e, EcalBDT_EnergyD_pass7_negative, EcalBDT_EnergyD_pass7_positive, EcalBDT_EnergyD_pion, EcalBDT_EnergyD_d;
    tantiproton->SetBranchAddress("EcalBDT_EnergyD",&EcalBDT_EnergyD_p);
    telectron->SetBranchAddress("EcalBDT_EnergyD",&EcalBDT_EnergyD_e);
    tpass7_negative->SetBranchAddress("EcalBDT_EnergyD",&EcalBDT_EnergyD_pass7_negative);
    tpass7_positive->SetBranchAddress("EcalBDT_EnergyD",&EcalBDT_EnergyD_pass7_positive);
    tpion->SetBranchAddress("EcalBDT_EnergyD",&EcalBDT_EnergyD_pion);
    tdeuteron->SetBranchAddress("EcalBDT_EnergyD",&EcalBDT_EnergyD_d);


    //float Richbetacut = atof(argv[1]);
    //string Richbetacut[16] = {"0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95"};
    //for (float Richbetacut = 0.990; Richbetacut<1.004; Richbetacut = Richbetacut + 0.001)
    //  {
    // Settings for mass^2 fit, still tuning.
    //float Richbetacut = 100;
    //float TrdLikelihoodcut = 0.9;
    //float TrdLikelihoodHeProtoncut = 0.01;
    //float TofBetacut = 0.94;
    //float TofBetacut = 100;
    //float ProtonCCMVABDcutT = 0.99;
    //float EcalBDT_EnergyDcut = -0.9;


    // Settings for tof fit
    float TrdLikelihoodcut = 0.8;
    float TrdLikelihoodHeProtoncut = 0.2;
    float TofBetacut = 100;
    float ProtonCCMVABDcutT = 0.9;
    float EcalBDT_EnergyDcut = -0.9;
    float Richbetacut = 0.997;

    // TRD_electron_proton: all, TRD_He_Proton: Antiproton Template, positivedata,  
        int entries = tantiproton->GetEntries();
        for( int i=0;i<entries;i++){
             tantiproton->GetEntry(i);
             if ( TrdLogLikelihoodRatioElectronProtonTracker_p > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_p < TrdLikelihoodHeProtoncut ){
//             template_antiproton->Fill((pow(Rigidity_p,2) * (1-pow(TofBeta_p,2))) / pow(TofBeta_p,2), Weight_p);
             Tof_template_antiproton->Fill(1./TofBeta_p);
             }
        }
        int entries2 = telectron->GetEntries();
        for( int i=0;i<entries2;i++){
             telectron->GetEntry(i);
             if (TrdLogLikelihoodRatioElectronProtonTracker_e > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_e < TrdLikelihoodHeProtoncut  ){
//             template_electron->Fill((pow(Rigidity_e,2) * (1-pow(TofBeta_e,2))) / pow(TofBeta_e,2), Weight_e);
             Tof_template_electron->Fill(1./TofBeta_e);
             }
        }
        int entries3 = tpass7_negative->GetEntries();
        for( int i=0;i<entries3;i++){
             tpass7_negative->GetEntry(i);
//             if (TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_pass7_negative < TrdLikelihoodHeProtoncut && TofBeta_pass7_negative < TofBetacut && ProtonCCMVABDT_pass7_negative>ProtonCCMVABDcutT && EcalBDT_EnergyD_pass7_negative < EcalBDT_EnergyDcut && ( RichBeta_pass7_negative > Richbetacut || RichBeta_pass7_negative == 0 )){
             if (TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_pass7_negative < TrdLikelihoodHeProtoncut && TofBeta_pass7_negative < TofBetacut && ProtonCCMVABDT_pass7_negative>ProtonCCMVABDcutT && EcalBDT_EnergyD_pass7_negative < EcalBDT_EnergyDcut && ( RichBeta_pass7_negative < Richbetacut )){
//             data_pass7_negative->Fill((pow(Rigidity_pass7_negative,2) * (1-pow(TofBeta_pass7_negative,2))) / pow(TofBeta_pass7_negative,2));
             Tof_data_pass7_negative->Fill(1./TofBeta_pass7_negative);
             }
        }
        int entries4 = tpass7_positive->GetEntries();
        for( int i=0;i<entries4;i++){
             tpass7_positive->GetEntry(i);
//             if (TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_pass7_positive < TrdLikelihoodHeProtoncut && TofBeta_pass7_positive < TofBetacut && ProtonCCMVABDT_pass7_positive>ProtonCCMVABDcutT && EcalBDT_EnergyD_pass7_positive < EcalBDT_EnergyDcut && ( RichBeta_pass7_positive > Richbetacut || RichBeta_pass7_positive == 0 ) ){
             if (TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_pass7_positive < TrdLikelihoodHeProtoncut && TofBeta_pass7_positive < TofBetacut && ProtonCCMVABDT_pass7_positive>ProtonCCMVABDcutT && EcalBDT_EnergyD_pass7_positive < EcalBDT_EnergyDcut && ( RichBeta_pass7_positive < Richbetacut ) ){
//             data_pass7_positive->Fill((pow(Rigidity_pass7_positive,2) * (1-pow(TofBeta_pass7_positive,2))) / pow(TofBeta_pass7_positive,2));
             Tof_data_pass7_positive->Fill(1./TofBeta_pass7_positive);
             }
        }
        int entries5 = tpion->GetEntries();
        for( int i=0;i<entries5;i++){
             tpion->GetEntry(i);
//             if (TrdLogLikelihoodRatioElectronProtonTracker_pion > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_pion < TrdLikelihoodHeProtoncut ){
//             template_pion->Fill((pow(Rigidity_pion,2) * (1-pow(TofBeta_pion,2))) / pow(TofBeta_pion,2),Weight_pion);
             Tof_template_pion->Fill(1./TofBeta_pion);
//             }
        }

/*
        int entries6 = tdeuteron->GetEntries();
        for( int i=0;i<entries6;i++){
             tdeuteron->GetEntry(i);
             if (TrdLogLikelihoodRatioElectronProtonTracker_d > TrdLikelihoodcut && TrdLogLikelihoodRatioProtonHeliumTracker_d < TrdLikelihoodHeProtoncut ){
             template_deuteron->Fill((pow(Rigidity_d,2) * (1-pow(TofBeta_d,2))) / pow(TofBeta_d,2),Weight_d);
             }
        }
*/

        //// Parametrisation
/*
        //// 1. Mass^2 fit: get antiproton from antiproton template, remove pions
        template_antiproton->Fit("gaus","","",0.5,1.2596); //0.938 for proton, mass^2=0.8798. Symmetric bin:0.5-1.2596. If Fit range is not given, then defalt is whole range.
        TF1 *antiprotonfit = template_antiproton->GetFunction("gaus");
        Utilities::ConvertToHistogram(antiprotonfit, *template_antiproton_parametrized);
        //// 2. Mass^2 fit: get deuteron from deuteron template, remove proton
        template_deuteron->Fit("gaus","","",2);  //1.875 for deuteron, mass^2=3.51. Only low limit is given.  If Fit range is not given, then defalt is whole range.
        TF1 *deuteronfit = template_deuteron->GetFunction("gaus");
        Utilities::ConvertToHistogram(deuteronfit, *template_deuteron_parametrized);
        //// 3. Mass^2 fit: get backround template from negative data with gaussion template for pion and electrons
        data_pass7_negative->Fit("gaus","","",-1.3,1.3);
        TF1 *data_pass7_negative_fit = data_pass7_negative->GetFunction("gaus");
        Utilities::ConvertToHistogram(data_pass7_negative_fit, *data_pass7_negative_fit_parametrized);
        //// 4. Mass^2 fit: get antiproton template from positive data.  
        data_pass7_positive->Fit("gaus","","",0.5,1.2596);  // 0.938 for proton, mass^2=0.8798. Symmetric bin:0.5-1.2596.
        TF1 *data_pass7_positive_fit = data_pass7_positive->GetFunction("gaus");
        Utilities::ConvertToHistogram(data_pass7_positive_fit, *data_pass7_positive_fit_parametrized);
*/

        //// 5. Tof fit: get antiproton tof_beta from antiproton MC, remove pion.
        Tof_template_antiproton->Fit("gaus","","",1.0,1.3); // to be tuned. 
        TF1 *Tof_template_antiprotonfit = Tof_template_antiproton->GetFunction("gaus");
        Utilities::ConvertToHistogram(Tof_template_antiprotonfit,*Tof_template_antiproton_parametrized);
        //// 6. Tof fit: get backround template from negative data with gaussion template for pion and electrons
        Tof_data_pass7_negative->Fit("gaus","","",0.8,1.1); // to be tuned. 
        TF1 *Tof_data_pass7_negative_fit = Tof_data_pass7_negative->GetFunction("gaus");
        Utilities::ConvertToHistogram(Tof_data_pass7_negative_fit, *Tof_data_pass7_negative_fit_parametrized); 
        //// 7. Tof fit: get antiproton tof_beta from positive data.
        Tof_data_pass7_positive->Fit("gaus","","",0.9,1.3); // to be tuned.
        TF1 *Tof_data_pass7_positive_fit = Tof_data_pass7_positive->GetFunction("gaus");
        Utilities::ConvertToHistogram(Tof_data_pass7_positive_fit, *Tof_data_pass7_positive_fit_parametrized);


        //// Template Fit
/*
        //// 1. Template Fit Negaitve
        MYUtilities::TemplateFitter templateFitter_negative(0);
    //    templateFitter_negative.AddTemplateHistogram(template_antiproton);
    //    templateFitter_negative.AddTemplateHistogram(template_antiproton_parametrized);
    //    templateFitter_negative.AddTemplateHistogram(data_pass7_positive_fit_parametrized);
        templateFitter_negative.AddTemplateHistogram(data_pass7_positive);
        templateFitter_negative.AddTemplateHistogram(data_pass7_negative_fit_parametrized);
    //    templateFitter_negative.AddTemplateHistogram(template_electron);
    //    templateFitter_negative.AddTemplateHistogram(template_pion);
        templateFitter_negative.SetDataHistogram(data_pass7_negative);
        templateFitter_negative.SetStartValue(0, 0.00002); //Start Values are relative! Must be between 0 and 1.
        templateFitter_negative.SetStartValue(1, 0.5);
    //    templateFitter_negative.SetStartValue(2, 0.5);
        templateFitter_negative.Fit(1); // if this fit range is smaller than whole range, then two vertical lines are drawed in fit figure.  If Fit range is not given, then defalt is whole range.
        //// 2. Template Fit Positive
        MYUtilities::TemplateFitter templateFitter_positive(0);
        templateFitter_positive.AddTemplateHistogram(template_antiproton);
    //    templateFitter_positive.AddTemplateHistogram(template_antiproton_parametrized);
        templateFitter_positive.AddTemplateHistogram(template_pion);
        templateFitter_positive.AddTemplateHistogram(template_deuteron_parametrized);
        templateFitter_positive.SetDataHistogram(data_pass7_positive);
        templateFitter_positive.Fit(1,-1.3,5.02); //If Fit range is not given, then defalt is whole range.
*/
        //// 3. 1./TOF Beta
        MYUtilities::TemplateFitter templateFitter_tof(0);
//        templateFitter_tof.AddTemplateHistogram(Tof_template_antiproton_parametrized);
//        templateFitter_tof.AddTemplateHistogram(Tof_template_antiproton);
        templateFitter_tof.AddTemplateHistogram(Tof_data_pass7_positive);
//        templateFitter_tof.AddTemplateHistogram(Tof_data_pass7_negative_fit_parametrized);
        templateFitter_tof.AddTemplateHistogram(Tof_template_electron);
        templateFitter_tof.AddTemplateHistogram(Tof_template_pion);
        templateFitter_tof.SetDataHistogram(Tof_data_pass7_negative);
        templateFitter_tof.Fit(1,0.8,1.4); //If Fit range is not given, then defalt is whole range.


         // plots
/*
         // 1. Mass2
         TCanvas * c1 = new TCanvas;
         template_antiproton->Draw("");
         gPad->SetLogy();
         c1->SaveAs((string("Time_Averaged_ratio_Low/mass2/template_antiproton_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c2 = new TCanvas;
         template_electron->Draw("");
         gPad->SetLogy();
         c2->SaveAs( (string("Time_Averaged_ratio_Low/mass2/template_electron_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str() );

         TCanvas * c21 = new TCanvas;
         template_pion->Draw("");
         gPad->SetLogy();
         c21->SaveAs( (string("Time_Averaged_ratio_Low/mass2/template_pion_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str() );

         TCanvas * c3 = new TCanvas;
         data_pass7_negative->Draw("");
         gPad->SetLogy();
         c3->SaveAs( (string("Time_Averaged_ratio_Low/mass2/negativedata_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str()); 

         TCanvas * c31 = new TCanvas;
         data_pass7_positive->Draw("");
         gPad->SetLogy();
         c31->SaveAs( (string("Time_Averaged_ratio_Low/mass2/positivedata_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c4 = templateFitter_negative.CreateResultDrawing("Fit_Result",800,600);
         c4->SaveAs( (string("Time_Averaged_ratio_Low/mass2/FitResult_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c5 = templateFitter_negative.CreateContourPlot(0,1,"ContourPlot",800,600);
         c5->SaveAs( (string("Time_Averaged_ratio_Low/mass2/ContourPlot_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c6 = templateFitter_positive.CreateResultDrawing("Fit_Result_positive",800,600);
         c6->SaveAs( (string("Time_Averaged_ratio_Low/mass2/FitResult_positive_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c7 = templateFitter_positive.CreateContourPlot(0,1,"ContourPlot_positive",800,600);
         c7->SaveAs( (string("Time_Averaged_ratio_Low/mass2/ContourPlot_positive_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());
*/

         //// 2. Tof 
         TCanvas * c8 = templateFitter_tof.CreateResultDrawing("Fit_Result_tof",800,600);
         c8->SaveAs( (string("Time_Averaged_ratio_Low/tof/FitResult_tof_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c9 = new TCanvas;
         Tof_template_electron->Draw("");
         gPad->SetLogy();
         c9->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_template_electron_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c10 = new TCanvas;
         Tof_template_antiproton_parametrized->Draw("");
         gPad->SetLogy();
         c10->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_template_antiproton_parametrized_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c101 = new TCanvas;
         Tof_template_antiproton->Draw("");
         gPad->SetLogy();
         c101->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_template_antiproton_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c11 = new TCanvas;
         Tof_template_pion->Draw("");
         gPad->SetLogy();
         c11->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_template_pion_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c12 = new TCanvas;
         Tof_data_pass7_negative->Draw("");
         gPad->SetLogy();
         c12->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_data_pass7_negative_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c13 = new TCanvas;
         Tof_data_pass7_negative_fit_parametrized->Draw("");
         gPad->SetLogy();
         c13->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_data_pass7_negative_fit_parametrized_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c14 = new TCanvas;
         Tof_data_pass7_positive->Draw("");
         gPad->SetLogy();
         c14->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_data_pass7_positive_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

         TCanvas * c15 = new TCanvas;
         Tof_data_pass7_positive_fit_parametrized->Draw("");
         gPad->SetLogy();
         c15->SaveAs( (string("Time_Averaged_ratio_Low/tof/Tof_data_pass7_positive_fit_parametrized_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".png")).c_str());

        // store those
/*
        // Mass2 fit negative
        vector<double> result, ResultError;
        assert(result.empty());
        assert(ResultError.empty());
        result.assign(templateFitter_negative.GetAbsoluteResult().begin(), templateFitter_negative.GetAbsoluteResult().end());
        ResultError.assign(templateFitter_negative.GetAbsoluteResultError().begin(), templateFitter_negative.GetAbsoluteResultError().end());
        double Chi2 = templateFitter_negative.Chi2(); 
        int NDF = templateFitter_negative.NDF();
        double CHI2dof = Chi2/NDF;

        // Mass2 fit positive
        vector<double> result_positive, ResultError_positive;
        assert(result_positive.empty());
        assert(ResultError_positive.empty());
        result_positive.assign(templateFitter_positive.GetAbsoluteResult().begin(), templateFitter_positive.GetAbsoluteResult().end());
        ResultError_positive.assign(templateFitter_positive.GetAbsoluteResultError().begin(), templateFitter_positive.GetAbsoluteResultError().end());
*/

        vector<double> result_tof, ResultError_tof;
        assert(result_tof.empty());
        assert(ResultError_tof.empty());
        result_tof.assign(templateFitter_tof.GetAbsoluteResult().begin(), templateFitter_tof.GetAbsoluteResult().end());
        ResultError_tof.assign(templateFitter_tof.GetAbsoluteResultError().begin(), templateFitter_tof.GetAbsoluteResultError().end());
        double Chi2_tof = templateFitter_tof.Chi2();
        int NDF_tof = templateFitter_tof.NDF();
        double CHI2dof_tof = Chi2_tof/NDF_tof;

 
        ///////// Choose ratio you use according to different cuts  //////////////// 
/*
        // Mass2 fit ratio       
//      v_ratio_result.push_back(result[0]/result_positive[0]*100000); proton from template fit
        v_ratio_result.push_back(result[0]/data_pass7_positive->GetEntries()*100000); // proton from total number after cut
//        v_ratio_result.push_back(result[0]/tpass7_positive->GetEntries()*100000); // proton from total number before cut
        v_ratio_result_with_effective.push_back(v_ratio_result.back()*Effective_Acceptance->GetY()[index]);   
        v_chi2.push_back(CHI2dof);
*/
        // Tof Fit ratio
//      v_ratio_result_tof.push_back(result_tof[0]/result_positive[0]*100000); proton from template fit
//      v_ratio_result_tof.push_back(result_tof[0]/data_pass7_positive->GetEntries()*100000); // tof fit, proton from total number after cut
        v_ratio_result_tof.push_back(result_tof[0]/tpass7_positive->GetEntries()*100000); // proton from total number before cut
        v_ratio_result_tof_with_effective.push_back(v_ratio_result_tof.back()*Effective_Acceptance->GetY()[index]);
        v_chi2_tof.push_back(CHI2dof_tof);
         
//        cout << "fit result:" << result << endl;
//        cout << "ResultError:" << ResultError << endl;
//        cout << "Chi2:" << Chi2 << endl;
//        cout << "NDF:" << NDF << endl;
//        cout << "CHI2dof:" <<CHI2dof << endl;
//        cout << "proton total number:" << data_pass7_positive->GetEntries() <<endl;
//        cout << "proton number from template fit:" << result_positive[0] <<endl;
//        cout << "ratio with proton from template fit:" << result[0]/result_positive[0] <<endl; // proton from template fit 
//        cout << "error with proton from template fit:" << ResultError[0]/result_positive[0] <<endl; // proton from template fit
//        cout << "ratio with total proton (positive data after cut):" << result[0]/data_pass7_positive->GetEntries() <<endl; // proton from total number after cut
//        cout << "error with total proton (positive data after cut):" << ResultError[0]/data_pass7_positive->GetEntries() <<endl; // proton from total number after cut
//        cout << "ratio with total proton (positive data before cut):" << result[0]/tpass7_positive->GetEntries() <<endl; // proton from total number before cut
//        cout << "error with total proton (positive data before cut):" << ResultError[0]/tpass7_positive->GetEntries() <<endl; // proton from total number before cut


//        cout << "raw negative events number:" << tpass7_negative->GetEntries() <<endl;
//        cout << "after cut negative events number:" << Tof_data_pass7_negative->GetEntries() <<endl;
        cout << "tof fit result:" << result_tof << endl;
//        cout << "proton total number from template fit:" << result_positive[0] << endl; // tof fit, proton from template fit   
        cout << "proton total number before cut:" << tpass7_positive->GetEntries()  << endl; // tof fit, proton from total number before cut
//        cout << "proton total number after cut:" << data_pass7_positive->GetEntries()  << endl; // tof fit, proton from total number after cut
        cout << "Chi2/dof_tof:" << CHI2dof_tof << endl;
        cout << "ratio from 1./tofbeta fit: "<< result_tof[0]/tpass7_positive->GetEntries()*100000  <<endl;


        TFile testroot((string("Time_Averaged_ratio_Low/averaged_ratio_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".root")).c_str(),"RECREATE");
/*
        // Mass2
        template_antiproton->Write("template_antiproton");
        template_antiproton_parametrized->Write("template_antiproton_parametrized");
        template_electron->Write("template_electron");
        template_pion->Write("template_pion");
        data_pass7_negative->Write("data_pass7_negative");
        data_pass7_positive->Write("data_pass7_positive");
        data_pass7_positive_fit_parametrized->Write("data_pass7_positive_fit_parametrized");
        template_deuteron->Write("template_deuteron");
        template_deuteron_parametrized->Write("template_deuteron_parametrized");
        //      trdhe_template_antiproton->Write("trdhe_template_antiproton");
        //      trdhe_template_electron->Write("trdhe_template_electron");
        //      trdhe_data_pass7_negative->Write("trdhe_data_pass7_negative");
        //      trdhe_data_pass7_positive->Write("trdhe_data_pass7_positive");
*/
        // Tof
        Tof_template_antiproton->Write("Tof_template_antiproton");
        Tof_template_antiproton_parametrized->Write("Tof_template_antiproton_parametrized");
        Tof_template_electron->Write("Tof_template_electron");
        Tof_template_pion->Write("Tof_template_pion");
        Tof_data_pass7_negative->Write("Tof_data_pass7_negative");
        Tof_data_pass7_negative_fit_parametrized->Write("Tof_data_pass7_negative_fit_parametrized");
        Tof_data_pass7_positive->Write("Tof_data_pass7_positive");
        Tof_data_pass7_positive_fit_parametrized->Write("Tof_data_pass7_positive_fit_parametrized");
        testroot.Close();

}


////////////////////////////////////////////////////////////////
///////////////       Results  Plots        ////////////////////
////////////////////////////////////////////////////////////////

auto substart = bincenter.begin() + stoi(rigidity_start);
auto subend = bincenter.begin() + stoi(rigidity_end);
std::vector<double> subbincenter(substart, subend);

auto lowstart = bincenter.begin()+1;  // first bin is used for underflow. 
auto lowend = bincenter.begin()+13;
std::vector<double> lowbincenter(lowstart, lowend);

auto first = AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin();
auto end = AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin()+12;
std::vector<double> published(first, end);
auto firsterror = AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin();
auto enderror = AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin()+12;
std::vector<double> publishederror(firsterror, enderror);

for (auto i=0; i<published.size(); i++){
published.at(i) = published.at(i)*100000;
publishederror.at(i) = publishederror.at(i)*100000;
}

cout<< lowbincenter <<endl;
cout<< published <<endl;

/*
// Mass2
TGraph *g_ratio = new TGraph(v_ratio_result.size(), subbincenter.data(), v_ratio_result.data());
TGraph *g_ratio_with_effective = new TGraph(v_ratio_result.size(), subbincenter.data(), v_ratio_result_with_effective.data());
TGraph *g_chi2 = new TGraph(v_chi2.size(), subbincenter.data(), v_chi2.data());
*/

// Tof
TGraph *g_ratio_tof = new TGraph(v_ratio_result_tof.size(), subbincenter.data(), v_ratio_result_tof.data());
TGraph *g_ratio_tof_with_effective = new TGraph(v_ratio_result_tof.size(), subbincenter.data(), v_ratio_result_tof_with_effective.data());
TGraph *g_chi2_tof = new TGraph(v_chi2_tof.size(), subbincenter.data(), v_chi2_tof.data());

TGraphErrors *g_ratio_published = new TGraphErrors(12, lowbincenter.data(), published.data(), 0, publishederror.data());


chdir("./plots");
/*
TCanvas cc1("cc1","cc1",1000,500);
g_ratio->SetMarkerStyle(15);
g_ratio->SetMarkerColor(2);
g_ratio_published->Draw("AP *");
g_ratio->Draw("same P");
g_ratio_published->SetTitle("");
g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");
cc1.Update();
cc1.SaveAs("low_ratio.pdf");

TCanvas cc11("cc11","cc11",1000,500);
g_ratio_with_effective->SetMarkerStyle(15);
g_ratio_with_effective->SetMarkerColor(2);
g_ratio_published->Draw("AP *");
g_ratio_with_effective->Draw("same P");
g_ratio_published->SetTitle("");
g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");
cc11.Update();
cc11.SaveAs("low_ratio_with_effective.pdf");

TCanvas cc2("cc2","cc2",1000,500);
g_chi2->Draw("AP *");
g_chi2->SetTitle("");
g_chi2->GetXaxis()->SetTitle("Rigidity (GV)");
g_chi2->GetYaxis()->SetTitle("Chi2/dof");
cc2.Update();
cc2.SaveAs("chi2.pdf");
*/
TCanvas cc3("cc3","cc3",1000,500);
g_ratio_tof->SetMarkerStyle(15);
g_ratio_tof->SetMarkerColor(2);
g_ratio_published->Draw("AP *");
g_ratio_tof->Draw("same P");
g_ratio_published->SetTitle("");
g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");
cc3.Update();
cc3.SaveAs("low_ratio_tof.pdf");

TCanvas cc33("cc33","cc33",1000,500);
g_ratio_tof_with_effective->SetMarkerStyle(15);
g_ratio_tof_with_effective->SetMarkerColor(2);
g_ratio_published->Draw("AP *");
g_ratio_tof_with_effective->Draw("same P");
g_ratio_published->SetTitle("");
g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");
cc33.Update();
cc33.SaveAs("low_ratio_tof_with_effective.pdf");

TCanvas cc4("cc4","cc4",1000,500);
g_chi2_tof->Draw("AP *");
g_chi2_tof->SetTitle("");
g_chi2_tof->GetXaxis()->SetTitle("Rigidity (GV)");
g_chi2_tof->GetYaxis()->SetTitle("Chi2/dof");
cc4.Update();
cc4.SaveAs("chi2_tof.pdf");

TCanvas cc5("cc5","cc5",1000,500);
Effective_Acceptance->RemovePoint(0); //    MC1042 is from 1 to 800, so delete first point 0.8-1.0.
Effective_Acceptance->RemovePoint(61); //    MC1042 is from 1 to 800, so delete first point 0.8-1.0.
Effective_Acceptance->Draw("AP *");
Effective_Acceptance->SetTitle("");
Effective_Acceptance->GetXaxis()->SetTitle("Rigidity (GV)");
Effective_Acceptance->GetYaxis()->SetTitle("Effective_Acceptance");
Effective_Acceptance->GetXaxis()->SetLimits(1.0,800);
gPad->SetLogx();
cc5.Update();
cc5.SaveAs("Effective_Acceptance.pdf");


TFile ratioroot((string("Ratio")+ string("_") +issversion + string(".root")).c_str(),"RECREATE");
Effective_Acceptance->Write("Effective_Acceptance");
//g_ratio->Write("ratio_mass2");
//g_ratio_with_effective->Write("ratio_mass2_with_effective");
//g_chi2->Write("chi2_mass2");
g_ratio_tof->Write("ratio_tof");
g_ratio_tof_with_effective->Write("ratio_tof_with_effective");
g_chi2_tof->Write("chi2_tof");
g_ratio_published->Write("published");
ratioroot.Close();

//  }

  return EXIT_SUCCESS;
}
