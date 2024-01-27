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

chdir( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/")).c_str());
//chdir( (string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v4.4") + string("/totalall/")).c_str());

////////////////////////////////////////////////////
vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());
vector<double> bincenter (AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins());
vector<double> v_ratio_result_tof, v_chi2_tof, v_ratio_result_tof_with_effective, v_error, v_error_relative;
vector<double> subbincenter(bincenter.begin()+stoi(rigidity_start), bincenter.begin()+stoi(rigidity_end));
vector<double> subbinedge(binning.begin()+stoi(rigidity_start), binning.begin()+stoi(rigidity_end)+1);

// Published result
int publishedresult_showpoint = 21;  //16points:1-17:1.08-5.635.  29points:1.08-17.3.  24points:1.08-11.5
std::vector<double> lowbincenter(bincenter.begin()+1, bincenter.begin()+publishedresult_showpoint+1); //first bin is used for underflow. 
std::vector<double> published(AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin()+publishedresult_showpoint);
std::vector<double> publishederror(AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin()+publishedresult_showpoint);
vector<double> v_StatisticError_relative(AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.begin()+stoi(rigidity_start), AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.begin()+stoi(rigidity_end));
vector<double> v_SystematicError_relative(AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.begin()+stoi(rigidity_start), AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.begin()+stoi(rigidity_end));
for (auto i=0; i<published.size(); i++){
  published.at(i) = published.at(i)*100000;
  publishederror.at(i) = publishederror.at(i)*100000;
}
TGraphErrors *g_ratio_published = new TGraphErrors(publishedresult_showpoint, lowbincenter.data(), published.data(), 0, publishederror.data());
TGraph *g_Published_Statistic_Error_Relative = new TGraph(subbincenter.size(), subbincenter.data(), v_StatisticError_relative.data());
TH1D h_Published_Statistic_Error_Relative = TH1D("", "", subbincenter.size(), subbinedge.data());
Utilities::ConvertToHistogram(g_Published_Statistic_Error_Relative, h_Published_Statistic_Error_Relative);
TGraph *g_Published_Systematic_Error_Relative = new TGraph(subbincenter.size(), subbincenter.data(), v_SystematicError_relative.data());
TH1D h_Published_Systematic_Error_Relative = TH1D("", "", subbincenter.size(), subbinedge.data());
Utilities::ConvertToHistogram(g_Published_Systematic_Error_Relative, h_Published_Systematic_Error_Relative);

//Intermediate Range result
string IntermediateFilename;
if (issversion == "pass7.8"){
    IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/binmerge1/intermediate_0124_free_pass7.8binmerge1.root");
}
else if (issversion == "2016paper"){
    IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/binmerge1/intermediate_0124_free_2016paperbinmerge1.root");
}
TFile *IntermediateFile = new TFile(IntermediateFilename.c_str());
TGraphErrors *g_IntermediateResult = (TGraphErrors*)IntermediateFile->Get("g_ratio_with_effective_acceptance");
for (int index = 0; index < g_IntermediateResult->GetN(); index++){
g_IntermediateResult->SetPoint(index, g_IntermediateResult->GetX()[index], g_IntermediateResult->GetY()[index]*100000);
g_IntermediateResult->SetPointError(index, 0, g_IntermediateResult->GetErrorY(index)*100000);
}

// Effective acceptance
TFile *f_eff_antiproton = new TFile("EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root");
TFile *f_eff_proton = new TFile("EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root");
TGraphAsymmErrors *Acceptance_antiproton = (TGraphAsymmErrors*)f_eff_antiproton->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
TGraphAsymmErrors *Acceptance_proton = (TGraphAsymmErrors*)f_eff_proton->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
int number =17; // To 5.9GV
Double_t xa[number],Aa[number],ea[number],xp[number],Ap[number],ep[number];
TGraphErrors *Effective_Acceptance = new TGraphErrors(number);
for (int p = 0; p < number; p++) {       
  Acceptance_antiproton->GetPoint(p,xa[p],Aa[p]);
  Acceptance_proton->GetPoint(p,xp[p],Ap[p]);
  ea[p] = Acceptance_antiproton->GetErrorY(p);
  ep[p] = Acceptance_proton->GetErrorY(p);
  Effective_Acceptance->SetPoint(p,xa[p],Ap[p]/Aa[p]);
  Effective_Acceptance->SetPointError(p, 0.0, sqrt(pow(ep[p],2)/pow(Aa[p],2)+pow(Ap[p],2)/pow(Aa[p],4)*pow(ea[p],2)));
  //cout << p << endl;
  //cout << "x:" << xa[p] << endl;
  //cout << "ratio of acceptance:" << Ap[p]/Aa[p] << endl;
  //cout << "Ap[p]" << Ap[p] << endl;
  //cout << "Aa[p]" << Aa[p] << endl;
}

// Load root files and apply cuts
for (int index=stoi(rigidity_start); index<stoi(rigidity_end); index++)
{   // rigidity loop start
    cout<< string("index now is:") + to_string(index) <<endl;
    string negativename;
    string positivename;
    string antiprotonname = string("B1042_antipr.pl1.1800_7.6_all_Tree_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    string electronname = string("B1091_el.pl1.0_25200_7.6_all_Tree_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    string pionname = string("B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    string deuteronname = string("B1128_d.pl1ph.021000_7.6_all_Tree_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    if (issversion == "pass7.8"){
        negativename = string("B1130_pass7_7.8_all_Tree_negative_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
        positivename = string("B1130_pass7_7.8_all_Tree_positive_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    }
    else if (issversion == "2016paper"){
        negativename = string("B1130_pass7_7.8_all_Tree_negative_May2015_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
        positivename = string("B1130_pass7_7.8_all_Tree_positive_May2015_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string(".root");
    }

   // Allcut:     TRD_electron_proton: all, TRD_He_Proton: Antiproton Template, positivedata,
   // float TrdLikelihoodcut = 0.8;
   // float TrdLikelihoodHeProtoncut = 0.2;
   // float TofBetacut = 100;
   // float ProtonCCMVABDcut = 0.9;
   // float EcalBDT_EnergyDcut = -0.9;
   // float Richbetacut = 1.01;
    std::string trackerpattern = "0";
    TCut PatternCut = (string("Pattern==") + trackerpattern).c_str();
    TCut TrdLikelihoodCut = "TrdLogLikelihoodRatioElectronProtonTracker > 0.8";    // Antiproton:1to1.2, Electron:0.5
    TCut TrdLikelihoodHeProtonCut = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.1";
    TCut ProtonCCMVABDCut = "ProtonCCMVABDT > 0.9";
    TCut EcalBDT_EnergyDCut = "EcalBDT_EnergyD < -0.9";
    TCut RicheBetaCut = "RichBeta < 0.997";
    TCut TofBetaCut = "TofBeta < 100";

    TCut AntiprotonCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut; 
    TCut ElectronCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut;
    //TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RicheBetaCut && TofBetaCut;
    TCut NegativeCut = "TrdLogLikelihoodRatioElectronProtonTracker>0.75*TofBeta && RichIsNaF==1";
    //TCut PositiveCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RicheBetaCut && TofBetaCut;
    TCut PositiveCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut;
    TCut PionCut; 
    TCut DeuteronCut;
    TCut AllCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RicheBetaCut;

    TCut ISSNegativeTemplate = "TrdLogLikelihoodRatioElectronProtonTracker<0.75*TofBeta";
    // apply cut
    TFile *fantiproton = new TFile(antiprotonname.c_str());
    TTree *tantiproton_all = (TTree*)fantiproton->Get("AntiprotonLowEnergyTree");
    gROOT->cd();
    TTree *tantiproton = tantiproton_all->CopyTree(AntiprotonCut);
    fantiproton->Close();
    
    TFile *felectron = new TFile(electronname.c_str());
    TTree *telectron_all = (TTree*)felectron->Get("AntiprotonLowEnergyTree");
    gROOT->cd();
    TTree *telectron = telectron_all->CopyTree(ElectronCut);
    felectron->Close();

    TFile *fpass7_negative = new TFile(negativename.c_str());
    TTree *tpass7_negative_all = (TTree*)fpass7_negative->Get("AntiprotonLowEnergyTree");
    gROOT->cd();
    TTree *tpass7_negative = tpass7_negative_all->CopyTree(NegativeCut);
    gROOT->cd();
    TTree *tpass7_negative_template = tpass7_negative_all->CopyTree(ISSNegativeTemplate);
    fpass7_negative->Close();

    TFile *fpass7_positive = new TFile(positivename.c_str());
    TTree *tpass7_positive_all = (TTree*)fpass7_positive->Get("AntiprotonLowEnergyTree");
    gROOT->cd();
    TTree *tpass7_positive = tpass7_positive_all->CopyTree(PositiveCut);
    //TTree *tpass7_positive = tpass7_positive_all->CopyTree("");
    fpass7_positive->Close();

    TFile *fpion = new TFile(pionname.c_str());
    TTree *tpion_all = (TTree*)fpion->Get("AntiprotonLowEnergyTree");
    gROOT->cd();
    TTree *tpion = tpion_all->CopyTree(PionCut);
    fpion->Close();

    TFile *fdeuteron = new TFile(deuteronname.c_str());
    TTree *tdeuteron_all = (TTree*)fdeuteron->Get("AntiprotonLowEnergyTree");
    gROOT->cd();
    TTree *tdeuteron = tdeuteron_all->CopyTree(DeuteronCut);
    fdeuteron->Close();


    // Fill histograms
    tantiproton->Draw("1./TofBeta>>th1f_antiproton(30,0.6,1.5)");
    TH1F *Tof_template_antiproton = (TH1F*)gDirectory->Get("th1f_antiproton");
    Tof_template_antiproton->SetTitle("Antiproton");

    telectron->Draw("1./TofBeta>>th1f_electron(30,0.6,1.5)");
    TH1F *Tof_template_electron = (TH1F*)gDirectory->Get("th1f_electron");
    Tof_template_electron->SetTitle("Electron");

    tpass7_negative->Draw("1./TofBeta>>th1f_negative(30,0.6,1.5)");
    TH1F *Tof_data_pass7_negative = (TH1F*)gDirectory->Get("th1f_negative");
    Tof_data_pass7_negative->SetTitle("ISSnegative");

    tpass7_positive->Draw("1./TofBeta>>th1f_positive(30,0.6,1.5)");
    TH1F *Tof_data_pass7_positive = (TH1F*)gDirectory->Get("th1f_positive");
    Tof_data_pass7_positive->SetTitle("ISSpositive");

    tpion->Draw("1./TofBeta>>th1f_pion(30,0.6,1.5)");
    TH1F *Tof_template_pion = (TH1F*)gDirectory->Get("th1f_pion");
    Tof_template_pion->SetTitle("Pion");

    tdeuteron->Draw("1./TofBeta>>th1f_deuteron(30,0.6,1.5)");
    TH1F *Tof_template_deuteron = (TH1F*)gDirectory->Get("th1f_deuteron");
    Tof_template_deuteron->SetTitle("Deuteron");

    tpass7_negative_template->Draw("1./TofBeta>>th1f_negative_template(30,0.6,1.5)");
    TH1F *Tof_data_pass7_negative_template = (TH1F*)gDirectory->Get("th1f_negative_template");
    Tof_data_pass7_negative_template->SetTitle("ISSnegative_template");

    // Tof_Beta template
    TH1F *Tof_template_antiproton_parametrized = new TH1F("Tof_Beta_p_parametrilized","",30,0.6,1.5);
    TH1F *Tof_data_pass7_negative_fit_parametrized = new TH1F("Tof_data_pass7_negative_fit_parametrized","",30,0.6,1.5);
    TH1F *Tof_data_pass7_positive_fit_parametrized = new TH1F("Tof_data_pass7_positive_fit_parametrized","",30,0.6,1.5);

    //string Richbetacut[16] = {"0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95"};
    //for (float Richbetacut = 0.990; Richbetacut<1.004; Richbetacut = Richbetacut + 0.001)

    // Parametrisation
    // 1. Tof fit: get antiproton tof_beta from antiproton MC, remove pion.
    Tof_template_antiproton->Fit("gaus","","",1.0,1.3); // to be tuned. 
    TF1 *Tof_template_antiprotonfit = Tof_template_antiproton->GetFunction("gaus");
    Utilities::ConvertToHistogram(Tof_template_antiprotonfit,*Tof_template_antiproton_parametrized);
    // 2. Tof fit: get backround template from negative data with gaussion template for pion and electrons
    Tof_data_pass7_negative->Fit("gaus","","",0.8,1.1); // to be tuned. 
    TF1 *Tof_data_pass7_negative_fit = Tof_data_pass7_negative->GetFunction("gaus");
    Utilities::ConvertToHistogram(Tof_data_pass7_negative_fit, *Tof_data_pass7_negative_fit_parametrized); 
    // 3. Tof fit: get antiproton tof_beta from positive data.
    Tof_data_pass7_positive->Fit("gaus","","",0.9,1.3); // to be tuned.
    TF1 *Tof_data_pass7_positive_fit = Tof_data_pass7_positive->GetFunction("gaus");
    Utilities::ConvertToHistogram(Tof_data_pass7_positive_fit, *Tof_data_pass7_positive_fit_parametrized);

    // Template Fit
    MYUtilities::TemplateFitter templateFitter_tof(0);
    // Signal
    //templateFitter_tof.AddTemplateHistogram(Tof_template_antiproton_parametrized);
    //templateFitter_tof.AddTemplateHistogram(Tof_template_antiproton);
    templateFitter_tof.AddTemplateHistogram(Tof_data_pass7_positive);
    // Background
    //templateFitter_tof.AddTemplateHistogram(Tof_data_pass7_negative_fit_parametrized);
    templateFitter_tof.AddTemplateHistogram(Tof_data_pass7_negative_template);
    //templateFitter_tof.AddTemplateHistogram(Tof_template_electron);
    //templateFitter_tof.AddTemplateHistogram(Tof_template_pion);
    // Data
    templateFitter_tof.SetDataHistogram(Tof_data_pass7_negative);
    // StartValue
    templateFitter_tof.SetStartValue(0, 500);

    // Fit
    templateFitter_tof.Fit(1,0.8,1.4); //If Fit range is not given, then defalt is whole range.

    // Template plots
    TCanvas * c8 = templateFitter_tof.CreateResultDrawing("Fit_Result_tof",1000,500);
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

    // Template Fit result
    vector<double> result_tof, ResultError_tof;
    assert(result_tof.empty());
    assert(ResultError_tof.empty());
    result_tof.assign(templateFitter_tof.GetAbsoluteResult().begin(), templateFitter_tof.GetAbsoluteResult().end());
    ResultError_tof.assign(templateFitter_tof.GetAbsoluteResultError().begin(), templateFitter_tof.GetAbsoluteResultError().end());
    double Chi2_tof = templateFitter_tof.Chi2();
    int NDF_tof = templateFitter_tof.NDF();
    double CHI2dof_tof = Chi2_tof/NDF_tof;

    ///////// Choose ratio you use according to different cuts  //////////////// 
    //v_ratio_result_tof.push_back(result_tof[0]/result_positive[0]*100000); proton from template fit
    //v_ratio_result_tof.push_back(result_tof[0]/data_pass7_positive->GetEntries()*100000); // tof fit, proton from total number after cut
    v_ratio_result_tof.push_back(result_tof[0]/tpass7_positive->GetEntries()*100000); // proton from total number before cut

    v_ratio_result_tof_with_effective.push_back(v_ratio_result_tof.back()*Effective_Acceptance->GetY()[index]);
    v_chi2_tof.push_back(CHI2dof_tof);
    v_error.push_back(ResultError_tof[0]/tpass7_positive->GetEntries()*100000);
    v_error_relative.push_back( (ResultError_tof[0]/tpass7_positive->GetEntries()) / (result_tof.at(0)/tpass7_positive->GetEntries()*Effective_Acceptance->GetY()[index]) *100);
     
    //cout << "raw negative events number:" << tpass7_negative->GetEntries() <<endl;
    //cout << "after cut negative events number:" << Tof_data_pass7_negative->GetEntries() <<endl;
    cout << "tof fit result:" << result_tof << endl;
    //cout << "proton total number from template fit:" << result_positive[0] << endl; // tof fit, proton from template fit   
    cout << "proton total number before cut:" << tpass7_positive->GetEntries()  << endl; // tof fit, proton from total number before cut
    //cout << "proton total number after cut:" << data_pass7_positive->GetEntries()  << endl; // tof fit, proton from total number after cut
    cout << "Chi2/dof_tof:" << CHI2dof_tof << endl;
    cout << "ratio from 1./tofbeta fit: "<< result_tof[0]/tpass7_positive->GetEntries()*100000  <<endl;

    // Save templates histograms in root file
    TFile testroot((string("Time_Averaged_ratio_Low/tof/averaged_ratio_") + doubleToString(binning[index]) + string("_") + doubleToString(binning[index+1]) + string("_") +issversion + string(".root")).c_str(),"RECREATE");
    Tof_template_antiproton->Write("Tof_template_antiproton");
    Tof_template_antiproton_parametrized->Write("Tof_template_antiproton_parametrized");
    Tof_template_electron->Write("Tof_template_electron");
    Tof_template_pion->Write("Tof_template_pion");
    Tof_data_pass7_negative->Write("Tof_data_pass7_negative");
    Tof_data_pass7_negative_fit_parametrized->Write("Tof_data_pass7_negative_fit_parametrized");
    Tof_data_pass7_positive->Write("Tof_data_pass7_positive");
    Tof_data_pass7_positive_fit_parametrized->Write("Tof_data_pass7_positive_fit_parametrized");
    testroot.Close();

    ////////////////////////////////////////////////////////////////
    ///////////////       Results  Plots        ////////////////////
    ////////////////////////////////////////////////////////////////
    TGraph *g_ratio_tof = new TGraph(v_ratio_result_tof.size(), subbincenter.data(), v_ratio_result_tof.data());
    TGraphErrors *g_ratio_tof_with_effective = new TGraphErrors(v_ratio_result_tof.size(), subbincenter.data(), v_ratio_result_tof_with_effective.data(), 0, v_error.data());
    TGraph *g_chi2_tof = new TGraph(v_chi2_tof.size(), subbincenter.data(), v_chi2_tof.data());
    TGraph *g_error_relative = new TGraph(v_error_relative.size(), subbincenter.data(), v_error_relative.data());
    TH1D h_error_relative = TH1D("", "", v_error_relative.size(), subbinedge.data());
    Utilities::ConvertToHistogram(g_error_relative,h_error_relative);

    TCanvas cc3("cc3","cc3",1000,500);
    g_ratio_tof->SetMarkerStyle(15);
    g_ratio_tof->SetMarkerColor(2);
    g_ratio_published->Draw("AP *");
    g_ratio_tof->Draw("same P");
    g_ratio_published->SetTitle("");
    g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
    g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");
    cc3.Update();
    cc3.SaveAs( (string("Time_Averaged_ratio_Low/plots/low_ratio_tof_") + issversion + string(".pdf")).c_str() );


    TCanvas cc33("cc33","cc33",1000,500);
    g_ratio_tof_with_effective->SetMarkerStyle(15);
    g_ratio_tof_with_effective->SetMarkerColor(2);
    g_IntermediateResult->SetMarkerStyle(15);
    g_IntermediateResult->SetMarkerColor(4);
    g_ratio_published->Draw("AP *");
    g_ratio_tof_with_effective->Draw("same P");
    g_IntermediateResult->Draw("same P");
    g_ratio_published->SetTitle("");
    g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
    g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");
    cc33.Update();
    cc33.SaveAs( (string("Time_Averaged_ratio_Low/plots/low_ratio_tof_with_effective_") + issversion + string(".pdf")).c_str() );


    TCanvas cc4("cc4","cc4",1000,500);
    g_chi2_tof->Draw("AP *");
    g_chi2_tof->SetTitle("");
    g_chi2_tof->GetXaxis()->SetTitle("Rigidity (GV)");
    g_chi2_tof->GetYaxis()->SetTitle("Chi2/dof");
    cc4.Update();
    cc4.SaveAs( (string("Time_Averaged_ratio_Low/plots/chi2_tof_") + issversion + string(".pdf")).c_str() );

    TCanvas ccerror("ccerror","ccerror",1000,500);
    h_error_relative.Draw("");
    h_error_relative.SetTitle("");
    h_error_relative.GetYaxis()->SetRangeUser(0,12);
    h_error_relative.GetXaxis()->SetTitle("Rigidity (GV)");
    h_error_relative.GetYaxis()->SetTitle("Relative Error (%)");
    h_error_relative.GetXaxis()->SetLabelSize(0.04);
    h_error_relative.GetXaxis()->SetTitleSize(0.04);
    h_error_relative.GetYaxis()->SetLabelSize(0.04);
    h_error_relative.GetYaxis()->SetTitleSize(0.04);
    h_error_relative.SetFillColor(0);
    h_error_relative.SetLineColor(1);
    h_Published_Statistic_Error_Relative.Draw("same");
    h_Published_Statistic_Error_Relative.SetLineColor(2);
    h_Published_Systematic_Error_Relative.Draw("same");
    h_Published_Systematic_Error_Relative.SetLineColor(4);
    TLegend * errorleg = new TLegend(0.4,0.75,0.7,0.85); //(xmin, ymin, xmax, ymax)
    errorleg->SetFillColor(0);
    if (issversion == "pass7.8"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (Full Range)","lp");}
    else if (issversion == "2016paper"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (Published Range)","lp");}
    errorleg->AddEntry(&h_Published_Statistic_Error_Relative,"Published_Statistic_Error_Relative","lp");
    errorleg->AddEntry(&h_Published_Systematic_Error_Relative,"Published_Systematic_Error_Relative","lp");
    gStyle->SetLegendTextSize(0.03);
    errorleg->Draw();
    ccerror.Update();
    ccerror.SaveAs((string("Time_Averaged_ratio_Low/plots/error_compare_") + issversion + string(".pdf")).c_str());

    // Save fit result in root file
    TFile ratioroot((string("Time_Averaged_ratio_Low/plots/Ratio")+ string("_") +issversion + string(".root")).c_str(),"RECREATE");
    Effective_Acceptance->Write("Effective_Acceptance");
    g_ratio_tof->Write("ratio_tof");
    g_ratio_tof_with_effective->Write("ratio_tof_with_effective");
    g_chi2_tof->Write("chi2_tof");
    g_ratio_published->Write("published");
    ratioroot.Close();

}  // rigidity loop end.

// Plot Effective_Acceptance, this Effective_Acceptance plot can not be put in rigidity loop, because it remove some points in Effective_Acceptance. Therefore has to be in last part.
TCanvas cc5("cc5","cc5",1000,500);
Effective_Acceptance->RemovePoint(0); //    MC1042 is from 1 to 800, so delete first point 0.8-1.0.
Effective_Acceptance->Draw("AP *");
Effective_Acceptance->SetTitle("");
Effective_Acceptance->GetXaxis()->SetTitle("Rigidity (GV)");
Effective_Acceptance->GetYaxis()->SetTitle("Effective_Acceptance");
Effective_Acceptance->GetXaxis()->SetLimits(1.0,800);
gPad->SetLogx();
cc5.Update();
cc5.SaveAs("Time_Averaged_ratio_Low/plots/Effective_Acceptance.pdf");

  return EXIT_SUCCESS;
}
