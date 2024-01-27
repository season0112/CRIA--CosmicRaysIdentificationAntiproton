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
#include <sstream>
#include <string>
#include <unistd.h>

#include "PredefinedBinnings.hh"
#include "TimeUtilities.hh"
#include "AMSUnixTimes.hh"

#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>

#define INFO_OUT_TAG "BartelRotation"
#include "debugging.hh"

#define BeginOfIssTime 1305800000.0
#define EndOfPass7Time 1527491432.0
#define EndOfPass7extTime 1546405413.0


int main(int argc, char* argv[]) {
(void) argc;
(void) argv;
/*
std::string pass7negative_name = "pass7_negative_";
std::string pass7positive_name = "pass7_positive_";
std::string rootname = ".root";

const char* pass7_negative = "B1130_pass7_7.7_all_Tree_negative.root";
const char* pass7_positive = "B1130_pass7_7.7_all_Tree_positive.root";

TFile *fpass7_negative = new TFile(pass7_negative);
TFile *fpass7_positive = new TFile(pass7_positive);

TTree *tpass7_negative = (TTree*)fpass7_negative->Get("AntiprotonIntermediateEnergyTree");
TTree *tpass7_positive = (TTree*)fpass7_positive->Get("AntiprotonIntermediateEnergyTree");

tpass7_negative->SetBranchStatus("*",0);
tpass7_negative->SetBranchStatus("RichBeta",1);
tpass7_negative->SetBranchStatus("TrdLogLikelihoodRatioElectronProtonTracker",1);
tpass7_negative->SetBranchStatus("TimeStamp",1);

//// rootfiles generation
chdir("./rootfiles");
for (int i=2426;i<2520;i++){
std::string pass7negative2 = pass7negative_name + std::to_string(i) + rootname;
std::string pass7positive2 = pass7positive_name + std::to_string(i) + rootname;
const char* pass7negative = (char*)pass7negative2.c_str();
const char* pass7positive = (char*)pass7positive2.c_str();
TFile f1(pass7negative, "RECREATE");
TTree *newtree = tpass7_negative->CloneTree(0);
newtree->Write("AntiprotonIntermediateEnergyTree");
f1.Close(); 
TFile f2(pass7positive, "RECREATE");
TTree *newtree2 = tpass7_negative->CloneTree(0);
newtree2->Write("AntiprotonIntermediateEnergyTree");
f2.Close();
}
chdir("..");



////  Filling TimeStampt
TH1D *timestamp_pass7_negative = new TH1D("TimeStamp_pass7_negative","",100, BeginOfIssTime, EndOfPass7Time);
TH1D *timestamp_pass7_positive = new TH1D("TimeStamp_pass7_positive","",100, BeginOfIssTime, EndOfPass7Time);

UInt_t TimeStamp_pass7_negative, TimeStamp_pass7_positive, TimeStamp_pass7_negative_finebin, TimeStamp_pass7_positive_finebin;
float TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative,TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive,TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative_finebin,TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive_finebin;
float RichBeta_pass7_negative, RichBeta_pass7_positive, RichBeta_pass7_negative_finebin, RichBeta_pass7_positive_finebin;

int BartelsRotationNumber_pass7_negative;
int BartelsRotationNumber_pass7_positive;

tpass7_negative->SetBranchAddress("TimeStamp",&TimeStamp_pass7_negative);
tpass7_positive->SetBranchAddress("TimeStamp",&TimeStamp_pass7_positive);
tpass7_negative->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative);
tpass7_positive->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive);
tpass7_negative->SetBranchAddress("RichBeta",&RichBeta_pass7_negative);
tpass7_positive->SetBranchAddress("RichBeta",&RichBeta_pass7_positive);

int entries = tpass7_positive->GetEntries();
for( int i=0;i<entries;i++){
//for( int i=0;i<2000;i++){
     tpass7_positive->GetEntry(i);
     timestamp_pass7_positive->Fill(TimeStamp_pass7_positive);

     BartelsRotationNumber_pass7_positive = Utilities::Time::BartelsRotationNumber(TimeStamp_pass7_positive);
     std::string p7_p_bartel2 = pass7positive_name + std::to_string(BartelsRotationNumber_pass7_positive) + rootname;
     const char* p7_p_bartel = (char*)p7_p_bartel2.c_str();
     chdir("./rootfiles");
     TFile f(p7_p_bartel,"UPDATE");
     TTree *t3 = (TTree*)f.Get("AntiprotonIntermediateEnergyTree");
     t3->SetBasketSize("*",4);

     t3->SetBranchAddress("TimeStamp",&TimeStamp_pass7_positive_finebin);
     t3->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive_finebin);
     t3->SetBranchAddress("RichBeta",&RichBeta_pass7_positive_finebin);
*/
/*
     std::cout << t3->GetEntries() << std::endl;
     t3->SetEntries(t3->GetEntries());
     auto TimeStampBranch = t3->GetBranch("TimeStamp");
     TimeStampBranch->Fill();
     auto TrdLogLikelihoodRatioElectronProtonTrackerBranch = t3->GetBranch("TrdLogLikelihoodRatioElectronProtonTracker");
     TrdLogLikelihoodRatioElectronProtonTrackerBranch->Fill();
     auto RichBetaBranch = t3->GetBranch("RichBeta");
     RichBetaBranch->Fill();
*/
/*
     TimeStamp_pass7_positive_finebin = TimeStamp_pass7_positive;
     TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive_finebin = TrdLogLikelihoodRatioElectronProtonTracker_pass7_positive;
     RichBeta_pass7_positive_finebin = RichBeta_pass7_positive;

     t3->Fill();
     t3->Write("AntiprotonIntermediateEnergyTree", TObject::kOverwrite);
     f.Close();
     chdir("..");
     std::cout<< i << std::endl;
}
*/
/*
int entries3 = tpass7_negative->GetEntries();
//for( int i=0;i<entries3;i++){
for( int i=0;i<2000;i++){
     tpass7_negative->GetEntry(i);
     timestamp_pass7_negative->Fill(TimeStamp_pass7_negative);

     BartelsRotationNumber_pass7_negative = Utilities::Time::BartelsRotationNumber(TimeStamp_pass7_negative);
     std::string p7_n_bartel2 = pass7negative_name + std::to_string(BartelsRotationNumber_pass7_negative) + rootname;
     const char* p7_n_bartel = (char*)p7_n_bartel2.c_str();
     chdir("./rootfiles");
     TFile f2(p7_n_bartel,"UPDATE");
     TTree *t4 = (TTree*)f2.Get("AntiprotonIntermediateEnergyTree");

//     t4->SetBranchAddress("TimeStamp",&TimeStamp_pass7_negative_finebin);
     t4->SetBranchAddress("TrdLogLikelihoodRatioElectronProtonTracker",&TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative_finebin);
     t4->SetBranchAddress("RichBeta",&RichBeta_pass7_negative_finebin);

//     auto TimeStampBranch2 = t4->GetBranch("TimeStamp");
//     TimeStampBranch2->Fill();
//     auto TrdLogLikelihoodRatioElectronProtonTrackerBranch2 = t4->GetBranch("TrdLogLikelihoodRatioElectronProtonTracker");
//     TrdLogLikelihoodRatioElectronProtonTrackerBranch2->Fill();
//     auto RichBetaBranch2 = t4->GetBranch("RichBeta");
//     RichBetaBranch2->Fill();

//     TimeStamp_pass7_negative_finebin = TimeStamp_pass7_negative;
     TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative_finebin = TrdLogLikelihoodRatioElectronProtonTracker_pass7_negative;
     RichBeta_pass7_negative_finebin = RichBeta_pass7_negative;

     t4->Fill();
     t4->Write("AntiprotonIntermediateEnergyTree", TObject::kOverwrite);
     f2.Close();
     chdir("..");
}
*/

/*
TCanvas * c1 = new TCanvas;
timestamp_pass7_positive->Draw("COLZ");
c1->SaveAs("bartel_positive.png");

TCanvas * c9 = new TCanvas;
timestamp_pass7_negative->Draw("COLZ");
c9->SaveAs("bartel_negative.png");


*/
  return EXIT_SUCCESS;
}








