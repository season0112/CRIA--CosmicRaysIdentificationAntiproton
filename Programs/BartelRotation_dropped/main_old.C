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

std::string pass7negative_name = "pass7_negative_";
std::string pass7positive_name = "pass7_positive_";
std::string rootname = ".root";

//Float_t new_v;

//Int_t split = 0;
//Int_t bsize = 64000;
//Event  *event = 0;

//// rootfiles generation
chdir("./rootfiles");
//const char* pass7_negative = "B1130_pass7_7.7_all_Tree_negative.root";
//TFile *fpass7_negative = new TFile(pass7_negative);
//TTree *tpass7_negative = (TTree*)fpass7_negative->Get("AntiprotonIntermediateEnergyTree");
for (int i=2426;i<2520;i++){
std::string pass7negative2 = pass7negative_name + std::to_string(i) + rootname;
std::string pass7positive2 = pass7positive_name + std::to_string(i) + rootname;
const char* pass7negative = (char*)pass7negative2.c_str();
const char* pass7positive = (char*)pass7positive2.c_str();
TFile f1(pass7negative, "RECREATE");
//TTree::TTree(const TTree & tt);
//TTree tree = new TTree(*tpass7_negative);
TTree *tree = new TTree();
//TBranch *newBranch = tree->Branch("new_v", &new_v, "new_v/F");
//newBranch->Fill();
tree->Write("tree");
//TFile f2(pass7positive, "NEW");
TFile f2(pass7positive, "RECREATE");
TTree *tree2 = new TTree();
tree2->Write("tree");
}
chdir("..");



////  Filling TimeStampt
const char* pass7_negative = "B1130_pass7_7.7_all_Tree_negative.root";
const char* pass7_positive = "B1130_pass7_7.7_all_Tree_positive.root";

TFile *fpass7_negative = new TFile(pass7_negative);
TFile *fpass7_positive = new TFile(pass7_positive);

TTree *tpass7_negative = (TTree*)fpass7_negative->Get("AntiprotonIntermediateEnergyTree");
TTree *tpass7_positive = (TTree*)fpass7_positive->Get("AntiprotonIntermediateEnergyTree");

TH1D *timestamp_pass7_negative = new TH1D("TimeStamp_pass7_negative","",100, BeginOfIssTime, EndOfPass7Time);
TH1D *timestamp_pass7_positive = new TH1D("TimeStamp_pass7_positive","",100, BeginOfIssTime, EndOfPass7Time);

UInt_t TimeStamp_pass7_negative,TimeStamp_pass7_positive;
//int BartelsRotationNumber_pass7_negative;
int BartelsRotationNumber_pass7_positive;

tpass7_negative->SetBranchAddress("TimeStamp",&TimeStamp_pass7_negative);
tpass7_positive->SetBranchAddress("TimeStamp",&TimeStamp_pass7_positive);

//int entries = tpass7_positive->GetEntries();
//for( int i=0;i<entries;i++){
for( int i=0;i<100;i++){
     tpass7_positive->GetEntry(i);
     timestamp_pass7_positive->Fill(TimeStamp_pass7_positive);

     BartelsRotationNumber_pass7_positive = Utilities::Time::BartelsRotationNumber(TimeStamp_pass7_positive);
     std::string p7_p_bartel2 = pass7positive_name + std::to_string(BartelsRotationNumber_pass7_positive) + rootname;
     const char* p7_p_bartel = (char*)p7_p_bartel2.c_str();
     chdir("./rootfiles");
     TFile f(p7_p_bartel,"UPDATE");
//     auto t3 = f.Get<TTree>("tree");
     TTree *t3 = (TTree*)f.Get("tree");
     std::cout << t3->GetEntries() << std::endl; 

     if (t3->GetEntries()==0)  
     {
         auto newBranch = t3->Branch("TimeStamp", &TimeStamp_pass7_positive, "TimeStamp/i");
         newBranch->Fill();
         t3->Fill();
     std::cout << "yes" << std::endl; 
     std::cout << i << std::endl;
     std::cout << TimeStamp_pass7_positive << std::endl;
    }
     else
     {
         t3->SetBranchAddress("TimeStamp",&TimeStamp_pass7_negative);
         t3->Fill();
     std::cout << "no" << std::endl;
     std::cout << i << std::endl;
     std::cout << TimeStamp_pass7_positive << std::endl;
     }
     t3->Write("tree", TObject::kOverwrite);
//     tree->Branch("event", "Event", &event, bsize,split);
//     tpass7_negative->Fill();
//     tpass7_negative->Fill();     
     chdir("..");
}

int entries3 = tpass7_negative->GetEntries();
for( int i=0;i<entries3;i++){
     tpass7_negative->GetEntry(i);
     timestamp_pass7_negative->Fill(TimeStamp_pass7_negative);
//     BartelsRotationNumber_pass7_negative=Utilities::Time::BartelsRotationNumber(TimeStamp_pass7_negative);
}

TCanvas * c1 = new TCanvas;
timestamp_pass7_positive->Draw("COLZ");
c1->SaveAs("bartel_positive.png");

TCanvas * c9 = new TCanvas;
timestamp_pass7_negative->Draw("COLZ");
c9->SaveAs("bartel_negative.png");



/*
 * std::cout << Utilities::Time::BartelsRotationNumber(1527491432) << std::endl;
 * std::cout << Utilities::Time::StartTimeOfBartelsRotation(2521) << std::endl;
 * std::cout << Utilities::AMSUnixTimes::BeginOfIssTime() << std::endl;
*/
/*
     std::cout << TimeStamp_pass7_negative << std::endl;
     std::cout << (float)(BartelsRotationNumber_pass7_negative-2426)/3 << std::endl;
*/

  return EXIT_SUCCESS;
}








