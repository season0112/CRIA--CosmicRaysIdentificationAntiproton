#include <iostream>
#include <fstream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"

void split(){

TFile file_antiproton("B1042_antipr.pl1.1800_7.6_all_Tree_all.root");
TTree *AntiprotonIntermediateEnergyTree_antiproton = (TTree*)file_antiproton.Get("AntiprotonIntermediateEnergyTree");
gROOT->cd();
TTree* Tree_antiproton = AntiprotonIntermediateEnergyTree_antiproton->CopyTree("Rigidity>-5.9");
Tree_antiproton->SaveAs("B1042_antipr.pl1.1800_7.6_all_Tree.root");

/*
TFile file_electron("B1091_el.pl1.0_25200_7.6_all_Tree_all.root");
TTree *AntiprotonIntermediateEnergyTree_electron = (TTree*)file_electron.Get("AntiprotonIntermediateEnergyTree");
gROOT->cd();
TTree* Tree_electron = AntiprotonIntermediateEnergyTree_electron->CopyTree("Rigidity>-5.9");
Tree_electron->SaveAs("B1091_el.pl1.0_25200_7.6_all_Tree.root");

TFile file_negative("B1130_pass7_7.7_all_Tree_negative_all.root");
TTree *AntiprotonIntermediateEnergyTree_negative = (TTree*)file_negative.Get("AntiprotonIntermediateEnergyTree");
gROOT->cd();
TTree* Tree_negative = AntiprotonIntermediateEnergyTree_negative->CopyTree("Rigidity>-5.9");
Tree_negative->SaveAs("B1130_pass7_7.7_all_Tree_negative.root");


TFile file_positive("B1130_pass7_7.7_all_Tree_positive_all.root");
TTree *AntiprotonIntermediateEnergyTree_positive = (TTree*)file_positive.Get("AntiprotonIntermediateEnergyTree");
gROOT->cd();
TTree* Tree_positive = AntiprotonIntermediateEnergyTree_positive->CopyTree("Rigidity<5.9");
Tree_positive->SaveAs("B1130_pass7_7.7_all_Tree_positive.root");
*/

}



