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

void timestamp_split(){

TFile file_negative("B1130_pass7_7.7_all_Tree_negative.root");
TTree *negativetree = (TTree*)file_negative.Get("AntiprotonIntermediateEnergyTree");
gROOT->cd();
TTree* Tree_negative = negativetree->CopyTree("TimeStamp<1432598400");
Tree_negative->SaveAs("B1130_pass7_7.7_all_Tree_negative_May2015.root");



TFile file_positive("B1130_pass7_7.7_all_Tree_positive.root");
TTree *positivetree = (TTree*)file_positive.Get("AntiprotonIntermediateEnergyTree");
gROOT->cd();
TTree* Tree_positive = positivetree->CopyTree("TimeStamp<1432598400");
Tree_positive->SaveAs("B1130_pass7_7.7_all_Tree_positive_May2015.root");

}
