#include "ExampleAnalysisTree.hh"

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
#include "TemplateFitter.hh"
#include <iostream>
#include <cassert>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include "TemplateFitter2D.hh"
#include <vector>
#include <TCanvas.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include "IterativeUnfolding.hh"
#include "AcceptanceUnfolding.hh"


using namespace std;

#define INFO_OUT_TAG "Unfolding_Iteration"
#include "debugging.hh"

int main(int argc, char* argv[]) {
//int main(char* argv[]) {
(void)argc;

//char* s2;
//s2 = argv[1];

TFile *f1 = new TFile("ratio_TH2D.root");
TFile *f2 = new TFile("Unfolding_MatricesTH2D.root");

TH2D *hMigrationMatrix = (TH2D*)f2->Get("Unfolding_Matrices");
TH1D *hUnfoldInputRate = (TH1D*)f1->Get(argv[1]);

Unfolding::IterativeUnfolding iterativeUnfolding(hMigrationMatrix, hUnfoldInputRate, 32, 1.0, 1e-3);
iterativeUnfolding.Run();
const TH1* hUnfoldedRate = iterativeUnfolding.UnfoldedRate();

iterativeUnfolding.Plot();

for (int i=1;i<20;i++){
cout << i << endl;
cout << hUnfoldedRate->GetBinContent(i) << endl;
}
  return EXIT_SUCCESS;
}







