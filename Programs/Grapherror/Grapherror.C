// ACsoft includes
#include "AnalysisEvent.hh"
#include "AcceptanceManager.hh"
#include "AntiprotonBinning.hh"
#include "BinningDefinition.hh"
#include "BinningFunctions.hh"
#include "BinningTools.hh"
#include "ConfigHandler.hh"
#include "CutFactory.hh"
#include "CutAttachment.hh"
#include "Environment.hh"
#include "EfficiencyHistograms.hh"
#include "EventFactory.hh"
#include "FileManagerController.hh"
#include "FileManager.hh"
#include "GlobalOptions.hh"
#include "MPIEnvironment.hh"
#include "McSpectrumScaler.hh"
#include "MeasuringTime.hh"
#include "ObjectManager.hh"
#include "PredefinedBinnings.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "Utilities.hh"
#include "ValueHistograms.hh"

// ROOT includes
#include <math.h>
#include <TFile.h>
#include <TH1.h>
#include "TreeFormula.hh"
#include "Utilities.hh"
#include <TROOT.h>
#include "TreeWriter.hh"
#include <TApplication.h>
#include <TProof.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TAttMarker.h>
#include <TAttLine.h>

#define INFO_OUT_TAG "FluxCalculation"
#include "debugging.hh"

//// Produce flux, right now only on 16.6-259GV, last two points binnings are different, need to be correctd. 


int main(int argc, char** argv) {
/*
TFile file("TEfficiency_effective_acceptance.root");
TGraphAsymmErrors *old = (TGraphAsymmErrors*)file.Get("eff_graph");
TGraphAsymmErrors * newgraph = new TGraphAsymmErrors();


Double_t AntiprotonBinning_double_largerange[47] = {4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09,
    7.76, 8.48, 9.26, 10.1, 11, 12, 13, 14.1, 15.3,
    16.6, 18, 19.5, 21.1, 22.8, 24.7, 26.7, 28.8, 31.1, 33.5, 36.1,
    38.9, 41.9, 45.1, 48.5, 52.2, 56.1, 60.3, 64.8, 69.7, 74.9,
    80.5, 93, 108, 125, 147, 175, 211, 259, 330, 525};
AntiprotonBinning_length = 46;

TH1D  h_old = TH1D("","",AntiprotonBinning_length, AntiprotonBinning_double_largerange);
Utilities::ConvertToHistogram ( old, h_old);

//h_old.Smooth();
//h_old.Smooth();

Double_t ax[46],ay[46];
for (int j=0; j<46; j++){
ax[j] = h_old.GetXaxis()->GetBinCenter(j+1);
ay[j] = h_old.GetBinContent(j+1);
//old->GetPoint(j,ax[j],ay[j]);
std::cout<< ax[j] << std::endl;
std::cout<< ay[j] << std::endl;
//std::cout<< j  << std::endl;
}

for (int i=0; i<46; i++){
newgraph->SetPoint(i, ax[i], 1./ay[i] );
//newgraph->SetPointError(i, old->GetErrorXlow(i), old->GetErrorXhigh(i), old->GetErrorYlow(i), old->GetErrorYhigh(i));
newgraph->SetPointError(i, 0.0, 0.0, old->GetErrorYlow(i), old->GetErrorYhigh(i));
std::cout << i << std::endl;
}

TFile file2("TEfficiency_effective_acceptance2.root","RECREATE");
old->Write("old");
newgraph->Write("final");
file2.Close();
*/
  return EXIT_SUCCESS;
}
