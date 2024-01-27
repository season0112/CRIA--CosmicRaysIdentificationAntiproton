#include "UnfoldingToyMc.hh"

#include "MigrationParameterization.hh"
#include "MigrationRandomGenerator.hh"

#include <BinningTools.hh>
#include <ProgressBar.hh>
#include <Statistics.hh>
#include <Utilities.hh>

#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#define INFO_OUT_TAG "UnfoldingToyMc"
#include <debugging.hh>

UnfoldingToyMc::UnfoldingToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, int randomSeed)
  : fEnergyBinning(energyBinning)
  , fParam(param)
{
  fRandom = new TRandom3(randomSeed);
  fRandomGenerator = new MigrationRandomGenerator(fEnergyBinning, fParam);
  fRandomGenerator->SetRandomGenerator(fRandom);
}


UnfoldingToyMc::~UnfoldingToyMc() {

  DeleteHistograms();

  delete fRandom;
  delete fRandomGenerator;
}


TGraphErrors* UnfoldingToyMc::MakeMeasuredOverTruthRatio() const {

  TGraphErrors* grToyMeasuredTruthRatio = new TGraphErrors;
  grToyMeasuredTruthRatio->SetName("grToyMeasuredTruthRatio");
  for (int ibin = 1; ibin <= hToyTruthFlux->GetNbinsX(); ++ibin) {

    double E = hToyTruthFlux->GetXaxis()->GetBinCenterLog(ibin);
    double truthFlux = hToyTruthFlux->GetBinContent(ibin);

    double toyMeasTruthRatio = hToyMeasFlux->GetBinContent(ibin) / truthFlux - 1.0;
    grToyMeasuredTruthRatio->SetPoint(ibin-1, E, toyMeasTruthRatio);
    grToyMeasuredTruthRatio->SetPointError(ibin-1, 0.0, hToyMeasFlux->GetBinError(ibin) / truthFlux);
  }

  return grToyMeasuredTruthRatio;
}

TGraphErrors* UnfoldingToyMc::MakeUnfoldedOverTruthRatio() const {

  if (!hToyUnfoldedFlux) {
    WARN_OUT << "Received nullptr for unfolded flux!" << std::endl;
    return nullptr;
  }

  TGraphErrors* grToyUnfoldedTruthRatio = new TGraphErrors;
  grToyUnfoldedTruthRatio->SetName("grToyUnfoldedTruthRatio");
  for (int ibin = 1; ibin <= hToyTruthFlux->GetNbinsX(); ++ibin) {

    double E = hToyTruthFlux->GetXaxis()->GetBinCenterLog(ibin);
    double truthFlux = hToyTruthFlux->GetBinContent(ibin);   
    double toyUnfoldedTruthRatio = hToyUnfoldedFlux->GetBinContent(ibin) / truthFlux - 1.0;
    grToyUnfoldedTruthRatio->SetPoint(ibin-1, E, toyUnfoldedTruthRatio);
    grToyUnfoldedTruthRatio->SetPointError(ibin-1, 0.0, hToyUnfoldedFlux->GetBinError(ibin) / truthFlux);
  }

  return grToyUnfoldedTruthRatio;
}

TGraphErrors* UnfoldingToyMc::MakeUnfoldedOverAverageRatio() const {

  if (!hToyUnfoldedFlux) {
    WARN_OUT << "Received nullptr for unfolded flux!" << std::endl;
    return nullptr;
  }

  TGraphErrors* grToyUnfoldedAverageRatio = new TGraphErrors;
  grToyUnfoldedAverageRatio->SetName("grToyUnfoldedAverageRatio");
  for (int ibin = 1; ibin <= hToyAverageFlux->GetNbinsX(); ++ibin) {

    double E = hToyAverageFlux->GetXaxis()->GetBinCenterLog(ibin);
    double averageFlux = hToyAverageFlux->GetBinContent(ibin);
    double toyUnfoldedAverageRatio = hToyUnfoldedFlux->GetBinContent(ibin) / averageFlux - 1.0;
    grToyUnfoldedAverageRatio->SetPoint(ibin-1, E, toyUnfoldedAverageRatio);
    grToyUnfoldedAverageRatio->SetPointError(ibin-1, 0.0, hToyUnfoldedFlux->GetBinError(ibin) / averageFlux);
  }

  return grToyUnfoldedAverageRatio;
}

TGraphErrors* UnfoldingToyMc::MakeTruthOverAverageRatio() const {

  if (!hToyTruthFlux) {
    WARN_OUT << "Received nullptr for truth flux!" << std::endl;
    return nullptr;
  }

  TGraphErrors* grToyTruthAverageRatio = new TGraphErrors;
  grToyTruthAverageRatio->SetName("grToyTruthAverageRatio");
  for (int ibin = 1; ibin <= hToyAverageFlux->GetNbinsX(); ++ibin) {

    double E = hToyAverageFlux->GetXaxis()->GetBinCenterLog(ibin);
    double averageFlux = hToyAverageFlux->GetBinContent(ibin);
    double toyTruthAverageRatio = hToyTruthFlux->GetBinContent(ibin) / averageFlux - 1.0;
    grToyTruthAverageRatio->SetPoint(ibin-1, E, toyTruthAverageRatio);
    grToyTruthAverageRatio->SetPointError(ibin-1, 0.0, hToyTruthFlux->GetBinError(ibin) / averageFlux);
  }

  return grToyTruthAverageRatio;

}

void UnfoldingToyMc::DeleteHistograms() {

  delete hToyAverageEventCounts;
  delete hToyTruthEventCounts;
  delete hToyMeasEventCounts;
  delete hToyUnfoldedEventCounts;
  delete hToyAverageFlux;
  delete hToyTruthFlux;
  delete hToyMeasFlux;
  delete hToyUnfoldedFlux;
}

