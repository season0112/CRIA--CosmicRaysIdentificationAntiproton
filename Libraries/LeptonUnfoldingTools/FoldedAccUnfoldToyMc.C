#include "FoldedAccUnfoldToyMc.hh"

#include "MigrationParameterization.hh"
#include "MigrationRandomGenerator.hh"
#include "FoldedAcceptanceUnfoldingForLeptons.hh"

#include <BinningTools.hh>
#include <ProgressBar.hh>
#include <Statistics.hh>
#include <Utilities.hh>

#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#define INFO_OUT_TAG "FoldedAccUnfoldToyMc"
#include <debugging.hh>

FoldedAccUnfoldToyMc::FoldedAccUnfoldToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, int randomSeed) :
  UnfoldingToyMc(energyBinning, param, randomSeed) {

}


FoldedAccUnfoldToyMc::~FoldedAccUnfoldToyMc() {
}

TH2* FoldedAccUnfoldToyMc::GenerateToyResponse(unsigned long long nToyResponse) {

  const double w = 1.0;

  double xlow = fEnergyBinning.Min();
  double xup = fEnergyBinning.Max();

  TH2D* toyResponse = Make<TH2D>("hToyResponseFoldedAcc", ";E_{meas};E_{truth}", fEnergyBinning, fEnergyBinning);

  Utilities::ProgressBar pb(nToyResponse);
  if (fVerbosity > 0) {
    INFO_OUT << "Fill response matrix..." << std::endl;
    pb.PrintScale();
  }
  for (unsigned long long int i = 0; i < nToyResponse; ++i){
    if (fVerbosity > 0)
      pb.Update(i);

    double E = Utilities::PowerLawRandomNumber(xlow, xup, -1.0, gRandom);
    double Emeas = fRandomGenerator->GetRandomSmearedEnergy(E);

    toyResponse->Fill(Emeas, E, w);
  }

  return toyResponse;
}

void FoldedAccUnfoldToyMc::SimulateUnfoldedMeasurement(const TH1* hAssumedFlux, const TH1* hMeasuringTime, const TH1* hTriggerEfficiency, const TH1* hAcceptance,
                                                       const TH2* toyResponse,
                                                       double fractionOfMeasuringTime) {

  // clean-up from previous run (if any)
  DeleteHistograms();

  const double w = 1.0;

  TH1* hMeasuringTimeForToymc = (TH1*)hMeasuringTime->Clone("hMeasuringTimeForToymc");
  hMeasuringTimeForToymc->Scale(fractionOfMeasuringTime);

  double meastime_max = hMeasuringTimeForToymc->GetBinContent(hMeasuringTimeForToymc->GetNbinsX());

  hToyAverageFlux = (TH1*)hAssumedFlux->Clone("hToyAverageFlux");
  hToyAverageEventCounts = (TH1*)hAssumedFlux->Clone("hToyAverageEventCounts");
  hToyAverageEventCounts->SetTitle("(before cutoff)");
  hToyAverageEventCounts->GetXaxis()->SetTitle("E_{truth} (GeV)");
  hToyAverageEventCounts->GetYaxis()->SetTitle("toy MC: true mean counts");

  hToyAverageEventCounts->Multiply(hAcceptance);
  hToyAverageEventCounts->Multiply(hTriggerEfficiency);
  hToyAverageEventCounts->Scale(meastime_max);
  Utilities::ScaleByBinWidth(*hToyAverageEventCounts);

  // set poisson errors
  for (int ibin = 1; ibin <= hToyAverageEventCounts->GetNbinsX(); ++ibin)
    hToyAverageEventCounts->SetBinError(ibin, TMath::Sqrt(hToyAverageEventCounts->GetBinContent(ibin)));

  // smear average counts to arrive at counts as a function of true energy
  // and set poisson errors
  hToyTruthEventCounts = Make<TH1D>("hToyTruthEventCounts", "(before cutoff);E_{truth} (GeV);toy MC: measured event counts", fEnergyBinning);
  for (int ibin = 1; ibin <= hToyTruthEventCounts->GetNbinsX(); ++ibin) {
    hToyTruthEventCounts->SetBinContent(ibin, fRandom->Poisson(hToyAverageEventCounts->GetBinContent(ibin)));
    hToyTruthEventCounts->SetBinError(ibin, TMath::Sqrt(hToyTruthEventCounts->GetBinContent(ibin)));
  }

  unsigned long long int nToySim = (unsigned long long int)(hToyTruthEventCounts->Integral());
  unsigned long long int iTotal = 0;

  Utilities::ProgressBar pb2(nToySim);
  if (fVerbosity > 0) {
    INFO_OUT << "Toy MC: simulate measurement, using " << nToySim << " events..." << std::endl;
    pb2.PrintScale();
  }

  hToyMeasEventCounts = Make<TH1D>("hToyMeasEventCounts", ";E_{meas} (GeV);toy MC: measured event counts", fEnergyBinning);

  for (int ebin = 1; ebin <= hToyAverageEventCounts->GetNbinsX(); ++ebin){

    int nEventsInBin = TMath::Nint(hToyAverageEventCounts->GetBinContent(ebin));
    double E = hToyAverageEventCounts->GetXaxis()->GetBinCenterLog(ebin);

    for (int iEvent = 0; iEvent < nEventsInBin; ++iEvent) {

      if (fVerbosity > 0)
        pb2.Update(iTotal++);

      double Emeas = fRandomGenerator->GetRandomSmearedEnergy(E);

      double cutoffProb = hMeasuringTimeForToymc->GetBinContent(hMeasuringTimeForToymc->GetXaxis()->FindFixBin(Emeas)) / meastime_max;
      if (fRandom->Rndm() < cutoffProb)
        hToyMeasEventCounts->Fill(Emeas, w);
    }
  }

  if (fVerbosity > 0)
    INFO_OUT << "Toy MC: unfolding..." << std::endl;

  FoldedAcceptanceUnfoldingForLeptons unfolding(hToyMeasEventCounts, hMeasuringTimeForToymc, hAcceptance, hTriggerEfficiency, toyResponse, fVerbosity);
  hToyUnfoldedFlux = unfolding.MakeUnfoldedFlux();
  hToyFoldedUnfoldedFlux = unfolding.MakeFoldedUnfoldedFlux();

  hToyUnfoldedEventCounts = nullptr;
  if (hToyUnfoldedFlux) {
    hToyUnfoldedEventCounts = (TH1*)hToyUnfoldedFlux->Clone("hToyUnfoldedEventCounts");
    hToyUnfoldedEventCounts->Multiply(hAcceptance);
    hToyUnfoldedEventCounts->Multiply(hTriggerEfficiency);
    hToyUnfoldedEventCounts->Scale(meastime_max);
    Utilities::ScaleByBinWidth(*hToyUnfoldedEventCounts);
  }

  // calculate fluxes, for easier comparison
  hToyTruthFlux = (TH1D*)hToyTruthEventCounts->Clone("hToyTruthFlux");
  hToyTruthFlux->Divide(hAcceptance);
  hToyTruthFlux->Divide(hTriggerEfficiency);
  hToyTruthFlux->Scale(1.0 / meastime_max);
  hToyTruthFlux->Scale(1.0, "width");

  hToyMeasFlux = (TH1D*)hToyMeasEventCounts->Clone("hToyMeasFlux");
  hToyMeasFlux->Divide(hAcceptance);
  hToyMeasFlux->Divide(hTriggerEfficiency);
  hToyMeasFlux->Divide(hMeasuringTimeForToymc);
  hToyMeasFlux->Scale(1.0, "width");

  delete hMeasuringTimeForToymc;
}

TGraph* FoldedAccUnfoldToyMc::MakeFoldedOverMeasuredRatio() const {

  if (!hToyFoldedUnfoldedFlux)
    return nullptr;

  TGraph* grFoldedMeasuredRatio = new TGraph;
  grFoldedMeasuredRatio->SetName("grFoldedMeasuredRatio");
  for (int ibin = 1; ibin <= hToyMeasFlux->GetNbinsX(); ++ibin) {

    double E = hToyMeasFlux->GetXaxis()->GetBinCenterLog(ibin);
    double measFlux = hToyMeasFlux->GetBinContent(ibin);

    double toyFoldedMeasRatio = hToyFoldedUnfoldedFlux->GetBinContent(ibin) / measFlux - 1.0;
    grFoldedMeasuredRatio->SetPoint(ibin-1, E, toyFoldedMeasRatio);
  }

  return grFoldedMeasuredRatio;
}


