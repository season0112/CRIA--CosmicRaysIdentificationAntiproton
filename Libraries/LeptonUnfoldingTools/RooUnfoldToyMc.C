#include "RooUnfoldToyMc.hh"

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

#define INFO_OUT_TAG "RooUnfoldToyMc"
#include <debugging.hh>

RooUnfoldToyMc::RooUnfoldToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, int randomSeed) :
  UnfoldingToyMc(energyBinning, param, randomSeed) {

}


RooUnfoldToyMc::~RooUnfoldToyMc() {
}

RooUnfoldResponse* RooUnfoldToyMc::GenerateToyResponseNoCutoff(unsigned long long nToyResponse) {

  const double w = 1.0;

  TH1D* hToyTruthShape = Make<TH1D>("hToyTruthShape", ";E_{truth}", fEnergyBinning);
  TH1D* hToyMeasShape = Make<TH1D>("hToyMeasShape", ";E_{meas}", fEnergyBinning);
  RooUnfoldResponse* toyResponse = new RooUnfoldResponse(hToyMeasShape, hToyTruthShape);

  delete hMigrationMatrixBeforeCutoff;
  hMigrationMatrixBeforeCutoff = Make<TH2D>("hMigrationMatrixBeforeCutoff", ";E_{meas};E_{truth}", fEnergyBinning, fEnergyBinning);

  double xlow = fEnergyBinning.Min();
  double xup = fEnergyBinning.Max();

  Utilities::ProgressBar pb(nToyResponse);
  if (fVerbosity > 0) {
    INFO_OUT << "Fill response matrix..." << std::endl;
    pb.PrintScale();
  }
  for (unsigned long long int i = 0; i < nToyResponse; ++i){
    if (fVerbosity > 0)
      pb.Update(i);

    double E = Utilities::PowerLawRandomNumber(xlow, xup, -1.0, fRandom);
    double Emeas = fRandomGenerator->GetRandomSmearedEnergy(E);

    hMigrationMatrixBeforeCutoff->Fill(Emeas, E);
    toyResponse->Fill(Emeas, E, w);
  }

  delete hToyMeasShape;
  delete hToyTruthShape;

  return toyResponse;
}

RooUnfoldResponse* RooUnfoldToyMc::GenerateToyResponseWithCutoff(unsigned long long nToyResponse, const TH1* hMeasuringTime) {

  const double w = 1.0;

  TH1D* hToyTruthShape = Make<TH1D>("hToyTruthShape", ";E_{truth}", fEnergyBinning);
  TH1D* hToyMeasShape = Make<TH1D>("hToyMeasShape", ";E_{meas}", fEnergyBinning);
  //TH2D* hToyMigrShape = Make<TH2D>("hToyMigrShape", ";E_{meas};E_{truth}", fEnergyBinning, fEnergyBinning); // TEST
  RooUnfoldResponse* toyResponse = new RooUnfoldResponse(hToyMeasShape, hToyTruthShape);

  delete hMigrationMatrixBeforeCutoff;
  hMigrationMatrixBeforeCutoff = Make<TH2D>("hMigrationMatrixBeforeCutoff", ";E_{meas};E_{truth}", fEnergyBinning, fEnergyBinning);

  double xlow = fEnergyBinning.Min();
  double xup = fEnergyBinning.Max();

  double meastime_max = hMeasuringTime->GetBinContent(hMeasuringTime->GetNbinsX());

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

    hMigrationMatrixBeforeCutoff->Fill(Emeas, E);

    double cutoffProb = hMeasuringTime->GetBinContent(hMeasuringTime->GetXaxis()->FindFixBin(Emeas)) / meastime_max;

    if (fRandom->Rndm() < cutoffProb)
      toyResponse->Fill(Emeas, E, w);
    else
      toyResponse->Miss(E, w);

    // TEST
//    if (fRandom->Rndm() < cutoffProb) {
//      hToyTruthShape->Fill(E, w);
//      hToyMeasShape->Fill(Emeas, w);
//      hToyMigrShape->Fill(Emeas, E, w);
//    }
//    else
//      hToyTruthShape->Fill(E, w);
  }

  // TEST
//  RooUnfoldResponse* toyResponse = new RooUnfoldResponse(hToyMeasShape, hToyTruthShape, hToyMigrShape);

  delete hToyMeasShape;
  delete hToyTruthShape;
//  delete hToyMigrShape; // TEST

  return toyResponse;
}

void RooUnfoldToyMc::SimulateUnfoldedMeasurement(const TH1* hAssumedFlux, const TH1* hMeasuringTime, const TH1* hTriggerEfficiency, const TH1* hAcceptance,
                                                 RooUnfoldResponse* toyResponse,
                                                 double fractionOfMeasuringTime, bool simulateCutoff) {

  // clean-up from previous run (if any)
  DeleteHistograms();

  const double w = 1.0;

  double meastime_max = hMeasuringTime->GetBinContent(hMeasuringTime->GetNbinsX());

  hToyAverageFlux = (TH1*)hAssumedFlux->Clone("hToyAverageFlux");
  hToyAverageEventCounts = (TH1*)hAssumedFlux->Clone("hToyAverageEventCounts");
  hToyAverageEventCounts->SetTitle("(before cutoff)");
  hToyAverageEventCounts->GetXaxis()->SetTitle("E_{truth} (GeV)");
  hToyAverageEventCounts->GetYaxis()->SetTitle("toy MC: true mean counts");

  hToyAverageEventCounts->Multiply(hAcceptance);
  hToyAverageEventCounts->Multiply(hTriggerEfficiency);
  if (simulateCutoff)
    hToyAverageEventCounts->Scale(fractionOfMeasuringTime * meastime_max);
  else {
    hToyAverageEventCounts->Multiply(hMeasuringTime);
    hToyAverageEventCounts->Scale(fractionOfMeasuringTime);
  }
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

      if (simulateCutoff){
        double cutoffProb = hMeasuringTime->GetBinContent(hMeasuringTime->GetXaxis()->FindFixBin(Emeas)) / meastime_max;
        if (fRandom->Rndm() < cutoffProb)
          hToyMeasEventCounts->Fill(Emeas, w);
      }
      else
        hToyMeasEventCounts->Fill(Emeas, w);
    }
  }

  if (fVerbosity > 0)
    INFO_OUT << "Toy MC: unfolding..." << std::endl;
  RooUnfoldBayes unfold (toyResponse, hToyMeasEventCounts, fBayesIterations, false);
  unfold.SetVerbose(fVerbosity);
  hToyUnfoldedEventCounts = (TH1D*)unfold.Hreco()->Clone("hToyUnfoldedEventCounts");

  // calculate fluxes, for easier comparison
  hToyTruthFlux = (TH1D*)hToyTruthEventCounts->Clone("hToyTruthFlux");
  hToyMeasFlux = (TH1D*)hToyMeasEventCounts->Clone("hToyMeasFlux");
  hToyUnfoldedFlux = (TH1D*)hToyUnfoldedEventCounts->Clone("hToyUnfoldedFlux");

  hToyTruthFlux->Divide(hAcceptance);
  hToyTruthFlux->Divide(hTriggerEfficiency);
  if (simulateCutoff)
    hToyTruthFlux->Scale(1.0 / (fractionOfMeasuringTime * meastime_max));
  else {
    hToyTruthFlux->Divide(hMeasuringTime);
    hToyTruthFlux->Scale(1.0/fractionOfMeasuringTime);
  }
  hToyTruthFlux->Scale(1.0, "width");

  hToyMeasFlux->Divide(hAcceptance);
  hToyMeasFlux->Divide(hTriggerEfficiency);
  hToyMeasFlux->Divide(hMeasuringTime);
  hToyMeasFlux->Scale(1.0 / fractionOfMeasuringTime);
  hToyMeasFlux->Scale(1.0, "width");

  hToyUnfoldedFlux->Divide(hAcceptance);
  hToyUnfoldedFlux->Divide(hTriggerEfficiency);
  if (simulateCutoff)
    hToyUnfoldedFlux->Scale(1.0 / (fractionOfMeasuringTime * meastime_max));
  else {
    hToyUnfoldedFlux->Divide(hMeasuringTime);
    hToyUnfoldedFlux->Scale(1.0/fractionOfMeasuringTime);
  }
  hToyUnfoldedFlux->Scale(1.0, "width");
}

