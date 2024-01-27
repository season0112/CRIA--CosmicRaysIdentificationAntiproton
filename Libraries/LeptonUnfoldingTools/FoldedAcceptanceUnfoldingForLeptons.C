#include "FoldedAcceptanceUnfoldingForLeptons.hh"

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>

#include <Utilities.hh>
#include <AcceptanceUnfolding.hh>
#include <CubicSplineFit.hh>

#define INFO_OUT_TAG "FoldedAcceptanceUnfoldingForLeptons"
#include <debugging.hh>

FoldedAcceptanceUnfoldingForLeptons::FoldedAcceptanceUnfoldingForLeptons(const TH1* hEventCounts, const TH1* hMeasuringTime, const TH1* hAcceptance, const TH1* hTriggerEfficiency,
                                                                         const TH2* hMigrationMatrix, int verbosity)
  : fMigrationMatrix(hMigrationMatrix)
  , fAcceptance(hAcceptance)
  , fEventCounts(hEventCounts)
  , fMeasuringTime(hMeasuringTime)
  , fTriggerEfficiency(hTriggerEfficiency)
{

  fRateToBeUnfolded = (TH1*)hEventCounts->Clone("hRateToBeUnfolded");
  fRateToBeUnfolded->SetTitle("");
  fRateToBeUnfolded->GetYaxis()->SetTitle("unfolded flux");
  fRateToBeUnfolded->Divide(hMeasuringTime);

  fAcceptanceUnfolding = new Unfolding::AcceptanceUnfolding(hMigrationMatrix, fRateToBeUnfolded,
                                                            Unfolding::AcceptanceUnfolding::RegularizationHandling::ApplyRegularization,
                                                            24, 4e-3, true);

  fAcceptanceForUnfolding = (TH1*)hAcceptance->Clone("hAcceptanceForUnfolding");
  fAcceptanceForUnfolding->Multiply(hTriggerEfficiency);
  fAcceptanceForUnfolding->SetMinimum();
  fAcceptanceForUnfolding->SetMaximum();
  fAcceptanceForUnfolding->SetTitle("");
  fAcceptanceUnfolding->SetEffectiveAcceptance(fAcceptanceForUnfolding);
  fAcceptanceUnfolding->SetRegularizationSplineMinPositionX(0.5);
  fAcceptanceUnfolding->SetRegularizationSplineMaxPositionX(1000.0);
  fAcceptanceUnfolding->SetNodes({0.6, 1.25, 2.5, 5.0, 10., 30., 100., 200., 700.});
  fAcceptanceUnfolding->SetConvergenceCheckMin(1.0);
  fAcceptanceUnfolding->SetConvergenceCheckMax(700.0);
  fAcceptanceUnfolding->SetSplineXValueHandling(Math::CubicSplineFit::ValueHandling::Logarithmic);
  fAcceptanceUnfolding->SetSplineYValueHandling(Math::CubicSplineFit::ValueHandling::Linear);
  fAcceptanceUnfolding->SetSplineXPositionsFlags(Math::CubicSplineFit::FixAllSplineXNodes|
                                                 Math::CubicSplineFit::FixLeftBeginValue|Math::CubicSplineFit::FixRightEndValue|
                                                 Math::CubicSplineFit::FixLeftDerivative|Math::CubicSplineFit::FixRightDerivative);
  fAcceptanceUnfolding->Run();

  if (verbosity >= 1 && !gROOT->IsBatch()) {
    fAcceptanceUnfolding->Plot();
    fIsPlotted = true;
  }
}

FoldedAcceptanceUnfoldingForLeptons::~FoldedAcceptanceUnfoldingForLeptons() {

  if (!fIsPlotted)
    delete fAcceptanceUnfolding;
  delete fAcceptanceForUnfolding;
  delete fRateToBeUnfolded;
}

TH1* FoldedAcceptanceUnfoldingForLeptons::MakeUnfoldedFlux() const {

  TH1* hFluxUnfolded = fAcceptanceUnfolding->Converged() ? (TH1*)fAcceptanceUnfolding->UnfoldedFluxAtLastIteration()->Clone(Form("hFluxUnfolded%d", ++fHistoCounter)) : nullptr;
  return hFluxUnfolded;
}

TH1* FoldedAcceptanceUnfoldingForLeptons::MakeFoldedUnfoldedFlux() const {

  TH1* hFluxUnfolded = MakeUnfoldedFlux();
  if (!hFluxUnfolded)
    return nullptr;

  TH2* hNormalizedMigrationMatrix = (TH2*)fMigrationMatrix->Clone("hNormalizedMigrationMatrix");
  Utilities::NormalizeHistogramYSlices(hNormalizedMigrationMatrix);

  double meastime_max = fMeasuringTime->GetBinContent(fMeasuringTime->GetNbinsX());

  TH1* hUnfoldedCountsBeforeCutoff = (TH1*)hFluxUnfolded->Clone("hUnfoldedCountsBeforeCutoff");
  hUnfoldedCountsBeforeCutoff->Reset();
  for (int ibin = 1; ibin <= hFluxUnfolded->GetNbinsX(); ++ibin) {
    double unfoldedCountsBeforeCutoff = hFluxUnfolded->GetBinContent(ibin) * fAcceptance->GetBinContent(ibin) * fTriggerEfficiency->GetBinContent(ibin) *
        hFluxUnfolded->GetXaxis()->GetBinWidth(ibin) * meastime_max;
    hUnfoldedCountsBeforeCutoff->SetBinContent(ibin, unfoldedCountsBeforeCutoff);
  }

  TH1* hMigratedFlux = (TH1*)hFluxUnfolded->Clone("hMigratedFluxFoldedAcc");
  hMigratedFlux->Reset();

  // migration matrix: x: rec, y: true
  // nSmeared(E_rec) = sum_(E_true) n(E_true) * prob(E_rec|E_true)
  for (int jbin = 1; jbin <= hNormalizedMigrationMatrix->GetNbinsY(); ++jbin) {

    double migratedCounts = 0.0;
    for (int ibin = 1; ibin <= hNormalizedMigrationMatrix->GetNbinsX(); ++ibin) {

      migratedCounts += hUnfoldedCountsBeforeCutoff->GetBinContent(ibin) * hNormalizedMigrationMatrix->GetBinContent(jbin, ibin);
    }
    hMigratedFlux->SetBinContent(jbin, migratedCounts);
  }

  hMigratedFlux->Scale(1.0/meastime_max);
  hMigratedFlux->Divide(fTriggerEfficiency);
  hMigratedFlux->Divide(fAcceptance);
  hMigratedFlux->Scale(1.0, "width");

  delete hUnfoldedCountsBeforeCutoff;
  delete hNormalizedMigrationMatrix;
  delete hFluxUnfolded;

  return hMigratedFlux;
}
