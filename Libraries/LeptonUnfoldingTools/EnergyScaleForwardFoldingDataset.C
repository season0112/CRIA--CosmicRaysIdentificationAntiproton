#include "EnergyScaleForwardFoldingDataset.hh"

#include "BinningTools.hh"
#include "FluxTools.hh"
#include "Utilities.hh"
#include "MatrixTools.hh"
#include "PredefinedBinnings.hh"

#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#define INFO_OUT_TAG "EnergyScaleForwardFoldingDataset"
#include "debugging.hh"

namespace Modelling {

EnergyScaleForwardFoldingDataset::EnergyScaleForwardFoldingDataset(const TH1* hEventCounts, const TH1* hAcceptance, const TH1* hTriggerEfficiency, const TH1* hMeasuringTime,
                                                                   TTree* migrationMatrixTree, bool useElectronWeights, const TF1* referenceModel)
  : Dataset()
  , fReferenceModel(referenceModel)
  , fEventCounts(hEventCounts)
  , fAcceptance(hAcceptance)
  , fTriggerEfficiency(hTriggerEfficiency)
  , fMeasuringTime(hMeasuringTime) {

  double mcEventWeight = 0.0;
  float mcGeneratedMomentum = 0.0f;
  float ecalEnergyBestMaximumShower = 0.0f;

  INFO_OUT << "Initializing migration matrix mini-tree..." << std::endl;
  migrationMatrixTree->SetBranchAddress("EcalEnergyBestMaximumShower", &ecalEnergyBestMaximumShower);
  migrationMatrixTree->SetBranchAddress("McGeneratedMomentum", &mcGeneratedMomentum);
  migrationMatrixTree->SetBranchAddress(useElectronWeights ? "McEventWeightElectron" : "McEventWeightPositron", &mcEventWeight);

  auto entries = migrationMatrixTree->GetEntries();
  for (decltype(entries) i = 0; i < entries; ++i) {
    migrationMatrixTree->GetEntry(i);
    fMcEventWeights.emplace_back(mcEventWeight);
    fMcGeneratedMomenta.emplace_back(mcGeneratedMomentum);
    fEcalEnergies.emplace_back(ecalEnergyBestMaximumShower);
  }

  INFO_OUT << "... done. Entries=" << migrationMatrixTree->GetEntries() << std::endl;

  if (!fEventCounts || fEventCounts->GetEntries() <= 0.0)
    FATAL_OUT << "Null pointer or counts histogram with zero entries used for dataset." << std::endl;
  if (!migrationMatrixTree || migrationMatrixTree->GetEntries() <= 0.0)
    FATAL_OUT << "Null pointer or empty migration matrix tree used for dataset." << std::endl;

  fBinning = Binning::Tools::FromTAxis(fEventCounts->GetXaxis());
  fName = std::string(fEventCounts->GetName());

  // test binning
  if (fAcceptance && fBinning != Binning::Tools::FromTAxis(fAcceptance->GetXaxis()))
    FATAL_OUT << "Binning mismatch with acceptance histogram." << std::endl;
  if (fTriggerEfficiency && fBinning != Binning::Tools::FromTAxis(fTriggerEfficiency->GetXaxis()))
    FATAL_OUT << "Binning mismatch with trigger efficiency histogram." << std::endl;
  if (fMeasuringTime && fBinning != Binning::Tools::FromTAxis(fMeasuringTime->GetXaxis()))
    FATAL_OUT << "Binning mismatch with MeasuringTime histogram." << std::endl;

  fXmin = fBinning.Min();
  fXmax = fBinning.Max();
}

EnergyScaleForwardFoldingDataset::~EnergyScaleForwardFoldingDataset() {

}

TH2* EnergyScaleForwardFoldingDataset::GenerateMigrationMatrix(const TF1* energyScaleDifferenceModel) const {

  const auto& analysisBinning = Binning::Predefined::AbsoluteEnergyBinning(); // Assumes this is the right analysis binning.
  TH2* migrationMatrix = Make<TH2D>("newMigrationMatrix", "", analysisBinning, analysisBinning);

  auto entries = fMcEventWeights.size();
  for (decltype(entries) i = 0; i < entries; ++i)
    migrationMatrix->Fill(fMcGeneratedMomenta[i], fEcalEnergies[i] * energyScaleDifferenceModel->Eval(fEcalEnergies[i]), fMcEventWeights[i]);

  migrationMatrix = MatrixTools::MakeMatrixStrippedOfOutermostBins(migrationMatrix);
  MatrixTools::Transpose(migrationMatrix);
  Utilities::NormalizeHistogramYSlices(migrationMatrix, "", true);

  if (fBinning != Binning::Tools::FromTAxis(migrationMatrix->GetXaxis()) || fBinning != Binning::Tools::FromTAxis(migrationMatrix->GetYaxis()))
    FATAL_OUT << "Binning mismatch with migration matrix." << std::endl;

  return migrationMatrix;
}

std::vector<double> EnergyScaleForwardFoldingDataset::ForwardFoldedCounts(const TF1* energyScaleDifferenceModel) const {

  std::vector<double> modelIntegralCounts(fBinning.NumberOfBins(), 0.0);
  for (unsigned int bin = 1; bin <= fBinning.NumberOfBins(); ++bin) {

    double X1 = fBinning.LowEdge(bin);
    double X2 = fBinning.UpEdge(bin);
    modelIntegralCounts[bin - 1] = IntegratePowerLaw(fReferenceModel, X1, X2);
  }

  double fullMeasurementTime = fMeasuringTime ? fMeasuringTime->GetBinContent(fMeasuringTime->GetNbinsX()) : 1.0;

  // Calculate counts before convolution with the migration matrix:
  // Multiply flux by acceptance, trigger efficiency and measuring time (if given).
  for (unsigned int bin = 1; bin <= fBinning.NumberOfBins(); ++bin) {

    if (fAcceptance)
      modelIntegralCounts[bin - 1] *= fAcceptance->GetBinContent(bin);
    if (fTriggerEfficiency)
      modelIntegralCounts[bin - 1] *= fTriggerEfficiency->GetBinContent(bin);
    if (fMeasuringTime)
      modelIntegralCounts[bin - 1] *= fullMeasurementTime;
  }

  // forward-folding of counts
  std::vector<double> forwardFoldedCounts(fBinning.NumberOfBins(), 0.0);

  auto* newMigrationMatrix = GenerateMigrationMatrix(energyScaleDifferenceModel);

  for (unsigned int xbin = 1; xbin <= fBinning.NumberOfBins(); ++xbin) {

    double migratedCounts = 0.0;
    for (unsigned int i = 1; i <= fBinning.NumberOfBins(); ++i)
      migratedCounts += modelIntegralCounts[i - 1] * newMigrationMatrix->GetBinContent(xbin, i);

    // account for geomagnetic cutoff
    double cutoffSurvivalProb = fMeasuringTime ? fMeasuringTime->GetBinContent(xbin) / fullMeasurementTime : 1.0;
    forwardFoldedCounts[xbin - 1] = migratedCounts * cutoffSurvivalProb;
  }

  delete newMigrationMatrix;
  return forwardFoldedCounts;
}

// re-implement Chi2 function for efficiency
double EnergyScaleForwardFoldingDataset::Chi2(const TF1* energyScaleDifferenceModel) const {

  std::vector<double> modelCounts = ForwardFoldedCounts(energyScaleDifferenceModel);

  fNPointsInFit = 0;
  double chi2 = 0.0;

  for (int i = FirstIndexForLoops(); i <= LastIndexForLoops(); ++i) {

    if (IsInRange(X(i))) {

      double dataCounts = fEventCounts->GetBinContent(i);
      double dataCountsError = fEventCounts->GetBinError(i);
      double pull = (dataCounts - modelCounts.at(i - 1)) / dataCountsError;

      chi2 += pull * pull;
      ++fNPointsInFit;
    }
  }

  return chi2;
}

double EnergyScaleForwardFoldingDataset::Pull(const TF1* energyScaleDifferenceModel, int bin) const {

  const std::vector<double>& modelCounts = ForwardFoldedCounts(energyScaleDifferenceModel);

  double dataCounts = fEventCounts->GetBinContent(bin);
  double dataCountsError = fEventCounts->GetBinError(bin);

  double pull = (dataCounts - modelCounts.at(bin - 1)) / dataCountsError;

  if (fVerbosity) {
    double X1 = fBinning.LowEdge(bin);
    double X2 = fBinning.UpEdge(bin);
    std::cout << "ds " << fIdentifier << ": bin= " << bin << " X1 = " << X1 << " X2 = " << X2
              << " data = " << dataCounts << " +- " << dataCountsError
              << " model= " << modelCounts.at(bin - 1)
              << " -> pull " << pull
              << " chi2 contrib " << std::pow(pull, 2) << std::endl;
  }

  return pull;
}

int EnergyScaleForwardFoldingDataset::LastIndexForLoops() const {

  return fBinning.NumberOfBins();
}

double EnergyScaleForwardFoldingDataset::X(int bin) const {

  return fBinning.BinCenterLog(bin);
}

double EnergyScaleForwardFoldingDataset::Value(double x) const {

  auto bin = fEventCounts->GetXaxis()->FindFixBin(x);
  double fluxValue = fEventCounts->GetBinContent(bin);
  if (fAcceptance)
    fluxValue /= fAcceptance->GetBinContent(bin);
  if (fTriggerEfficiency)
    fluxValue /= fTriggerEfficiency->GetBinContent(bin);
  if (fMeasuringTime)
    fluxValue /= fMeasuringTime->GetBinContent(bin);
  fluxValue /= fEventCounts->GetXaxis()->GetBinWidth(bin);

  return fluxValue;
}

}
