#ifndef UnfoldingToyMc_hh
#define UnfoldingToyMc_hh

#include <BinningDefinition.hh>

class MigrationRandomGenerator;
class MigrationParameterization;
class TRandom;
class TGraphErrors;
class TH1;

class UnfoldingToyMc {

public:
  UnfoldingToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, int randomSeed);
  virtual ~UnfoldingToyMc();

  virtual void SetVerbosity(int val) { fVerbosity = val; }

  TH1* ToyAverageEventCounts() { return hToyAverageEventCounts; }
  TH1* ToyTruthEventCounts() { return hToyTruthEventCounts; }
  TH1* ToyMeasEventCounts() { return hToyMeasEventCounts; }
  TH1* ToyUnfoldedEventCounts() { return hToyUnfoldedEventCounts; }

  TH1* ToyAverageFlux() { return hToyAverageFlux; }
  TH1* ToyTruthFlux() { return hToyTruthFlux; }
  TH1* ToyMeasFlux() { return hToyMeasFlux; }
  TH1* ToyUnfoldedFlux() { return hToyUnfoldedFlux; }

  TGraphErrors* MakeMeasuredOverTruthRatio() const;
  TGraphErrors* MakeUnfoldedOverTruthRatio() const;
  TGraphErrors* MakeUnfoldedOverAverageRatio() const;
  TGraphErrors* MakeTruthOverAverageRatio() const;

protected:
  virtual void DeleteHistograms();

protected:

  TRandom* fRandom;
  Binning::Definition fEnergyBinning;
  const MigrationParameterization* fParam = nullptr;
  MigrationRandomGenerator* fRandomGenerator = nullptr;
  int fVerbosity = 1;

  // results, filled by SimulateUnfoldedMeasurement()

  /** Average expected event counts for given measurement time, i.e. no smearing applied, before cutoff, as a function of true energy. */
  TH1* hToyAverageEventCounts = nullptr;
  /** Simulated event counts, before cutoff, as a function of true energy. */
  TH1* hToyTruthEventCounts = nullptr;
  /** Simulated event counts, after cutoff, as a function of reconstructed energy. */
  TH1* hToyMeasEventCounts = nullptr;
  /** Simulated event counts, after unfolding. */
  TH1* hToyUnfoldedEventCounts = nullptr;

  TH1* hToyAverageFlux = nullptr;
  TH1* hToyTruthFlux = nullptr;
  TH1* hToyMeasFlux = nullptr;
  TH1* hToyUnfoldedFlux = nullptr;
};

#endif
