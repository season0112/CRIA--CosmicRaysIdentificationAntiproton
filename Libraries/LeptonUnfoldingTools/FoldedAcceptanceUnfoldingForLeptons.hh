#ifndef FoldedAcceptanceUnfoldingForLeptons_hh_
#define FoldedAcceptanceUnfoldingForLeptons_hh_

class TH1;
class TH2;

namespace Unfolding {
class AcceptanceUnfolding;
}

class FoldedAcceptanceUnfoldingForLeptons {

public:

  /** Perform FoldedAcceptance unfolding for lepton fluxes.
    *
    * \param hEventCounts Histogram of event counts as a function of energy.
    * \param hMeasuringTime Histogram of live-time weighted exposure time, after accounting for geomagnetic cutoff, in same binning as \p hEventCounts histogram.
    * \param hAcceptance Effective acceptance, as a function of true energy, in same binning as \p hEventCounts histogram.
    * \param hTriggerEfficiency Trigger efficiency, as a function of true energy, in same binning as \p hEventCounts histogram.
    * \param hMigrationMatrix Migration matrix, with identical binning on both axes, same as \p hEventCounts histogram. Binning convention: reconstructed energy on x-axis, true energy on y-axis.
    * \param verbosity Verbosity level for unfolding (and plotting).
    *
    * \attention Take note of the binning convention for the migration matrix: reconstructed energy on x-axis, true energy on y-axis!
    */
  FoldedAcceptanceUnfoldingForLeptons(const TH1* hEventCounts, const TH1* hMeasuringTime, const TH1* hAcceptance, const TH1* hTriggerEfficiency,
                                      const TH2* hMigrationMatrix, int verbosity = 0);

  ~FoldedAcceptanceUnfoldingForLeptons();

  /** Make a copy of the unfolded flux histogram. */
  TH1* MakeUnfoldedFlux() const;

  /** Make a copy of the flux that was unfolded and convoluted again with the migration matrix. For closure test. */
  TH1* MakeFoldedUnfoldedFlux() const;

private:

  Unfolding::AcceptanceUnfolding* fAcceptanceUnfolding = nullptr;
  TH1* fAcceptanceForUnfolding = nullptr;
  TH1* fRateToBeUnfolded = nullptr;

  const TH2* fMigrationMatrix = nullptr;
  const TH1* fAcceptance = nullptr;
  const TH1* fEventCounts = nullptr;
  const TH1* fMeasuringTime = nullptr;
  const TH1* fTriggerEfficiency = nullptr;

  bool fIsPlotted = false;
  mutable int fHistoCounter = 0;
};

#endif
