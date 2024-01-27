#ifndef BayesUnfoldingWithCutoff_hh_
#define BayesUnfoldingWithCutoff_hh_

class TH1;
class TH1D;
class TH2;
class TH2D;
class RooUnfoldBayes;
class RooUnfoldResponse;

class BayesUnfoldingWithCutoff {

public:

  /** Perform Bayesian unfolding.
    *
    * \param hEventCounts Histogram of event counts as a function of energy.
    * \param hMeasuringTime Histogram of live-time weighted exposure time, after accounting for geomagnetic cutoff, in same binning as \p hEventCounts histogram.
    * \param hAcceptance Effective acceptance, as a function of true energy, in same binning as \p hEventCounts histogram.
    * \param hTriggerEfficiency Trigger efficiency, as a function of true energy, in same binning as \p hEventCounts histogram.
    * \param hMigrationMatrix Migration matrix, with identical binning on both axes, same as \p hEventCounts histogram. Binning convention: reconstructed energy on x-axis, true energy on y-axis.
    * \param nBayesIter Number of iterations for unfolding.
    * \param verbosity Verbosity level for unfolding.
    *
    * \attention Take note of the binning convention for the migration matrix: reconstructed energy on x-axis, true energy on y-axis!
    */
  BayesUnfoldingWithCutoff(const TH1* hEventCounts, const TH1* hMeasuringTime, const TH1* hAcceptance, const TH1* hTriggerEfficiency,
                           const TH2* hMigrationMatrix, double MaxMeasuringTime,
                           int nBayesIter = 4, int verbosity = 0);

  ~BayesUnfoldingWithCutoff();

  TH1D* UnfoldedEventCounts() const;
  TH2* MigrationMatrixAfterCutoff() const { return fMigrationMatrixAfterCutoff; }

  /** Make a copy of the unfolded flux histogram. */
  TH1D* MakeUnfoldedFlux() const;

  TH2D* MakeCountsCovarianceMatrixHistogram() const;

  /** Make a copy of the flux that was unfolded and convoluted again with the migration matrix. For closure test. */
  TH1D* MakeFoldedUnfoldedFlux() const;

  RooUnfoldResponse* Response() const { return fResponse; }

private:
  RooUnfoldResponse* fResponse = nullptr;
  RooUnfoldBayes* fUnfold = nullptr;
  double fMaximumMeasuringTime = 0.0;

  TH2* fMigrationMatrixAfterCutoff = nullptr;
  const TH2* fMigrationMatrix = nullptr;
  const TH1* fAcceptance = nullptr;
  const TH1* fEventCounts = nullptr;
  const TH1* fMeasuringTime = nullptr;
  const TH1* fTriggerEfficiency = nullptr;
};

#endif
