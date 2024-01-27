#ifndef RooUnfoldToyMc_hh
#define RooUnfoldToyMc_hh

#include <UnfoldingToyMc.hh>

class MigrationRandomGenerator;
class MigrationParameterization;
class RooUnfoldResponse;
class TH1;
class TH2D;

class RooUnfoldToyMc : public UnfoldingToyMc {

public:
  RooUnfoldToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, int randomSeed);
  virtual ~RooUnfoldToyMc();

  RooUnfoldResponse* GenerateToyResponseNoCutoff(unsigned long long int nToyResponse);
  RooUnfoldResponse* GenerateToyResponseWithCutoff(unsigned long long nToyResponse, const TH1* hMeasuringTime);

  void SimulateUnfoldedMeasurement(const TH1* hAssumedFlux, const TH1* hMeasuringTime, const TH1* hTriggerEfficiency, const TH1* hAcceptance,
                                   RooUnfoldResponse* toyResponse, double fractionOfMeasuringTime, bool simulateCutoff);

  void SetBayesIterations(int val) { fBayesIterations = val; }

  TH2D* MigrationMatrixBeforeCutoff() { return hMigrationMatrixBeforeCutoff; }

private:

  int fBayesIterations = 10;

  /** Migration matrix, filled in parallel to RooResponse object, but ignoring cutoff effects. x: Emeas, y: Etrue. */
  TH2D* hMigrationMatrixBeforeCutoff = nullptr;

};

#endif
