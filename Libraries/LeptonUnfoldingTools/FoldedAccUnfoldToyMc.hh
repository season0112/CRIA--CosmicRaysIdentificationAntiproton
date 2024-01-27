#ifndef FoldedAccUnfoldToyMc_hh
#define FoldedAccUnfoldToyMc_hh

#include <UnfoldingToyMc.hh>

class MigrationRandomGenerator;
class MigrationParameterization;
class TGraph;
class TH1;
class TH2;

class FoldedAccUnfoldToyMc : public UnfoldingToyMc {

public:
  FoldedAccUnfoldToyMc(const Binning::Definition& energyBinning, const MigrationParameterization* param, int randomSeed);
  virtual ~FoldedAccUnfoldToyMc();

  TH2* GenerateToyResponse(unsigned long long nToyResponse);

  void SimulateUnfoldedMeasurement(const TH1* hAssumedFlux, const TH1* hMeasuringTime, const TH1* hTriggerEfficiency, const TH1* hAcceptance,
                                   const TH2* toyResponse, double fractionOfMeasuringTime);

  TGraph* MakeFoldedOverMeasuredRatio() const;

protected:
  TH1* hToyFoldedUnfoldedFlux = nullptr;
};

#endif
