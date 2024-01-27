#ifndef EnergyScaleForwardFoldingDataset_hh
#define EnergyScaleForwardFoldingDataset_hh

#include "Dataset.hh"

class TF1;
class TH1;
class TH2;
class TTree;

namespace Modelling {

class EnergyScaleForwardFoldingDataset : public Dataset {

public:
  EnergyScaleForwardFoldingDataset(const TH1* hEventCounts, const TH1* hAcceptance, const TH1* hTriggerEfficiency, const TH1* hMeasuringTime,
                                   TTree* migrationMatrixTree, bool useElectronWeights, const TF1* referenceModel);
  virtual ~EnergyScaleForwardFoldingDataset();

  virtual double Pull(const TF1* energyScaleDifferenceModel, int bin) const;
  virtual double Chi2(const TF1* energyScaleDifferenceModel) const;

  virtual int FirstIndexForLoops() const { return 1; }
  virtual int LastIndexForLoops() const;
  virtual double X(int bin) const;
  virtual double Value(double x) const;

  TH2* GenerateMigrationMatrix(const TF1* energyScaleDifferenceModel) const;

private:
  std::vector<double> ForwardFoldedCounts(const TF1* energyScaleDifferenceModel) const;

  const TF1* fReferenceModel;
  const TH1* fEventCounts;
  const TH1* fAcceptance;
  const TH1* fTriggerEfficiency;
  const TH1* fMeasuringTime;

  Binning::Definition fBinning;

  std::vector<double> fMcEventWeights;
  std::vector<float> fMcGeneratedMomenta;
  std::vector<float> fEcalEnergies;
};

}

#endif
