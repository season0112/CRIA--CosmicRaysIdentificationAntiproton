#ifndef EnergyScaleDeviationModel_hh
#define EnergyScaleDeviationModel_hh

#include <Model.hh>
#include <ModelParameter.hh>

class TF1;

namespace Modelling {

class EnergyScaleDeviationModel : public Modelling::Model {
public:
  EnergyScaleDeviationModel(double _offset, double _slope1, double _slope2, double _slope3, double _break1, double _break2);
  virtual ~EnergyScaleDeviationModel();

  static EnergyScaleDeviationModel* DefaultModel();

  virtual TF1* GetPredictionForIdentifier(int) const { return Deviation; }

  double Function(double* x, double*);
  TF1* Deviation;

  Modelling::ModelParameter offset;
  Modelling::ModelParameter slope1;
  Modelling::ModelParameter slope2;
  Modelling::ModelParameter slope3;
  Modelling::ModelParameter break1;
  Modelling::ModelParameter break2;
};

};

#endif
