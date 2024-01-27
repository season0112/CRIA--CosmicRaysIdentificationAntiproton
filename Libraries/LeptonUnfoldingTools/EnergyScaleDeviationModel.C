#include "EnergyScaleDeviationModel.hh"

#include <TF1.h>

#include <cmath>

#define INFO_OUT_TAG "EnergyScaleDeviationModel"
#include "debugging.hh"

namespace Modelling {

EnergyScaleDeviationModel::EnergyScaleDeviationModel(double _offset, double _slope1, double _slope2, double _slope3, double _break1, double _break2)
  : Modelling::Model() {

  DefineParameter(offset = Modelling::ModelParameter("offset", _offset, 0.001));
  DefineParameter(slope1 = Modelling::ModelParameter("slope1", _slope1, 0.001));
  DefineParameter(slope2 = Modelling::ModelParameter("slope2", _slope2, 0.001));
  DefineParameter(slope3 = Modelling::ModelParameter("slope3", _slope3, 0.001));
  DefineParameter(break1 = Modelling::ModelParameter("break1", _break1, 0.01 * _break1));
  DefineParameter(break2 = Modelling::ModelParameter("break2", _break2, 0.01 * _break2));

  Deviation = StyleTF1(new TF1("EnergyScaleDeviation", this, &EnergyScaleDeviationModel::Function, 0.01, 10000.0, 0, "EnergyScaleDeviationModel", "Function"), kRed);
}

EnergyScaleDeviationModel::~EnergyScaleDeviationModel() {

  delete Deviation;
}

double EnergyScaleDeviationModel::Function(double* x, double*) {

  if (x[0] < break1)
    return slope1 * std::log10(x[0]) + (slope2 * std::log10(break1) + offset - slope1 * std::log10(break1));

  if (x[0] >= break1 && x[0] <= break2)
    return slope2 * std::log10(x[0]) + offset;

  if (x[0] > break2)
    return slope3 * std::log10(x[0]) + (slope2 * std::log10(break2) + offset - slope3 * std::log10(break2));

  return 0.0;
}

EnergyScaleDeviationModel* EnergyScaleDeviationModel::DefaultModel() {

  static Modelling::EnergyScaleDeviationModel* sModel = nullptr;
  if (!sModel)
    sModel = new Modelling::EnergyScaleDeviationModel(1.00905e+00,  1.01575e-02, -5.43975e-05, -1.98256e-02, 3.00000e+01, 1.50000e+02);
  return sModel;
}

}
