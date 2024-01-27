#include "ModelFits.hh"

#include <BrokenPowerLawModel.hh>
#include <ModellingData.hh>
#include <Model.hh>
#include <ModelFunctions.hh>
#include <FitFunction.hh>
#include <ModelAnalysis.hh>
#include <ModelParameter.hh>
#include <Dataset.hh>
#include <ParticleId.hh>

#include "EnergyScaleDeviationModel.hh"

#define INFO_OUT_TAG "ModelFits"
#include <debugging.hh>

Modelling::BrokenPowerLawModel* RunBrokenPowerLawFit(const char* species_symbol, Modelling::Dataset* ds, bool drawAuxCanvases, bool timedependentCase, int verbosity) {

  Modelling::Data* dataFit = new Modelling::Data;
  dataFit->AddDataset(0, ds);

  Modelling::BrokenPowerLawModel* model = nullptr;
  if (std::string(species_symbol) == "e-") {
    model = new Modelling::BrokenPowerLawModel(ParticleId::Mass(ParticleId::Electron), 2.75,
    { {5.0, 3.65, 0.25}, {30.0, 3.25, 0.3} } );
  }
  else {
    model = new Modelling::BrokenPowerLawModel(ParticleId::Mass(ParticleId::Positron), 2.38,
    { {3.55, 3.15, 0.25}, {27.0, 2.74, 0.3}, {500.0, 3.74, 0.4} } );
    model->phi.SetValue(0.33);
    model->BreakEnergy(3).Fix();
    model->BreakSmoothness(3).Fix();
  }
  model->BreakSmoothness(1).Fix();
  model->BreakSmoothness(2).Fix();

  if (timedependentCase) {

    if (std::string(species_symbol) == "e-") {
      model->BreakEnergy(2).Fix();
    }

    if (std::string(species_symbol) == "e+") {
      model->BreakEnergy(3).Fix();
      model->BreakSmoothness(3).Fix();
      model->SpectralIndexAfterBreak(3).Fix();
    }
  }

  model->AutosetStartValueForAmplitude(ds);
  Modelling::FitFunction* fitf = new Modelling::FitFunction(dataFit, model);
  Modelling::ModelAnalysis analysis(fitf, verbosity);
  // TEMP
  analysis.RunAnalysis();

  if (drawAuxCanvases) {
    static int sCtr = 0;
    analysis.MakePullCanvas(Form("pullCanvas_%d", sCtr));
    analysis.MakeCumulativeChi2Canvas(Form("cumulativeChi2Canvas_%d", sCtr));
    ++sCtr;
  }

  return model;
}

Modelling::EnergyScaleDeviationModel* RunEnergyScaleDeviationModelFit(Modelling::Dataset* ds, bool drawAuxCanvases, int verbosity) {

  Modelling::Data* dataFit = new Modelling::Data;
  dataFit->AddDataset(0, ds);

  Modelling::EnergyScaleDeviationModel* model = new Modelling::EnergyScaleDeviationModel(1.01, 0.01, 0.0, -0.02, 30.0, 150.0);
  model->break1.SetLimits(1, 30);
  model->break2.SetLimits(30, 300);

  model->break1.Fix();
  model->break2.Fix();

  Modelling::FitFunction* fitf = new Modelling::FitFunction(dataFit, model);
  Modelling::ModelAnalysis analysis(fitf, verbosity);
  analysis.RunAnalysis();

  if (drawAuxCanvases) {
    static int sCtr = 0;
    analysis.MakePullCanvas(Form("pullCanvasEneScale_%d", sCtr));
    analysis.MakeCumulativeChi2Canvas(Form("cumulativeChi2CanvasEneScale_%d", sCtr));
    ++sCtr;
  }

  return model;
}
