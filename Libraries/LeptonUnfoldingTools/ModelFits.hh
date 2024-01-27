#ifndef ModelFits_hh_
#define ModelFits_hh_

namespace Modelling {
class BrokenPowerLawModel;
class EnergyScaleDeviationModel;
class Dataset;
}

Modelling::BrokenPowerLawModel* RunBrokenPowerLawFit(const char* species_symbol, Modelling::Dataset* ds, bool drawAuxCanvases, bool timedependentCase,
                                                     int verbosity = 1);

Modelling::EnergyScaleDeviationModel* RunEnergyScaleDeviationModelFit(Modelling::Dataset* ds, bool drawAuxCanvases, int verbosity = 1);

#endif
