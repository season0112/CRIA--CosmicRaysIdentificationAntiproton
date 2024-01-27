#include "PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff.hh"

#include <TF1.h>

#include "FitFunction.hh"
#include "HistogramDataset.hh"
#include "ModelAnalysis.hh"
#include "ModelFunctions.hh"
#include "ModellingData.hh"
#include "Utilities.hh"

PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff(double _Epower)
  : Modelling::Model()
  , Epower(_Epower) {

  DefineParameter(phiPlus = Modelling::ModelParameter("#phi_{+}", 1.1, 0.01));

  DefineParameter(C_d = Modelling::ModelParameter("C_{d}", 6.51e-2, 0.01));
  DefineParameter(E_1 = Modelling::ModelParameter("E_{1}", 7.0, 0.1));
  DefineParameter(Gamma_d = Modelling::ModelParameter("#gamma_{d}", -4.07, 0.01));

  DefineParameter(C_s = Modelling::ModelParameter("C_{s}", 6.8e-5, 0.01));
  DefineParameter(E_2 = Modelling::ModelParameter("E_{2}", 60.0, 0.1));
  DefineParameter(invE_s = Modelling::ModelParameter("invE_{s}", 1.23 / 1000.0, 0.1));
  DefineParameter(Gamma_s = Modelling::ModelParameter("#gamma_{s}", -2.58, 0.01));

  E_1.Fix();
  E_2.Fix();

  double Emin = 0.5;
  double Emax = 1500.;

  PosiFlux = StyleTF1(new TF1("PosiFlux", this, &PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxF, Emin, Emax, 0, "PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff", "PosiFluxF"), kGreen + 3);
  PosiFluxPowerLawD = StyleTF1(new TF1("PosiFluxPowerLawD", this, &PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawDF, Emin, Emax, 0, "PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff", "PosiFluxPowerLawDF"), kGray);
  PosiFluxPowerLawS = StyleTF1(new TF1("PosiFluxPowerLawS", this, &PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawSF, Emin, Emax, 0, "PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff", "PosiFluxPowerLawSF"), kMagenta - 3);
  PosiFluxEpower = StyleTF1(new TF1("PosiFluxEpower", this, &PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxEpowerF, Emin, Emax, 0, "PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff", "PosiFluxEpowerF"), kGreen + 3);
  PosiFluxPowerLawDEpower = StyleTF1(new TF1("PosiFluxPowerLawDEpower", this, &PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawDEpowerF, Emin, Emax, 0, "PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff", "PosiFluxPowerLawDEpowerF"), kGray);
  PosiFluxPowerLawSEpower = StyleTF1(new TF1("PosiFluxPowerLawSEpower", this, &PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawSEpowerF, Emin, Emax, 0, "PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff", "PosiFluxPowerLawSEpowerF"), kMagenta - 3);

  PosiFlux->SetNpx(1e4);
  PosiFluxPowerLawD->SetNpx(1e4);
  PosiFluxPowerLawS->SetNpx(1e4);
  PosiFluxEpower->SetNpx(1e4);
  PosiFluxPowerLawDEpower->SetNpx(1e4);
  PosiFluxPowerLawSEpower->SetNpx(1e4);
}

double PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxF(double* x, double* par) {

  return PosiFluxPowerLawDF(x, par) + PosiFluxPowerLawSF(x, par);
}

double PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawDF(double* x, double*) {

  double E = x[0];
  double Ehat = E + phiPlus;
  return Modelling::SolarModulationTerm(E, phiPlus, M) * C_d * std::pow(Ehat / E_1, Gamma_d);
}

double PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawSF(double* x, double*) {

  double E = x[0];
  double Ehat = E + phiPlus;
  return Modelling::SolarModulationTerm(E, phiPlus, M) * C_s * std::pow(Ehat / E_2, Gamma_s) * std::exp(-Ehat * invE_s);
}

double PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * PosiFluxF(x, par);
}

double PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawDEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * PosiFluxPowerLawDF(x, par);
}

double PositronModelTwoPowerLawsWithSolarModulationAndExpCutoff::PosiFluxPowerLawSEpowerF(double* x, double* par) {

  return std::pow(x[0], Epower) * PosiFluxPowerLawSF(x, par);
}
