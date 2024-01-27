#ifndef SystematicUncertainties_hh
#define SystematicUncertainties_hh

#include <string>
#include <vector>

#include "BinningDefinition.hh"

class TF1;
class TH1D;
class TH2D;
class TGraph;
class TGraphAsymmErrors;
class TSpline5;

TGraph* AdditionalLowEnergySystematicUncertaintyGraph(TGraph*, TGraphAsymmErrors*);
void AddAdditionalRelativeSystematicUncertaintyToGraph(TGraphAsymmErrors*, TGraph* additionalSystematicUncertainty);

// Shared code between PlotFluxUncertainties.hh / PlotFluxUncertaintiesForBartelsRotation.hh
enum SpeciesType {
  ElectronFlux,
  PositronFlux,
  AllElectronFlux,
  PositronFraction,
  PositronElectronRatio
};

std::string SpeciesTitleForType(SpeciesType type);
std::string SpeciesNameForType(SpeciesType type);

void DumpTableEntry(unsigned int bin, double x, double y, double ey);
void DumpTableEntryInvalid(unsigned int bin, double x, double y, double ey, const std::string& problem);

TH1D* DrawUncertainty(SpeciesType speciesType, TGraph* uncertaintyGraph, bool same = false, double yLower = 1e-5, double yUpper = 1);
TH1D* DrawUncertaintyInBartelsRotation(SpeciesType speciesType, TGraph* uncertaintyGraph, bool same = false, double yLower = 1e-5, double yUpper = 1);
void DrawSingleUncertainty(TGraph* uncertainty, SpeciesType speciesType, const std::string& systName, const std::string& systTitle, std::string fileSuffix = "", double yLower = 1e-5, double yUpper = 1);

TGraph* ExtractTotalRelativeUncertaintyFromGraph(TGraphAsymmErrors* graph);
TGraph* CombineRelativeUncertainties(TGraphAsymmErrors* graph, const std::vector<TGraph*>& uncertainties);

void VerifyRelativeUncertaintyGraphIsNull(TGraph* graph, const std::string& description);
void VerifyRelativeUncertaintyGraphsAreEqual(TGraph* graph1, TGraph* graph2, const std::string& description, double tolerance = 1e-4);

TGraph* DetermineMigrationSystematicUncertainty(SpeciesType speciesType, TH2D* migrationMatrix, TF1* electronFluxModel, TF1* positronFluxModel, TF1* measuringTimeFunction, TGraphAsymmErrors* measuringTime, TSpline5* mcEffectiveAcceptanceFunction, TSpline5* triggerEfficiencyFunction);

#endif
