#ifndef UnfoldingSystematicErrors_hh
#define UnfoldingSystematicErrors_hh

#include <string>

class TH2;
class TGraphAsymmErrors;

TGraphAsymmErrors* ConstructUnfoldedAllElectronEventCount(TGraphAsymmErrors* electronEventCountUnfolded, TGraphAsymmErrors* positronEventCountUnfolded);
TGraphAsymmErrors* UnfoldEventCounts(const std::string& inputFileSuffix, TGraphAsymmErrors* eventCounts, TGraphAsymmErrors* measuringTime, TH2* migrationMatrix, const std::string& speciesSymbol, bool isTimeAveragedFlux, TGraphAsymmErrors* effectiveAcceptance = nullptr, TGraphAsymmErrors* triggerEfficiency = nullptr, TGraphAsymmErrors* ecalBDTEfficiency = nullptr);

#endif
