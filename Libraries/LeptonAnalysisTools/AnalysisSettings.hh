#ifndef AnalysisSettings_hh
#define AnalysisSettings_hh

#include <string>
#include <vector>

#include "BinningTools.hh"
#include "Environment.hh"
#include "PredefinedBinnings.hh"

class TCanvas;

extern const long gAMSFirstEvent;
extern const long gAMSLastEvent;
extern const long gAMSLastEventPass4;
extern const long gAMSLastEventPass6;
extern const char* gAMSDataPeriod;
extern const char* gAMSDataPeriodPass4;
extern const char* gAMSDataPeriodPass6;

extern const double gAMSMinEnergy;
extern const double gAMSMaxEnergy;

extern const double gAMSStartShowEnergy; // >= gAMSMinEnergy (for visual purposes)
extern const double gAMSStopShowEnergy; // <= gAMSMaxEnergy (for visual purposes)

extern const double gAMSTimeDependentFluxesMaxEnergy;
extern const double gAMSTimeDependentFluxesMinEnergy;
extern const double gAMSLaffertyWyattSpectralIndex;

static const unsigned int gTrackerPattern = 6;
extern const int gTrackerPatternColors[gTrackerPattern];

extern const int gElectronColor;
extern const int gCCElectronColor;
extern const int gPositronColor;
extern const int gCCPositronColor;
extern const int gProtonColor;
extern const int gCCProtonColor;
extern const int gCCColor;
extern const int gFitResultColor;

class AnalysisSettings {
public:
  static void Initialize();

  static void SwitchToQuadraticCanvas();
  static void SwitchToRectangularCanvas();
};

// FIXME: Move these to their own file.
void SaveCanvas(TCanvas*, std::string fileName);
void MoveTitleBoxToCenter();
int TransparentColor(int color, float opacity = 0.9);

bool IsBlacklistedBartelsRotation(unsigned int bartelsRotation);
bool IsBlacklistedBartelsRotation(const std::string& inputFileSuffix);

int ClassifyEnergyIntoBin(double energy);
int ClassifyEnergyIntoBin(double energy, const Binning::Definition& binning);

#endif // AnalysisSettings_hh
