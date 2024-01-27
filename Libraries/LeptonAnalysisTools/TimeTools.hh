#ifndef TimeTools_hh
#define TimeTools_hh

#include <cassert>
#include <ctime>
#include <cstdlib>
#include <sstream>

#include <TColor.h>
#include <TDatime.h>
#include <TLine.h>
#include <TROOT.h>
#include <TPad.h>
#include <TPaveText.h>

#include "AnalysisSettings.hh"

unsigned int BartelsRotationToMonth(unsigned int bartelsRotation) {

  const auto& binning = Binning::Tools::BartelsRotationTimeBinning(gAMSFirstEvent, gAMSLastEvent);
  assert(bartelsRotation <= binning.NumberOfBins());

  time_t startTime_as_time_t = binning.LowEdge(bartelsRotation + 1);
  struct std::tm* timeInfo = std::gmtime(&startTime_as_time_t);

  return timeInfo->tm_mon + 1;
}

unsigned int BartelsRotationToYear(unsigned int bartelsRotation) {

  const auto& binning = Binning::Tools::BartelsRotationTimeBinning(gAMSFirstEvent, gAMSLastEvent);
  assert(bartelsRotation <= binning.NumberOfBins());

  time_t startTime_as_time_t = binning.LowEdge(bartelsRotation + 1);
  struct std::tm* timeInfo = std::gmtime(&startTime_as_time_t);

  return timeInfo->tm_year + 1900;
}

std::string BartelsRotationToMonthDayString(unsigned int bartelsRotation) {

  const auto& binning = Binning::Tools::BartelsRotationTimeBinning(gAMSFirstEvent, gAMSLastEvent);
  assert(bartelsRotation <= binning.NumberOfBins());

  time_t startTime_as_time_t = binning.LowEdge(bartelsRotation + 1);
  struct std::tm* timeInfo = std::gmtime(&startTime_as_time_t);

  char buffer[80];
  std::strftime(buffer, 80, "%b %d", timeInfo);
  return std::string(buffer);
}

std::string BartelsRotationToDateTimeString(unsigned int bartelsRotation) {

  const auto& binning = Binning::Tools::BartelsRotationTimeBinning(gAMSFirstEvent, gAMSLastEvent);
  assert(bartelsRotation <= binning.NumberOfBins());

  time_t startTime_as_time_t = binning.LowEdge(bartelsRotation + 1);
  struct std::tm* timeInfo = std::gmtime(&startTime_as_time_t);

  char buffer[80];
  std::strftime(buffer, 80, "%c", timeInfo);
  return std::string(buffer);
}

std::string BartelsRotationToMonthYearString(unsigned int bartelsRotation) {

  const auto& binning = Binning::Tools::BartelsRotationTimeBinning(gAMSFirstEvent, gAMSLastEvent);
  assert(bartelsRotation <= binning.NumberOfBins());

  time_t startTime_as_time_t = binning.LowEdge(bartelsRotation + 1);
  struct std::tm* timeInfo = std::gmtime(&startTime_as_time_t);

  char buffer[80];
  std::strftime(buffer, 80, "%B %Y", timeInfo);
  return std::string(buffer);
}

int BartelsRotationToUnixTime(unsigned int bartelsRotation) {

  const auto& binning = Binning::Tools::BartelsRotationTimeBinning(gAMSFirstEvent, gAMSLastEvent);
  assert(bartelsRotation <= binning.NumberOfBins());
  return binning.LowEdge(bartelsRotation + 1);
}

// Expects a string in the form "Bartels_XY"
bool BartelsRotationStringToBeginAndEndTimes(const std::string& bartelsRotationStringRaw, int& bartelsRotationBeginTime, int& bartelsRotationEndTime) {

  std::string bartelsRotationString = bartelsRotationStringRaw.substr(8);

  unsigned int bartelsRotation = 0;
  std::stringstream conversionStream;
  conversionStream << bartelsRotationString;
  if (!(conversionStream >> bartelsRotation))
    return false;

  bartelsRotationBeginTime = BartelsRotationToUnixTime(bartelsRotation);
  bartelsRotationEndTime = BartelsRotationToUnixTime(bartelsRotation + 1);
  return true;
}

int ElectronColorForBartelsRotation(unsigned int firstBartelsRotation, unsigned int lastBartelsRotation) {

  const unsigned int sNumberOfStops = 2;
  double sRedElectron[sNumberOfStops]   = { 1.0, 1.0 };
  double sGreenElectron[sNumberOfStops] = { 0.0, 1.0 };
  double sBlueElectron[sNumberOfStops]  = { 0.0, 0.0 };

  double sStops[sNumberOfStops] = { 0.0, 1.0 };

  // This assumes first/lastBartelsRotation never change.
  static int sPalette = -1;
  if (sPalette == -1) {
    int numberOfColors = (lastBartelsRotation - firstBartelsRotation) + 1;
    sPalette = TColor::CreateGradientColorTable(sNumberOfStops, sStops, sRedElectron, sGreenElectron, sBlueElectron, numberOfColors);
  }

  return sPalette;
}

int PositronColorForBartelsRotation(unsigned int firstBartelsRotation, unsigned int lastBartelsRotation) {

  const unsigned int sNumberOfStops = 2;
  double sRedPositron[sNumberOfStops]   = { 0.0, 0.0 };
  double sGreenPositron[sNumberOfStops] = { 0.0, 1.0 };
  double sBluePositron[sNumberOfStops]  = { 1.0, 1.0 };

  double sStops[sNumberOfStops] = { 0.0, 1.0 };

  // This assumes first/lastBartelsRotation never change.
  static int sPalette = -1;
  if (sPalette == -1) {
    int numberOfColors = (lastBartelsRotation - firstBartelsRotation) + 1;
    sPalette = TColor::CreateGradientColorTable(sNumberOfStops, sStops, sRedPositron, sGreenPositron, sBluePositron, numberOfColors);
  }

  return sPalette;
}

void CreateYearIndicators(const char* text, int startBartelsRotation, int bartelsRotationsToSpan, int totalBartelsRotations, double yOffset = 0.0) {

  static double textBoxHeight = 0.06;
  double textBoxStart = gPad->GetLeftMargin();
  double textBoxEnd = 1.0 - gPad->GetRightMargin();
  double distance = (textBoxEnd - textBoxStart) / totalBartelsRotations;
  double textBoxY1 = gPad->GetBottomMargin() / 2.0 - 0.005 + yOffset;
  double textBoxY2 = 0.0;
  TPaveText* textBox = new TPaveText(textBoxStart + startBartelsRotation * distance + 0.001, textBoxY1, textBoxStart + (startBartelsRotation + bartelsRotationsToSpan) * distance - 0.001, textBoxY2, "brNDC");
  textBox->AddText(text);
  textBox->SetBorderSize(0);
  textBox->SetFillColor(kWhite);
  textBox->SetTextColor(kBlack);
  textBox->SetTextSize(textBoxHeight);
  textBox->SetTextAlign(22);
  textBox->Draw();

  gPad->Update();

  TLine* line = new TLine(textBoxStart + (startBartelsRotation + bartelsRotationsToSpan) * distance, 1.0 - gPad->GetTopMargin(), textBoxStart + (startBartelsRotation + bartelsRotationsToSpan) * distance, textBoxY1);
  line->SetNDC();
  line->Draw();
}

void DrawYearIndicatorsAndFormatTitle(unsigned int totalBartelsRotationsToProcess, double yOffset = 0.0) {

  if (TPaveText* text = dynamic_cast<TPaveText*>(gPad->FindObject("title"))) {
    MoveTitleBoxToCenter();
    text->SetBorderSize(0);
  }

  CreateYearIndicators("2011", 0, 9, totalBartelsRotationsToProcess, yOffset);
  CreateYearIndicators("2012", 9, 14, totalBartelsRotationsToProcess, yOffset);
  CreateYearIndicators("2013", 23, 13, totalBartelsRotationsToProcess, yOffset);
  CreateYearIndicators("2014", 36, 14, totalBartelsRotationsToProcess, yOffset);
  CreateYearIndicators("2015", 50, 13, totalBartelsRotationsToProcess, yOffset);
  CreateYearIndicators("2016", 63, 14, totalBartelsRotationsToProcess, yOffset);
  CreateYearIndicators("2017", 77, 11, totalBartelsRotationsToProcess, yOffset);
  gPad->Update();
}

#endif // TimeTools_hh
