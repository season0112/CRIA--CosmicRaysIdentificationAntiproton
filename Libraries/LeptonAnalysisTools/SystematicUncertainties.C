#include "SystematicUncertainties.hh"

#include <cassert>
#include <cmath>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TSpline.h>

#include "AnalysisSettings.hh"
#include "FluxTools.hh"
#include "GraphTools.hh"
#include "PositronRatioTools.hh"
#include "ResolutionModels.hh"
#include "SmoothingTools.hh"
#include "Utilities.hh"

#define INFO_OUT_TAG "SystematicUncertainties"
#include "debugging.hh"

int RootFindBestMatchingPoint(TGraphErrors *grIn, double X, double Delta) {
  for (int i=0; i<grIn->GetN(); i++) {
    if (std::abs(X-grIn->GetX()[i])<Delta) return i;
  }
  return -1;
}

TGraph* AdditionalLowEnergySystematicUncertaintyGraph(TGraph* gr_SystError, TGraphAsymmErrors* grEleE3Rwth) {

  TFile* referenceFile = TFile::Open("$LEPTONANALYSIS/ReferenceFiles/LeptonFlux_MIT_20140624_NewINFNBinning.root", "OPEN");
  assert(referenceFile);

  TGraphErrors* grEleE3MIT = dynamic_cast<TGraphErrors*>(referenceFile->Get("g_ElectronFlux"));
  assert(grEleE3MIT);

  TGraphErrors* grEleE3MIT_ErrSyst = dynamic_cast<TGraphErrors*>(referenceFile->Get("g_ElectronFlux_err_sys"));
  assert(grEleE3MIT_ErrSyst);

  TGraph* gr_SystError_LowEnergy = new TGraph(grEleE3Rwth->GetN());

  for (int i=0; i<grEleE3Rwth->GetN(); i++) {
    double  X     = grEleE3Rwth->GetX()[i];
    double  Y1    = grEleE3Rwth->GetY()[i];
    double  ey    = Y1*gr_SystError->GetY()[i];
    int     ip    = RootFindBestMatchingPoint(grEleE3MIT,X,1E-4);
    if (ip < 0)
      FATAL_OUT << "No matching point found. " << std::endl;
    double  Y2     = grEleE3MIT->GetY()[ip];

    double  ey_Mit = Y2*grEleE3MIT_ErrSyst->GetY()[ip];
    double  systErr = 0.0;
    if (ey_Mit>ey && X<10.0)
      systErr = std::sqrt(std::pow(ey_Mit,2)-std::pow(ey,2));
    gr_SystError_LowEnergy->SetPoint(i,X,systErr/Y1);
  }

  return gr_SystError_LowEnergy;
}

void AddAdditionalRelativeSystematicUncertaintyToGraph(TGraphAsymmErrors* graph, TGraph* additionalSystematicUncertainty) {

  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);

    double xSystematic, ySystematic;
    additionalSystematicUncertainty->GetPoint(i, xSystematic, ySystematic);

    assert(std::abs(x - xSystematic) < 1e-5);
    assert(std::abs(graph->GetErrorYlow(i) - graph->GetErrorYhigh(i)) < 1e-5);
    assert(!std::isinf(ySystematic));
    assert(!std::isnan(ySystematic));

    double eySystematic = ySystematic * y;
    double ey = std::sqrt(std::pow(graph->GetErrorY(i), 2) + std::pow(eySystematic, 2));
    assert(!std::isinf(ey));
    assert(!std::isnan(ey));

    graph->SetPointError(i, 0, 0, ey, ey);
  }
}

void DumpTableEntry(unsigned int bin, double x, double y, double ey) {

  INFO_OUT << "bin "  <<  std::setw(2) << std::right << bin          << std::fixed << std::setprecision(2)
           << " ("    <<  std::setw(7) << std::right << x            << " GeV) "
           << "y="    << std::setw(16) << std::right << y            << " "
           << "+/-"   << std::setw(14) << std::right << ey           << " "
           << "---> " <<  std::setw(8) << std::right << ey / y * 100 << " %" << std::endl;
}

void DumpTableEntryInvalid(unsigned int bin, double x, double y, double ey, const std::string& problem) {

  WARN_OUT << "bin "  <<  std::setw(2) << std::right << bin          << std::fixed << std::setprecision(2)
           << " ("    <<  std::setw(7) << std::right << x            << " GeV) "
           << "y="    << std::setw(16) << std::right << y            << " "
           << "+/-"   << std::setw(14) << std::right << ey           << " "
           << "---> " << problem       << std::endl;
}

std::string SpeciesTitleForType(SpeciesType type) {

  switch (type) {
  case ElectronFlux:
    return "Electron flux";
  case PositronFlux:
    return "Positron flux";
  case AllElectronFlux:
    return "All-electron flux";
  case PositronFraction:
    return "Positron fraction";
   case PositronElectronRatio:
    return "Positron/electron ratio";
  }

  assert(false);
  return "";
}

std::string SpeciesNameForType(SpeciesType type) {

  switch (type) {
  case ElectronFlux:
    return "ElectronFlux";
  case PositronFlux:
    return "PositronFlux";
  case AllElectronFlux:
    return "AllElectronFlux";
  case PositronFraction:
    return "PositronFraction";
   case PositronElectronRatio:
    return "PositronElectronRatio";
  }

  assert(false);
  return "";
}

template<typename T>
TH1D* HistogramFromGraph(T* graph, const Binning::Definition& binning) {

  TH1D* histogram = Make<TH1D>(Form("%sHistogram", graph->GetName()), graph->GetTitle(), binning);
  histogram->SetLineColor(graph->GetMarkerColor());
  for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
    histogram->SetBinContent(i + 1, y ? (graph->GetY()[i]) : -1.0);
  }
  return histogram;
}

TH1D* DrawUncertainty(SpeciesType speciesType, const Binning::Definition& energyBinning, TGraph* uncertaintyGraph, bool same, double yLower, double yUpper, double xLow, double xHigh) {

  TGraph* uncertaintyGraphClone = dynamic_cast<TGraph*>(uncertaintyGraph->Clone());

  if (speciesType == ElectronFlux) {
    if (gAMSDataPeriod == gAMSDataPeriodPass4)
      RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 700.0);
    else
      RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 1000.0);
  } else if (speciesType == PositronFlux || speciesType == PositronFraction || speciesType == PositronElectronRatio) {
    if (gAMSDataPeriod == gAMSDataPeriodPass4)
      RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 600.0);
    else
      RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 1000.0);
  } else {
    assert(speciesType == AllElectronFlux);
    if (gAMSDataPeriod == gAMSDataPeriodPass4)
      RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 1000.0);
    else
      RemoveGraphPointsAboveEnergy(uncertaintyGraphClone, 1500.0);
  }

  TH1D* histogram = HistogramFromGraph(uncertaintyGraphClone, energyBinning);
  histogram->SetLineWidth(4);
  histogram->GetXaxis()->SetRangeUser(xLow, xHigh);
  histogram->GetXaxis()->SetMoreLogLabels(kTRUE);
  histogram->GetXaxis()->SetNoExponent(kTRUE);
  histogram->GetXaxis()->SetTitle("Energy / GeV");
  histogram->GetYaxis()->SetTitle("Relative uncertainty");
  histogram->GetYaxis()->SetRangeUser(yLower, yUpper);
  histogram->SetTitle("");
  histogram->Draw(same ? "][.SAME" : "][");
  return histogram;
}

TH1D* DrawUncertainty(SpeciesType speciesType, TGraph* uncertaintyGraph, bool same, double yLower, double yUpper) {

  const Binning::Definition& energyBinning = Binning::Predefined::AbsoluteEnergyBinning();
  return DrawUncertainty(speciesType, energyBinning, uncertaintyGraph, same, yLower, yUpper, gAMSStartShowEnergy, gAMSStopShowEnergy);
}

TH1D* DrawUncertaintyInBartelsRotation(SpeciesType speciesType, TGraph* uncertaintyGraph, bool same, double yLower, double yUpper) {

  Binning::Definition energyBinning = Binning::Tools::SubRange(Binning::Predefined::AbsoluteEnergyBinning(),
                                                               gAMSTimeDependentFluxesMinEnergy, gAMSTimeDependentFluxesMaxEnergy, true);

  return DrawUncertainty(speciesType, energyBinning, uncertaintyGraph, same, yLower, yUpper, gAMSTimeDependentFluxesMinEnergy, gAMSTimeDependentFluxesMaxEnergy);
}

void DrawSingleUncertainty(TGraph* uncertainty, SpeciesType speciesType, const std::string& systName, const std::string& systTitle, std::string fileSuffix, double yLower, double yUpper) {

  std::string speciesTitle = SpeciesTitleForType(speciesType);
  std::string speciesName = SpeciesNameForType(speciesType);

  TCanvas* canvas = new TCanvas(Form("canvasUncertainty_%s_%s", speciesName.c_str(), systName.c_str()),
                                Form("%s) Single uncertainty '%s'", speciesTitle.c_str(), systTitle.c_str()));
  canvas->cd();
  gPad->SetLeftMargin(0.13);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.14);
  gPad->SetLogx();
  gPad->SetLogy();

  TH1D* uncertaintyHistogram = nullptr;
  if (!fileSuffix.empty())
    uncertaintyHistogram = DrawUncertaintyInBartelsRotation(speciesType, uncertainty, false, yLower, yUpper);
  else
    uncertaintyHistogram = DrawUncertainty(speciesType, uncertainty, false, yLower, yUpper);

  TLegend* legend = new TLegend(gPad->GetLeftMargin() + 0.03, 0.87, 1.0 - gPad->GetRightMargin() - 0.15, 1.0 - gPad->GetTopMargin() - 0.05, nullptr, "brNDC");
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->AddEntry(uncertaintyHistogram, systTitle.c_str(), "l");

  legend->Draw();
  if (!fileSuffix.empty()) {
    canvas->SaveAs(Form("%s_%s.png", canvas->GetName(), fileSuffix.c_str()));
    canvas->SaveAs(Form("%s_%s.pdf", canvas->GetName(), fileSuffix.c_str()));
  } else {
    canvas->SaveAs(Form("%s.png", canvas->GetName()));
    canvas->SaveAs(Form("%s.pdf", canvas->GetName()));
  }
}

TGraph* ExtractTotalRelativeUncertaintyFromGraph(TGraphAsymmErrors* graph) {

  TGraph* uncertaintyGraph = new TGraph(graph->GetN());
  for (int point = 0; point < graph->GetN(); ++point) {
    double x, y;
    graph->GetPoint(point, x, y);

    assert(std::abs(graph->GetErrorYlow(point) - graph->GetErrorYhigh(point)) < 1e-5);
    double ey = graph->GetErrorY(point);
    uncertaintyGraph->SetPoint(point, x, y ? ey / y : 0.0);
  }

  return uncertaintyGraph;
}

TGraph* CombineRelativeUncertainties(TGraphAsymmErrors* graph, const std::vector<TGraph*>& uncertainties) {

  assert(uncertainties.size() > 0);
  TGraph* uncertaintyGraph = new TGraph(uncertainties[0]->GetN());
  for (int point = 0; point < uncertainties[0]->GetN(); ++point) {
    double x = 0.0;
    double ySquareSum = 0.0;
    for (unsigned int i = 0; i < uncertainties.size(); ++i) {
      double xNew, y;
      uncertainties[i]->GetPoint(point, xNew, y);
      ySquareSum += std::pow(y, 2);

      if (i > 0)
        assert(std::abs(x - xNew) < 1e-5);
      else
        x = xNew;
    }

    // If the graph value is 0.0 (invalid graph point), set all uncertainties also to 0.0, by convention.
    double ey = std::sqrt(ySquareSum);
    if (graph->GetY()[point] == 0.0)
      ey = 0.0;

    uncertaintyGraph->SetPoint(point, x, ey);
  }

  return uncertaintyGraph;
}

void VerifyRelativeUncertaintyGraphIsNull(TGraph* graph, const std::string& description) {

  for (int i = 0; i < graph->GetN(); ++i) {
    double x1, y1;
    graph->GetPoint(i, x1, y1);

    // Warn about any inconsistencies.
    if (y1 != 0.0)
      FATAL_OUT << "bin=" << i << " x=" << x1 << " (GeV) -> y should be null, but is " << y1 << " in " << description << "." << std::endl;
  }
}

void VerifyRelativeUncertaintyGraphsAreEqual(TGraph* graph1, TGraph* graph2, const std::string& description, double tolerance) {

  assert(graph1->GetN() == graph2->GetN());
  for (int i = 0; i < graph1->GetN(); ++i) {
    double x1, y1;
    graph1->GetPoint(i, x1, y1);

    double x2, y2;
    graph2->GetPoint(i, x2, y2);
    assert(std::abs(x1 - x2) < 1e-5);

    // Warn about any inconsistencies.
    if (std::abs(y1 - y2) > tolerance)
      FATAL_OUT << "bin=" << i << " x=" << x1 << " (GeV) -> Deviation of " << std::abs(y1 - y2) * 1000.0 << " per-mille in " << description << "." << std::endl;
  }
}

void MultiplyMigrationMatrixWithCounts(TH2D* migrationMatrix, TGraphErrors* countsBeforeMigration, TGraphErrors* countsAfterMigration, TF1* measuringTimeFunction) {

  const auto& analysisBinning = Binning::Predefined::AbsoluteEnergyBinning();
  assert(countsBeforeMigration->GetN() == int(analysisBinning.NumberOfBins()));

  for (unsigned int reconstructedEnergyBin = 1; reconstructedEnergyBin <= analysisBinning.NumberOfBins(); ++reconstructedEnergyBin) {
    double reconstructedEnergyLow = analysisBinning.LowEdge(reconstructedEnergyBin);
    double reconstructedEnergyUp = analysisBinning.UpEdge(reconstructedEnergyBin);
    double reconstructedEnergyCenter = analysisBinning.Value(reconstructedEnergyBin);
    double reconstructedEnergy = Modelling::LaffertyWyatt(reconstructedEnergyLow, reconstructedEnergyUp, gAMSLaffertyWyattSpectralIndex);

    double xMigration = migrationMatrix->GetXaxis()->GetBinCenter(reconstructedEnergyBin);
    assert(std::abs(xMigration - reconstructedEnergyCenter) < 1e-5);

    double xCounts = countsBeforeMigration->GetX()[reconstructedEnergyBin - 1];
    assert(std::abs(xCounts - reconstructedEnergy) < 1e-5);

    double migratedCounts = 0.0;
    for (unsigned int trueEnergyBin = 1; trueEnergyBin <= analysisBinning.NumberOfBins(); ++trueEnergyBin)
      migratedCounts += countsBeforeMigration->GetY()[trueEnergyBin - 1] * migrationMatrix->GetBinContent(trueEnergyBin, reconstructedEnergyBin);

    // Take geomagnetic cut-off into account.
    migratedCounts *= measuringTimeFunction->Eval(reconstructedEnergy < 0.5 ? 0.5 : reconstructedEnergy); /* normalized to 1 above geomagnetic cut-off */

    countsAfterMigration->SetPoint(reconstructedEnergyBin - 1, reconstructedEnergy, migratedCounts);
    countsAfterMigration->SetPointError(reconstructedEnergyBin - 1, 0.0, std::sqrt(migratedCounts));
  }
}

TGraph* ConstructFluxFromEventCounts(TGraphErrors* eventCount, TGraphAsymmErrors* measuringTime, TSpline5* mcEffectiveAcceptanceFunction, TSpline5* triggerEfficiencyFunction) {

  const auto& analysisBinning = Binning::Predefined::AbsoluteEnergyBinning();
  assert(eventCount->GetN() == int(analysisBinning.NumberOfBins()));

  TGraph* fluxE3 = new TGraph;
  for (unsigned int bin = 1; bin <= analysisBinning.NumberOfBins(); ++bin) {
    double trueEnergy, eventCounts;
    eventCount->GetPoint(bin - 1, trueEnergy, eventCounts);

    double energyLow = analysisBinning.LowEdge(bin);
    double energyUp = analysisBinning.UpEdge(bin);
    double binWidth = energyUp - energyLow;

    double triggerEfficiency = triggerEfficiencyFunction->Eval(trueEnergy < 0.5 ? 0.5 : trueEnergy); // Not defined below 0.5 GeV, assume it stays constant.
    double effectiveAcceptance = (mcEffectiveAcceptanceFunction->Eval(trueEnergy) / 10000.0 /* cm^2 -> m^2 */) * triggerEfficiency;
    double denominator = binWidth * triggerEfficiency * effectiveAcceptance * measuringTime->Eval(trueEnergy);

    if (denominator == 0)
      fluxE3->SetPoint(bin - 1, trueEnergy, 0);
    else
      fluxE3->SetPoint(bin - 1, trueEnergy, eventCounts / denominator * std::pow(trueEnergy, 3));
  }

  return fluxE3;
}

TGraph* DetermineMigrationSystematicUncertainty(SpeciesType speciesType, TH2D* migrationMatrix, TF1* electronFluxModel, TF1* positronFluxModel, TF1* measuringTimeFunction, TGraphAsymmErrors* measuringTime, TSpline5* mcEffectiveAcceptanceFunction, TSpline5* triggerEfficiencyFunction) {

  std::cout << "\n\n";
  INFO_OUT << "Summary of bin-to-bin migration systematic uncertainties on " << SpeciesTitleForType(speciesType) << " :" << std::endl;
  const auto& analysisBinning = Binning::Predefined::AbsoluteEnergyBinning();

  TGraphErrors* electronEventCountBeforeMigrationWithoutCutOff = new TGraphErrors;
  TGraphErrors* electronEventCountBeforeMigration = new TGraphErrors;

  TGraphErrors* positronEventCountBeforeMigrationWithoutCutOff = new TGraphErrors;
  TGraphErrors* positronEventCountBeforeMigration = new TGraphErrors;

  double fullMeasuringTime = measuringTime->Eval(100);
  for (unsigned int bin = 1; bin <= analysisBinning.NumberOfBins(); ++bin) {
    double trueEnergyLow = analysisBinning.LowEdge(bin);
    double trueEnergyUp = analysisBinning.UpEdge(bin);

    double trueElectronFluxIntegral = electronFluxModel->Integral(trueEnergyLow, trueEnergyUp, 1e-6);
    double truePositronFluxIntegral = positronFluxModel->Integral(trueEnergyLow, trueEnergyUp, 1e-6);

    double trueEnergy = Modelling::LaffertyWyatt(trueEnergyLow, trueEnergyUp, gAMSLaffertyWyattSpectralIndex);
    double triggerEfficiency = triggerEfficiencyFunction->Eval(trueEnergy < 0.5 ? 0.5 : trueEnergy); // Not defined below 0.5 GeV, assume it stays constant.
    double weightBeforeMigration = (mcEffectiveAcceptanceFunction->Eval(trueEnergy) / 10000.0 /* cm^2 -> m^2 */) * triggerEfficiency;

    double electronCountBeforeMigrationWithoutCutOff = trueElectronFluxIntegral * weightBeforeMigration * fullMeasuringTime;
    double positronCountBeforeMigrationWithoutCutOff = truePositronFluxIntegral * weightBeforeMigration * fullMeasuringTime;
    double electronCountBeforeMigration = trueElectronFluxIntegral * weightBeforeMigration * measuringTime->Eval(trueEnergy);
    double positronCountBeforeMigration = truePositronFluxIntegral * weightBeforeMigration * measuringTime->Eval(trueEnergy);

    electronEventCountBeforeMigrationWithoutCutOff->SetPoint(bin - 1, trueEnergy, electronCountBeforeMigrationWithoutCutOff);
    electronEventCountBeforeMigrationWithoutCutOff->SetPointError(bin - 1, 0.0, std::sqrt(electronCountBeforeMigrationWithoutCutOff));

    electronEventCountBeforeMigration->SetPoint(bin - 1, trueEnergy, electronCountBeforeMigration);
    electronEventCountBeforeMigration->SetPointError(bin - 1, 0.0, std::sqrt(electronCountBeforeMigration));

    positronEventCountBeforeMigrationWithoutCutOff->SetPoint(bin - 1, trueEnergy, positronCountBeforeMigrationWithoutCutOff);
    positronEventCountBeforeMigrationWithoutCutOff->SetPointError(bin - 1, 0.0, std::sqrt(positronCountBeforeMigrationWithoutCutOff));

    positronEventCountBeforeMigration->SetPoint(bin - 1, trueEnergy, positronCountBeforeMigration);
    positronEventCountBeforeMigration->SetPointError(bin - 1, 0.0, std::sqrt(positronCountBeforeMigration));
  }

  Utilities::NormalizeHistogramXSlices(migrationMatrix);
  assert(migrationMatrix->GetNbinsX() == electronEventCountBeforeMigration->GetN());
  assert(migrationMatrix->GetNbinsY() == electronEventCountBeforeMigration->GetN());

  TGraphErrors* electronEventCountAfterMigration = new TGraphErrors;
  TGraphErrors* positronEventCountAfterMigration = new TGraphErrors;

  MultiplyMigrationMatrixWithCounts(migrationMatrix, electronEventCountBeforeMigrationWithoutCutOff, electronEventCountAfterMigration, measuringTimeFunction);
  MultiplyMigrationMatrixWithCounts(migrationMatrix, positronEventCountBeforeMigrationWithoutCutOff, positronEventCountAfterMigration, measuringTimeFunction);

  TGraphErrors* migrationRatioGraph = new TGraphErrors;
  migrationRatioGraph->SetTitle(Form("After / before migration ratio for '%s'", SpeciesTitleForType(speciesType).c_str()));

  for (unsigned int bin = 1; bin <= analysisBinning.NumberOfBins(); ++bin) {
    double yElectronEventCount = electronEventCountBeforeMigration->GetY()[bin - 1];
    double yPositronEventCount = positronEventCountBeforeMigration->GetY()[bin - 1];

    double yElectronEventCountSmeared = electronEventCountAfterMigration->GetY()[bin - 1];
    double yPositronEventCountSmeared = positronEventCountAfterMigration->GetY()[bin - 1];

    double yAllElectronEventCount = yElectronEventCount + yPositronEventCount;
    double yAllElectronEventCountSmeared = yElectronEventCountSmeared + yPositronEventCountSmeared;

    double positronRatio = 0.0;
    double positronRatioUncertainty = 0.0;
    if (speciesType == PositronFraction)
      CalculatePositronFraction(positronRatio, positronRatioUncertainty, yElectronEventCount, 0.0, yPositronEventCount, 0.0);
    else if (speciesType == PositronElectronRatio)
      CalculatePositronElectronRatio(positronRatio, positronRatioUncertainty, yElectronEventCount, 0.0, yPositronEventCount, 0.0);

    double positronRatioSmeared = 0.0;
    double positronRatioSmearedUncertainty = 0.0;
    if (speciesType == PositronFraction)
      CalculatePositronFraction(positronRatioSmeared, positronRatioSmearedUncertainty, yElectronEventCountSmeared, 0.0, yPositronEventCountSmeared, 0.0);
    else if (speciesType == PositronElectronRatio)
      CalculatePositronElectronRatio(positronRatioSmeared, positronRatioSmearedUncertainty, yElectronEventCountSmeared, 0.0, yPositronEventCountSmeared, 0.0);

    double ratio = 0.0;
    if (speciesType == ElectronFlux)
      ratio = yElectronEventCount > 0 ? yElectronEventCountSmeared / yElectronEventCount : 1.0;
    else if (speciesType == PositronFlux)
      ratio = yPositronEventCount > 0 ? yPositronEventCountSmeared / yPositronEventCount : 1.0;
    else if (speciesType == AllElectronFlux)
      ratio = yAllElectronEventCount > 0 ? yAllElectronEventCountSmeared / yAllElectronEventCount : 1.0;
    else if (speciesType == PositronFraction || speciesType == PositronElectronRatio)
      ratio = positronRatio > 0 ? positronRatioSmeared / positronRatio : 1.0;

    double xBinLow = analysisBinning.LowEdge(bin);
    double xBinUp = analysisBinning.UpEdge(bin);
    double xTrue = Modelling::LaffertyWyatt(xBinLow, xBinUp, gAMSLaffertyWyattSpectralIndex);

    double ratioUncertainty = std::abs(ratio - 1.0);
    assert(std::abs(electronEventCountBeforeMigration->GetX()[bin - 1] - xTrue) < 1e-5);
    assert(std::abs(positronEventCountBeforeMigration->GetX()[bin - 1] - xTrue) < 1e-5);
    assert(std::abs(electronEventCountAfterMigration->GetX()[bin - 1] - xTrue) < 1e-5);
    assert(std::abs(positronEventCountAfterMigration->GetX()[bin - 1] - xTrue) < 1e-5);
    migrationRatioGraph->SetPoint(bin - 1, xTrue, ratio);
    migrationRatioGraph->SetPointError(bin - 1, 0.0, ratioUncertainty);
  }

  // Draw electron flux before & after migration
  TGraph* electronFluxE3BeforeMigration = ConstructFluxFromEventCounts(electronEventCountBeforeMigration, measuringTime, mcEffectiveAcceptanceFunction, triggerEfficiencyFunction);
  TGraph* electronFluxE3AfterMigration = ConstructFluxFromEventCounts(electronEventCountAfterMigration, measuringTime, mcEffectiveAcceptanceFunction, triggerEfficiencyFunction);

  TCanvas* electronFluxMigrationCanvas = new TCanvas(Form("electronFluxMigrationCanvas_%s", SpeciesNameForType(speciesType).c_str()), Form("%s) Electron flux migration effect", SpeciesNameForType(speciesType).c_str()));
  electronFluxMigrationCanvas->cd();
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.10);
  gPad->SetGrid();
  gPad->SetLogx();
  electronFluxE3AfterMigration->SetLineColor(kRed);
  electronFluxE3AfterMigration->SetMarkerColor(kRed);
  electronFluxE3AfterMigration->SetMarkerStyle(kFullCircle);
  electronFluxE3AfterMigration->SetMarkerSize(0.8);
  electronFluxE3BeforeMigration->SetLineColor(kBlue);
  electronFluxE3BeforeMigration->SetMarkerColor(kBlue);
  electronFluxE3BeforeMigration->SetMarkerStyle(kOpenCircle);
  electronFluxE3BeforeMigration->SetMarkerSize(1.3);
  electronFluxE3BeforeMigration->GetXaxis()->SetTitle("Energy / GeV");
  electronFluxE3BeforeMigration->GetYaxis()->SetTitle("Electron flux");
  electronFluxE3BeforeMigration->GetYaxis()->SetRangeUser(1, 3e2);
  electronFluxE3BeforeMigration->GetYaxis()->SetNoExponent(kTRUE);
  electronFluxE3BeforeMigration->GetXaxis()->SetRangeUser(gAMSMinEnergy, gAMSMaxEnergy);
  electronFluxE3BeforeMigration->Draw("ap");
  electronFluxE3AfterMigration->Draw("p.same");

  TLegend* electronFluxMigrationLegend = new TLegend(0.3202811, 0.1220657, 0.8056627, 0.2597653, 0, "brNDC");
  electronFluxMigrationLegend->SetTextSize(0.05);
  electronFluxMigrationLegend->SetFillColor(kWhite);
  electronFluxMigrationLegend->AddEntry(electronFluxE3BeforeMigration, "e^{-} - before migration", "pl");
  electronFluxMigrationLegend->AddEntry(electronFluxE3AfterMigration, "e^{-} - after migration", "pl");
  electronFluxMigrationLegend->Draw();

  electronFluxMigrationCanvas->SaveAs(Form("%s.png", electronFluxMigrationCanvas->GetName()));
  electronFluxMigrationCanvas->SaveAs(Form("%s.pdf", electronFluxMigrationCanvas->GetName()));

  // Draw positron flux before & after migration
  TGraph* positronFluxE3BeforeMigration = ConstructFluxFromEventCounts(positronEventCountBeforeMigration, measuringTime, mcEffectiveAcceptanceFunction, triggerEfficiencyFunction);
  TGraph* positronFluxE3AfterMigration = ConstructFluxFromEventCounts(positronEventCountAfterMigration, measuringTime, mcEffectiveAcceptanceFunction, triggerEfficiencyFunction);

  TCanvas* positronFluxMigrationCanvas = new TCanvas(Form("positronFluxMigrationCanvas_%s", SpeciesNameForType(speciesType).c_str()), Form("%s) Positron flux migration effect", SpeciesNameForType(speciesType).c_str()));
  positronFluxMigrationCanvas->cd();
  gPad->SetLeftMargin(0.1);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.10);
  gPad->SetGrid();
  gPad->SetLogx();
  positronFluxE3AfterMigration->SetLineColor(kGreen + 2);
  positronFluxE3AfterMigration->SetMarkerColor(kGreen + 2);
  positronFluxE3AfterMigration->SetMarkerStyle(kFullCircle);
  positronFluxE3AfterMigration->SetMarkerSize(0.8);
  positronFluxE3BeforeMigration->SetLineColor(kOrange - 2);
  positronFluxE3BeforeMigration->SetMarkerColor(kOrange - 2);
  positronFluxE3BeforeMigration->SetMarkerStyle(kOpenCircle);
  positronFluxE3BeforeMigration->SetMarkerSize(1.3);
  positronFluxE3BeforeMigration->GetXaxis()->SetTitle("Energy / GeV");
  positronFluxE3BeforeMigration->GetYaxis()->SetTitle("Positron flux");
  positronFluxE3BeforeMigration->GetYaxis()->SetLabelSize(0.04);
  positronFluxE3BeforeMigration->GetYaxis()->SetTitleOffset(0.8);
  positronFluxE3BeforeMigration->GetYaxis()->SetTitleSize(0.06);
  positronFluxE3BeforeMigration->GetYaxis()->SetRangeUser(1, 25);
  positronFluxE3BeforeMigration->GetYaxis()->SetNoExponent(kTRUE);
  positronFluxE3BeforeMigration->GetXaxis()->SetRangeUser(gAMSMinEnergy, gAMSMaxEnergy);
  positronFluxE3BeforeMigration->Draw("ap");
  positronFluxE3AfterMigration->Draw("p.same");

  TLegend* positronFluxMigrationLegend = new TLegend(0.3202811, 0.1220657, 0.8056627, 0.2597653, 0, "brNDC");
  positronFluxMigrationLegend->SetTextSize(0.05);
  positronFluxMigrationLegend->SetFillColor(kWhite);
  positronFluxMigrationLegend->AddEntry(positronFluxE3BeforeMigration, "e^{+} - before migration", "pl");
  positronFluxMigrationLegend->AddEntry(positronFluxE3AfterMigration, "e^{+} - after migration", "pl");
  positronFluxMigrationLegend->Draw();

  positronFluxMigrationCanvas->SaveAs(Form("%s.png", positronFluxMigrationCanvas->GetName()));
  positronFluxMigrationCanvas->SaveAs(Form("%s.pdf", positronFluxMigrationCanvas->GetName()));

  // Draw after / before migration ratio
  TCanvas* migrationRatioGraphCanvas = new TCanvas(Form("migrationRatioGraphCanvas_%s", SpeciesNameForType(speciesType).c_str()), Form("%s) Smearing ratio", SpeciesNameForType(speciesType).c_str()));
  gStyle->SetOptTitle(1);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  migrationRatioGraphCanvas->Divide(1, 2, 1e-5, 1e-5);
  migrationRatioGraphCanvas->cd(1);
  gPad->SetLeftMargin(0.075);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(1e-5);
  gPad->SetPad(0, 0.3364486, 1, 1);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetLogy();
  positronEventCountAfterMigration->SetLineColor(kGreen + 2);
  positronEventCountAfterMigration->SetMarkerColor(kGreen + 2);
  positronEventCountAfterMigration->SetMarkerStyle(kFullCircle);
  positronEventCountAfterMigration->SetMarkerSize(0.8);
  positronEventCountBeforeMigration->SetLineColor(kOrange - 2);
  positronEventCountBeforeMigration->SetMarkerColor(kOrange - 2);
  positronEventCountBeforeMigration->SetMarkerStyle(kOpenCircle);
  positronEventCountBeforeMigration->SetMarkerSize(1.3);
  electronEventCountAfterMigration->SetLineColor(kRed);
  electronEventCountAfterMigration->SetMarkerColor(kRed);
  electronEventCountAfterMigration->SetMarkerStyle(kFullCircle);
  electronEventCountAfterMigration->SetMarkerSize(0.8);
  electronEventCountBeforeMigration->SetLineColor(kBlue);
  electronEventCountBeforeMigration->SetMarkerColor(kBlue);
  electronEventCountBeforeMigration->SetMarkerStyle(kOpenCircle);
  electronEventCountBeforeMigration->SetMarkerSize(1.3);
  electronEventCountBeforeMigration->GetYaxis()->SetTitle("Event counts");
  electronEventCountBeforeMigration->GetYaxis()->SetLabelOffset(0.003);
  electronEventCountBeforeMigration->GetYaxis()->SetLabelSize(0.04);
  electronEventCountBeforeMigration->GetYaxis()->SetTitleOffset(0.5);
  electronEventCountBeforeMigration->GetYaxis()->SetTitleSize(0.07);
  electronEventCountBeforeMigration->GetYaxis()->SetRangeUser(1, 5e7);
  electronEventCountBeforeMigration->GetXaxis()->SetRangeUser(gAMSMinEnergy, gAMSMaxEnergy);
  electronEventCountBeforeMigration->Draw("ap");
  positronEventCountBeforeMigration->Draw("p.same");
  electronEventCountAfterMigration->Draw("p.same");
  positronEventCountAfterMigration->Draw("p.same");

  TLegend* legend = new TLegend(0.2202811, 0.1220657, 0.5056627, 0.3597653, 0, "brNDC");
  legend->SetTextSize(0.05);
  legend->SetFillColor(kWhite);
  legend->AddEntry(electronEventCountBeforeMigration, "e^{-} - before migration", "pl");
  legend->AddEntry(electronEventCountAfterMigration, "e^{-} - after migration", "pl");
  legend->AddEntry(positronEventCountBeforeMigration, "e^{+} - before migration", "pl");
  legend->AddEntry(positronEventCountAfterMigration, "e^{+} - after migration", "pl");
  legend->Draw();

  migrationRatioGraphCanvas->cd(2);
  gPad->SetPad(0, 0, 1, 0.3364486);
  gPad->SetLeftMargin(0.07429719);
  gPad->SetRightMargin(0.02008032);
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.2136285);
  gPad->SetGrid();
  gPad->SetLogx();
  gPad->SetTickx();
  migrationRatioGraph->SetFillColor(kOrange + 7);
  migrationRatioGraph->SetFillStyle(3001);
  migrationRatioGraph->SetMarkerColor(kRed);
  migrationRatioGraph->SetMarkerStyle(20);
  migrationRatioGraph->SetMarkerSize(0.8);

  migrationRatioGraph->GetXaxis()->SetLimits(gAMSMinEnergy, gAMSMaxEnergy);
  migrationRatioGraph->GetXaxis()->SetLabelOffset(0.011);
  migrationRatioGraph->GetXaxis()->SetLabelSize(0.1);
  migrationRatioGraph->GetXaxis()->SetTitleOffset(1.31);
  migrationRatioGraph->GetXaxis()->SetTitleSize(0.08);
  migrationRatioGraph->GetYaxis()->SetLabelOffset(0.003);
  migrationRatioGraph->GetYaxis()->SetLabelSize(0.1);
  migrationRatioGraph->GetYaxis()->SetTitleOffset(0.31);
  migrationRatioGraph->GetYaxis()->SetTitleSize(0.11);
  migrationRatioGraph->GetYaxis()->SetTitle("Ratio");
  migrationRatioGraph->GetYaxis()->SetRangeUser(0.8, 1.05);

  migrationRatioGraph->Draw("AP");
  migrationRatioGraph->Draw("CE3.SAME");

  MoveTitleBoxToCenter();

  TPaveText* titleBox = dynamic_cast<TPaveText*>(gPad->FindObject("title"));
  titleBox->SetY1NDC(0.86);
  titleBox->SetBorderSize(0);

  migrationRatioGraphCanvas->SaveAs(Form("%s.png", migrationRatioGraphCanvas->GetName()));
  migrationRatioGraphCanvas->SaveAs(Form("%s.pdf", migrationRatioGraphCanvas->GetName()));

  TGraph* uncertaintyGraph = new TGraph;
  for (unsigned int i = 1; i <= analysisBinning.NumberOfBins(); ++i) {
    double x = electronEventCountBeforeMigration->GetX()[i - 1];
    double xTrue = Modelling::LaffertyWyatt(analysisBinning.LowEdge(i), analysisBinning.UpEdge(i), gAMSLaffertyWyattSpectralIndex);
    assert(std::abs(x - xTrue) < 1e-5);

    double y = 1.0;
    assert(y >= 0.0);

    double ey = std::abs(migrationRatioGraph->GetY()[i - 1] - 1.0);
    if (y == 0.0) {
      ey = 0.0;
      y = 1.0;
      DumpTableEntryInvalid(i - 1, xTrue, 0.0, ey, "Entries are null!");
    } else
      DumpTableEntry(i - 1, xTrue, y, ey);

    uncertaintyGraph->SetPoint(i - 1, xTrue, ey / y);
  }

  return SmoothLogarithmicGraph(uncertaintyGraph);
}
