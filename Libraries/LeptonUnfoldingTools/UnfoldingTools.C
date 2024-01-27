#include "UnfoldingTools.hh"

#include <cmath>

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <TH2.h>
#include <TLegend.h>
#include <TROOT.h>

#include <BinningTools.hh>
#include <BrokenPowerLawModel.hh>
#include <ForwardFoldingDataset.hh>
#include <FluxTools.hh>
#include <HistogramDataset.hh>
#include <PredefinedBinnings.hh>
#include <Utilities.hh>

#include "AnalysisSettings.hh"
#include "BayesUnfoldingWithCutoff.hh"
#include "MatrixTools.hh"
#include "ModelFits.hh"

#define INFO_OUT_TAG "UnfoldingTools"
#include <debugging.hh>

static const bool sVerboseUnfoldingLog = false;
static const bool sDebugUnfoldingAndForwardFolding = true;

TH1D* ConvertToHistogram(const Binning::Definition& energyBinning, TGraphAsymmErrors* flux) {

  TH1D* histogram = Make<TH1D>(Form("%s_histogram", flux->GetName()), "", energyBinning);
  for (unsigned int energyBin = 0; energyBin < energyBinning.NumberOfBins(); ++energyBin) {
    double xTrue = Modelling::LaffertyWyatt(energyBinning.LowEdge(energyBin + 1),
                                            energyBinning.UpEdge(energyBin + 1),
                                            gAMSLaffertyWyattSpectralIndex);

    double yFlux = flux->GetY()[energyBin];
    assert(std::abs(flux->GetX()[energyBin] - xTrue) < 1e-5);
    histogram->SetBinContent(energyBin + 1, yFlux);
    histogram->SetBinError(energyBin + 1, flux->GetErrorY(energyBin));
  }

  return histogram;
}

TGraphAsymmErrors* ConstructUnfoldedAllElectronEventCount(TGraphAsymmErrors* electronEventCountUnfolded, TGraphAsymmErrors* positronEventCountUnfolded) {

  TGraphAsymmErrors* allElectronEventCountUnfolded = new TGraphAsymmErrors;

  for (int i = 0; i < electronEventCountUnfolded->GetN(); ++i) {
    double xElectron, yElectron;
    electronEventCountUnfolded->GetPoint(i, xElectron, yElectron);

    double eyElectronLow = electronEventCountUnfolded->GetEYlow()[i];
    double eyElectronHigh = electronEventCountUnfolded->GetEYhigh()[i];

    double xPositron, yPositron;
    positronEventCountUnfolded->GetPoint(i, xPositron, yPositron);

    double eyPositronLow = positronEventCountUnfolded->GetEYlow()[i];
    double eyPositronHigh = positronEventCountUnfolded->GetEYhigh()[i];

    assert(std::abs(xElectron - xPositron) < 1e-5);

    double yAllElectron = yElectron + yPositron;
    double eyAllElectronLow = std::sqrt(std::pow(eyElectronLow, 2) + std::pow(eyPositronLow, 2));
    double eyAllElectronHigh = std::sqrt(std::pow(eyElectronHigh, 2) + std::pow(eyPositronHigh, 2));

    allElectronEventCountUnfolded->SetPoint(i, xElectron, yAllElectron);
    allElectronEventCountUnfolded->SetPointError(i, 0.0, 0.0, eyAllElectronLow, eyAllElectronHigh);
  }

  return allElectronEventCountUnfolded;
}

TGraphAsymmErrors* UnfoldEventCounts(const std::string& inputFileSuffix, TGraphAsymmErrors* eventCounts, TGraphAsymmErrors* measuringTime, TH2* migrationMatrix, const std::string& speciesSymbol, bool isTimeAveragedFlux, TGraphAsymmErrors* effectiveAcceptance, TGraphAsymmErrors* triggerEfficiency, TGraphAsymmErrors* ecalBDTEfficiency) {

  const Binning::Definition& binning = Binning::Predefined::AbsoluteEnergyBinning();

  TH1* hEventCounts = ConvertToHistogram(binning, eventCounts);
  TH1* hMeasuringTime = ConvertToHistogram(binning, measuringTime);
  TH1* hAcceptance = effectiveAcceptance ? ConvertToHistogram(binning, effectiveAcceptance) : nullptr;
  TH1* hTriggerEfficiency = triggerEfficiency ? ConvertToHistogram(binning, triggerEfficiency) : nullptr;
  TH1* hEcalBDTEfficiency = ecalBDTEfficiency ? ConvertToHistogram(binning, ecalBDTEfficiency) : nullptr;

  if (hAcceptance) {
    hAcceptance->Scale(1.0 / 10000.0 /* cm^2 -> m^2 */);

    for (int bin = 1; bin <= hAcceptance->GetNbinsX(); ++bin) {
      double ecalBDTEfficiency = hEcalBDTEfficiency->GetBinContent(bin);
      assert(ecalBDTEfficiency >= 0.0);
      hAcceptance->SetBinContent(bin, hAcceptance->GetBinContent(bin) * ecalBDTEfficiency);
      hAcceptance->SetBinError(bin, hAcceptance->GetBinError(bin) * ecalBDTEfficiency);
    }
  }

  unsigned int binningOffset = 1;
  Binning::Definition reducedBinning = binning;
  if (isTimeAveragedFlux)
    reducedBinning = Binning::Tools::SubRange(binning, 0.5, 1000, true);
  else
    reducedBinning = Binning::Tools::SubRange(binning, gAMSTimeDependentFluxesMinEnergy, gAMSTimeDependentFluxesMaxEnergy, true);

  hEventCounts = MatrixTools::MakeStrippedHistogram1D(hEventCounts, reducedBinning, MatrixTools::UnderOverflowHandling::ClosestBin);
  hMeasuringTime = MatrixTools::MakeStrippedHistogram1D(hMeasuringTime, reducedBinning, MatrixTools::UnderOverflowHandling::ClosestBin);
  if (hAcceptance)
    hAcceptance = MatrixTools::MakeStrippedHistogram1D(hAcceptance, reducedBinning, MatrixTools::UnderOverflowHandling::ClosestBin);
  if (hTriggerEfficiency)
    hTriggerEfficiency = MatrixTools::MakeStrippedHistogram1D(hTriggerEfficiency, reducedBinning, MatrixTools::UnderOverflowHandling::ClosestBin);
  if (!isTimeAveragedFlux) {
    // migrationMatrix is already in the range of 0.5 - 1000 GeV for time-averaged flux.
    migrationMatrix = MatrixTools::MakeStrippedHistogram2D(migrationMatrix, reducedBinning, true /* setUnderOverflows */);
    binningOffset++;
  }

  const int bayesIterations = 4;
  const int verbosity = sVerboseUnfoldingLog ? -1 : 0;
  BayesUnfoldingWithCutoff unfolding(hEventCounts, hMeasuringTime, hAcceptance, hTriggerEfficiency, migrationMatrix, bayesIterations, verbosity);
  TH1* hUnfoldedEventCounts = unfolding.UnfoldedEventCounts();

  if (isTimeAveragedFlux) {
    TGraphAsymmErrors* grEventCounts = new TGraphAsymmErrors;
    TGraphAsymmErrors* grEventCountsWithoutCutoff = new TGraphAsymmErrors;
    TGraphAsymmErrors* grUnfoldedEventCounts = new TGraphAsymmErrors;
    TGraphAsymmErrors* grUnfoldedEventCountsToEventCountsRatio = new TGraphAsymmErrors;
    TGraphAsymmErrors* grUnfoldedEventCountsWithoutCutoff = new TGraphAsymmErrors;
    for (int bin = 1; bin <= hEventCounts->GetNbinsX(); ++bin) {
      double xTrue = Modelling::LaffertyWyatt(reducedBinning.LowEdge(bin), reducedBinning.UpEdge(bin), gAMSLaffertyWyattSpectralIndex);

      grEventCounts->SetPoint(bin - 1, xTrue, hEventCounts->GetBinContent(bin));
      grEventCounts->SetPointError(bin - 1, 0.0, 0.0, hEventCounts->GetBinError(bin), hEventCounts->GetBinError(bin));

      double cutoffWeight = hMeasuringTime->GetBinContent(bin) / hMeasuringTime->GetBinContent(hMeasuringTime->GetNbinsX());
      grEventCountsWithoutCutoff->SetPoint(bin - 1, xTrue, hEventCounts->GetBinContent(bin) / cutoffWeight);
      grEventCountsWithoutCutoff->SetPointError(bin - 1, 0.0, 0.0, hEventCounts->GetBinError(bin) / cutoffWeight, hEventCounts->GetBinError(bin) / cutoffWeight);

      grUnfoldedEventCounts->SetPoint(bin - 1, xTrue, hUnfoldedEventCounts->GetBinContent(bin));
      grUnfoldedEventCounts->SetPointError(bin - 1, 0.0, 0.0, hUnfoldedEventCounts->GetBinError(bin), hUnfoldedEventCounts->GetBinError(bin));

      double unfoldingRatio = hUnfoldedEventCounts->GetBinContent(bin) / hEventCounts->GetBinContent(bin);
      double unfoldingRatioError = std::sqrt(std::pow(hUnfoldedEventCounts->GetBinError(bin) / hUnfoldedEventCounts->GetBinContent(bin), 2) +
                                             std::pow(hEventCounts->GetBinError(bin) / hEventCounts->GetBinContent(bin), 2)) * unfoldingRatio;
      grUnfoldedEventCountsToEventCountsRatio->SetPoint(bin - 1, xTrue, unfoldingRatio);
      grUnfoldedEventCountsToEventCountsRatio->SetPointError(bin - 1, 0.0, 0.0, unfoldingRatioError, unfoldingRatioError);

      grUnfoldedEventCountsWithoutCutoff->SetPoint(bin - 1, xTrue, hUnfoldedEventCounts->GetBinContent(bin) / cutoffWeight);
      grUnfoldedEventCountsWithoutCutoff->SetPointError(bin - 1, 0.0, 0.0, hUnfoldedEventCounts->GetBinError(bin) / cutoffWeight, hUnfoldedEventCounts->GetBinError(bin) / cutoffWeight);
    }

    // Plot: event counts vs. event counts without cutoff
    TCanvas* canvasEventCountsWithoutCutoff = new TCanvas("canvasEventCountsWithoutCutoff", Form("%s) event counts vs. event counts without cutoff", speciesSymbol.c_str()));
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetLogx();
    gPad->SetLogy();

    grEventCountsWithoutCutoff->SetMarkerStyle(kOpenCircle);
    grEventCountsWithoutCutoff->SetMarkerSize(2.0);
    grEventCountsWithoutCutoff->SetMarkerColor(kBlue);
    grEventCountsWithoutCutoff->GetXaxis()->SetLimits(gAMSStartShowEnergy, gAMSStopShowEnergy);
    grEventCountsWithoutCutoff->GetXaxis()->SetTitle("ECAL Energy / GeV");
    grEventCountsWithoutCutoff->GetXaxis()->SetNoExponent(kTRUE);
    grEventCountsWithoutCutoff->GetXaxis()->SetMoreLogLabels();
    grEventCountsWithoutCutoff->GetYaxis()->SetRangeUser(3, 5e7);
    grEventCountsWithoutCutoff->GetYaxis()->SetTitle("Event counts");
    grEventCountsWithoutCutoff->Draw("ap");

    grEventCounts->SetMarkerStyle(kFullCircle);
    grEventCounts->SetMarkerSize(2.0);
    grEventCounts->SetMarkerColor(kRed);
    grEventCounts->Draw("p.same");
    canvasEventCountsWithoutCutoff->SaveAs(Form("%s_%s_%s.png", canvasEventCountsWithoutCutoff->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    canvasEventCountsWithoutCutoff->SaveAs(Form("%s_%s_%s.pdf", canvasEventCountsWithoutCutoff->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(canvasEventCountsWithoutCutoff);

    // Plot: event counts without cutoff vs. unfolded event counts without cutoff
    TCanvas* canvasUnfoldedEventCountsWithoutCutoff = new TCanvas("canvasUnfoldedEventCountsWithoutCutoff", Form("%s) event counts without cutoff vs. unfolded event counts without cutoff", speciesSymbol.c_str()));
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetLogx();

    auto* grEventCountsWithoutCutoffClone = static_cast<TGraphAsymmErrors*>(grEventCountsWithoutCutoff->DrawClone("ap"));
    if (speciesSymbol == "e-")
      grEventCountsWithoutCutoffClone->GetYaxis()->SetRangeUser(1e5, 2.5e7);
    else if (speciesSymbol == "e+") {
      TGaxis::SetMaxDigits(2);
      grEventCountsWithoutCutoffClone->GetYaxis()->SetRangeUser(1e4, 2.7e6);
    }
    grEventCountsWithoutCutoffClone->GetXaxis()->SetLimits(gAMSStartShowEnergy, 10);
    grUnfoldedEventCountsWithoutCutoff->SetMarkerStyle(kFullCircle);
    grUnfoldedEventCountsWithoutCutoff->SetMarkerSize(2.0);
    grUnfoldedEventCountsWithoutCutoff->SetMarkerColor(kRed);
    grUnfoldedEventCountsWithoutCutoff->Draw("p.same");
    canvasUnfoldedEventCountsWithoutCutoff->SaveAs(Form("%s_%s_%s.png", canvasUnfoldedEventCountsWithoutCutoff->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    canvasUnfoldedEventCountsWithoutCutoff->SaveAs(Form("%s_%s_%s.pdf", canvasUnfoldedEventCountsWithoutCutoff->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(canvasUnfoldedEventCountsWithoutCutoff);

    // Plot: event counts vs. unfolded event counts
    TCanvas* canvasUnfoldedEventCounts = new TCanvas("canvasUnfoldedEventCounts", Form("%s) event counts vs. unfolded event counts", speciesSymbol.c_str()));
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetLogx();
    gPad->SetLogy();

    auto* grEventCountsClone = static_cast<TGraphAsymmErrors*>(grEventCounts->Clone());
    grEventCountsClone->SetMarkerStyle(kOpenCircle);
    grEventCountsClone->SetMarkerSize(2.0);
    grEventCountsClone->SetMarkerColor(kBlue);
    grEventCountsClone->GetXaxis()->SetLimits(gAMSStartShowEnergy, gAMSStopShowEnergy);
    grEventCountsClone->GetXaxis()->SetTitle("Energy / GeV");
    grEventCountsClone->GetXaxis()->SetNoExponent(kTRUE);
    grEventCountsClone->GetXaxis()->SetMoreLogLabels();
    grEventCountsClone->GetYaxis()->SetTitle("Event counts");
    grEventCountsClone->GetYaxis()->SetRangeUser(3, 5e6);
    grEventCountsClone->Draw("ap");

    grUnfoldedEventCounts->SetMarkerStyle(kFullCircle);
    grUnfoldedEventCounts->SetMarkerSize(2.0);
    grUnfoldedEventCounts->SetMarkerColor(kRed);
    grUnfoldedEventCounts->Draw("p.same");
    canvasUnfoldedEventCounts->SaveAs(Form("%s_%s_%s.png", canvasUnfoldedEventCounts->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    canvasUnfoldedEventCounts->SaveAs(Form("%s_%s_%s.pdf", canvasUnfoldedEventCounts->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(canvasUnfoldedEventCounts);

    // Plot: unfolded event counts / event counts
    TCanvas* canvasUnfoldedEventCountsRatio = new TCanvas("canvasUnfoldedEventCountsRatio", Form("%s) ratio unfolded event counts / event counts", speciesSymbol.c_str()));
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetLogx();

    grUnfoldedEventCountsToEventCountsRatio->SetMarkerStyle(kFullCircle);
    grUnfoldedEventCountsToEventCountsRatio->SetMarkerSize(2.0);
    grUnfoldedEventCountsToEventCountsRatio->SetMarkerColor(kRed);
    grUnfoldedEventCountsToEventCountsRatio->GetXaxis()->SetLimits(gAMSStartShowEnergy, gAMSStopShowEnergy);
    grUnfoldedEventCountsToEventCountsRatio->GetXaxis()->SetTitle("Energy / GeV");
    grUnfoldedEventCountsToEventCountsRatio->GetXaxis()->SetNoExponent(kTRUE);
    grUnfoldedEventCountsToEventCountsRatio->GetXaxis()->SetMoreLogLabels();
    grUnfoldedEventCountsToEventCountsRatio->GetYaxis()->SetTitle("Event counts");
    grUnfoldedEventCountsToEventCountsRatio->GetYaxis()->SetRangeUser(0.62, 1.38);
    grUnfoldedEventCountsToEventCountsRatio->Draw("ap");

    TLine* oneIndicatorLine = new TLine(gAMSStartShowEnergy, 1, gAMSStopShowEnergy, 1);
    oneIndicatorLine->SetLineStyle(7);
    oneIndicatorLine->SetLineWidth(4);
    oneIndicatorLine->Draw();

    canvasUnfoldedEventCountsRatio->SaveAs(Form("%s_%s_%s.png", canvasUnfoldedEventCountsRatio->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    canvasUnfoldedEventCountsRatio->SaveAs(Form("%s_%s_%s.pdf", canvasUnfoldedEventCountsRatio->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(canvasUnfoldedEventCountsRatio);
  }

  if (sDebugUnfoldingAndForwardFolding && hAcceptance) {
    double Emin = isTimeAveragedFlux ? 1.0 : reducedBinning.Min();
    double Emax = isTimeAveragedFlux ? 1000.0 : reducedBinning.Max();

    TH1* hFlux = static_cast<TH1*>(hEventCounts->Clone());
    TH1* hUnfoldedFlux = static_cast<TH1*>(hEventCounts->Clone());

    for (int bin = 0; bin <= hFlux->GetNbinsX() + 1; ++bin) {
      double rawDenominator = hEventCounts->GetBinWidth(bin) * hAcceptance->GetBinContent(bin) * hTriggerEfficiency->GetBinContent(bin) * hMeasuringTime->GetBinContent(bin);
      double rawFlux = rawDenominator > 0 ? hEventCounts->GetBinContent(bin) / rawDenominator : 0.0;
      double rawFluxErr = rawDenominator > 0 ? hEventCounts->GetBinError(bin) / rawDenominator : 0.0;
      hFlux->SetBinContent(bin, rawFlux);
      hFlux->SetBinError(bin, rawFluxErr);
    }

    for (int bin = 0; bin <= hUnfoldedFlux->GetNbinsX() + 1; ++bin) {
      double unfoldedDenominator = hUnfoldedEventCounts->GetBinWidth(bin) * hAcceptance->GetBinContent(bin) * hTriggerEfficiency->GetBinContent(bin) * hMeasuringTime->GetBinContent(bin);
      double unfoldedFlux = unfoldedDenominator > 0 ? hUnfoldedEventCounts->GetBinContent(bin) / unfoldedDenominator : 0.0;
      double unfoldedFluxErr = unfoldedDenominator > 0 ? hUnfoldedEventCounts->GetBinError(bin) / unfoldedDenominator : 0.0;
      hUnfoldedFlux->SetBinContent(bin, unfoldedFlux);
      hUnfoldedFlux->SetBinError(bin, unfoldedFluxErr);
    }

    // add systematic uncertainty
    const double relSystErr = 0.01;
    for (int ibin = 1; ibin <= hFlux->GetNbinsX(); ++ibin)
      hFlux->SetBinError(ibin, Utilities::QuadraticSum(hFlux->GetBinError(ibin), relSystErr*hFlux->GetBinContent(ibin)));
    for (int ibin = 1; ibin <= hEventCounts->GetNbinsX(); ++ibin)
      hEventCounts->SetBinError(ibin, Utilities::QuadraticSum(hEventCounts->GetBinError(ibin), relSystErr*hEventCounts->GetBinContent(ibin)));
    for (int ibin = 1; ibin <= hUnfoldedFlux->GetNbinsX(); ++ibin)
      hUnfoldedFlux->SetBinError(ibin, Utilities::QuadraticSum(hUnfoldedFlux->GetBinError(ibin), relSystErr*hUnfoldedFlux->GetBinContent(ibin)));

    std::cout << "Dumping fit input data (speciesSymbol='" << speciesSymbol << "' isTimeAveragedFlux=" << isTimeAveragedFlux << ") ..." << std::endl;

    for (int bin = 0; bin <= hFlux->GetNbinsX() + 1; ++bin)
      std::cout << "bin=" << std::setw(2) << bin << " xLow=" << std::setw(7) << Form("%.2f", hFlux->GetXaxis()->GetBinLowEdge(bin)) << " xUp=" << std::setw(7) << Form("%.2f", hFlux->GetXaxis()->GetBinUpEdge(bin)) << " rawFlux=" << std::setw(20) << Form("%.4f +/- %.4f", hFlux->GetBinContent(bin), hFlux->GetBinError(bin)) << " unfoldedFlux=" << std::setw(20) << Form("%.4f +/- %.4f", hUnfoldedFlux->GetBinContent(bin), hUnfoldedFlux->GetBinError(bin)) << std::endl;

    std::cout << std::endl;
    std::cout << "Start raw flux fit..." << std::endl;

    Modelling::HistogramDataset* dsRaw = new Modelling::HistogramDataset(hFlux);
    dsRaw->SetXRange(Emin, Emax);
    Modelling::BrokenPowerLawModel* modelRaw = RunBrokenPowerLawFit(speciesSymbol.c_str(), dsRaw, true, !isTimeAveragedFlux);

    std::cout << std::endl;
    std::cout << "Start unfolded flux fit..." << std::endl;

    Modelling::HistogramDataset* dsUnfolded = new Modelling::HistogramDataset(hUnfoldedFlux);
    dsUnfolded->SetXRange(Emin, Emax);
    Modelling::BrokenPowerLawModel* modelBayesUnfolded = RunBrokenPowerLawFit(speciesSymbol.c_str(), dsUnfolded, true, !isTimeAveragedFlux);

    std::cout << std::endl;
    std::cout << "Start forward-folding fit..." << std::endl;

    auto* dsForwardFolding = new Modelling::ForwardFoldingDataset(hEventCounts, hAcceptance, hTriggerEfficiency, hMeasuringTime, migrationMatrix);
    dsForwardFolding->SetXRange(Emin, Emax);
    Modelling::BrokenPowerLawModel* modelForwardFolding = RunBrokenPowerLawFit(speciesSymbol.c_str(), dsForwardFolding, true, !isTimeAveragedFlux);

    TGraphAsymmErrors* grRawFluxE3 = new TGraphAsymmErrors;
    TGraphAsymmErrors* grBayesUnfoldedFluxE3 = new TGraphAsymmErrors;
    for (int bin = 1; bin <= hEventCounts->GetNbinsX(); ++bin) {
      double xTrue = Modelling::LaffertyWyatt(reducedBinning.LowEdge(bin), reducedBinning.UpEdge(bin), gAMSLaffertyWyattSpectralIndex);

      grRawFluxE3->SetPoint(bin - 1, xTrue, hFlux->GetBinContent(bin));
      grRawFluxE3->SetPointError(bin - 1, 0.0, 0.0, hFlux->GetBinError(bin), hFlux->GetBinError(bin));

      grBayesUnfoldedFluxE3->SetPoint(bin - 1, xTrue, hUnfoldedFlux->GetBinContent(bin));
      grBayesUnfoldedFluxE3->SetPointError(bin - 1, 0.0, 0.0, hUnfoldedFlux->GetBinError(bin), hUnfoldedFlux->GetBinError(bin));
    }

    Utilities::MultiplyByPower(grRawFluxE3, 3.0);
    Utilities::MultiplyByPower(grBayesUnfoldedFluxE3, 3.0);

    // Plot: Forward-folding / unfolding comparison
    TCanvas* fluxE3Canvas = new TCanvas("fluxE3CanvasFF", Form("%s) E3", speciesSymbol.c_str()));
    fluxE3Canvas->cd();
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.14);

    Color_t bayesColor = kBlue;
    Color_t rawColor = kBlack;

    double ymax = 250.0;
    if (speciesSymbol == "e+")
      ymax = 30.0;
    Utilities::DrawEmptyPlot("", "E (GeV)", "E^{3} #Phi_{e-} (m^{-2} s^{-1} sr^{-1} GeV^{2})", 0.5, 1500., 0., ymax, true);

    TLegend* fluxE3Leg = new TLegend(0.69, 0.68, 0.91, 0.88);

    modelForwardFolding->FluxEpower->Draw("L SAME");
    fluxE3Leg->AddEntry(modelForwardFolding->FluxEpower, "forward-folding fit", "L");

    modelBayesUnfolded->FluxEpower->SetLineColor(bayesColor);
    modelBayesUnfolded->FluxEpower->SetLineStyle(7);
    modelBayesUnfolded->FluxEpower->Draw("L SAME");
    fluxE3Leg->AddEntry(modelBayesUnfolded->FluxEpower, "fit to Bayes-unfolded flux", "L");

    modelRaw->FluxEpower->SetLineColor(rawColor);
    modelRaw->FluxEpower->SetLineStyle(2);
    modelRaw->FluxEpower->SetLineWidth(1);
    modelRaw->FluxEpower->Draw("L SAME");
    fluxE3Leg->AddEntry(modelRaw->FluxEpower, "fit to raw flux", "L");

    grRawFluxE3->SetMarkerColor(rawColor);
    grRawFluxE3->SetLineColor(rawColor);
    grRawFluxE3->SetMarkerStyle(kOpenSquare);
    grRawFluxE3->SetMarkerSize(0.7);
    grRawFluxE3->Draw("P");
    fluxE3Leg->AddEntry(grRawFluxE3, "raw flux", "P");

    grBayesUnfoldedFluxE3->SetMarkerColor(bayesColor);
    grBayesUnfoldedFluxE3->SetLineColor(bayesColor);
    grBayesUnfoldedFluxE3->SetMarkerStyle(kFullCircle);
    grBayesUnfoldedFluxE3->Draw("P");
    fluxE3Leg->AddEntry(grBayesUnfoldedFluxE3, Form("unfolded %s flux (Bayes)", speciesSymbol.c_str()), "P");

    fluxE3Leg->Draw();
    gPad->SetLogx();

    fluxE3Canvas->SaveAs(Form("unfolding_debugging_%s_%s.png", speciesSymbol.c_str(), inputFileSuffix.c_str()));
    fluxE3Canvas->SaveAs(Form("unfolding_debugging_%s_%s.pdf", speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(fluxE3Canvas);

    // Plot: Unfolding / Forward-folding ratio minus one
    TCanvas* dsRatioCanvas = new TCanvas("dsRatioCanvasFF", Form("%s) FF/Unf ratio", speciesSymbol.c_str()));
    gPad->SetLeftMargin(0.13);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.14);
    gPad->SetLogx();

    Utilities::DrawEmptyPlot("", "E (GeV)", "Unfolding / FF - 1", Emin - 0.01, Emax + 0.01, -0.1, 0.1, true);

    TGraph* grRatio_FF_Unfolding = new TGraph;
    Binning::Definition ratioPointsBinning = Binning::Tools::Logarithmic(100, Emin, Emax);
    for (unsigned int i = 0; i < ratioPointsBinning.NumberOfBins(); ++i) {
      double x = ratioPointsBinning.Value(i + 1);
      double fluxFF = modelForwardFolding->Flux->Eval(x);
      double fluxUnf = modelBayesUnfolded->Flux->Eval(x);
      grRatio_FF_Unfolding->SetPoint(i, x, fluxUnf / fluxFF - 1.0);
    }
    grRatio_FF_Unfolding->GetYaxis()->SetRangeUser(-0.1, 0.1);
    grRatio_FF_Unfolding->GetYaxis()->SetNoExponent(kTRUE);
    grRatio_FF_Unfolding->SetMarkerStyle(kFullCircle);
    grRatio_FF_Unfolding->SetMarkerColor(kRed);
    grRatio_FF_Unfolding->SetMarkerSize(0.5);
    grRatio_FF_Unfolding->Draw("AP");

    dsRatioCanvas->SaveAs(Form("unfolding_debugging_ratio_%s_%s.png", speciesSymbol.c_str(), inputFileSuffix.c_str()));
    dsRatioCanvas->SaveAs(Form("unfolding_debugging_ratio_%s_%s.pdf", speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(dsRatioCanvas);

    // Plot: Normalized migration matrix
    TLine* diagonalLine = new TLine(0.40, 0.40, 1300, 1300);
    diagonalLine->SetLineStyle(2);
    diagonalLine->SetLineWidth(4);
    diagonalLine->SetLineColor(kBlack);

    TCanvas* canvasMigrationMatrixNormalized = new TCanvas("canvasMigrationMatrixNormalized", Form("%s) Normalized migration matrix", speciesSymbol.c_str()));
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.18);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.14);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();

    auto* normalizedMigrationMatrix = static_cast<TH2*>(migrationMatrix->Clone());
    Utilities::NormalizeHistogramYSlices(normalizedMigrationMatrix, "", true);

    normalizedMigrationMatrix->GetXaxis()->SetNoExponent(kTRUE);
    normalizedMigrationMatrix->GetXaxis()->SetLabelOffset(normalizedMigrationMatrix->GetXaxis()->GetLabelOffset() + 0.01);

    normalizedMigrationMatrix->GetYaxis()->SetNoExponent(kTRUE);

    normalizedMigrationMatrix->GetZaxis()->SetTitle("Normalized entries");
    normalizedMigrationMatrix->GetZaxis()->SetRangeUser(1e-3, 1);
    normalizedMigrationMatrix->GetZaxis()->SetLabelSize(normalizedMigrationMatrix->GetYaxis()->GetLabelSize());
    normalizedMigrationMatrix->GetZaxis()->SetLabelOffset(normalizedMigrationMatrix->GetYaxis()->GetLabelOffset());
    normalizedMigrationMatrix->GetZaxis()->SetTitleSize(normalizedMigrationMatrix->GetYaxis()->GetTitleSize());
    normalizedMigrationMatrix->GetZaxis()->SetTitleOffset(normalizedMigrationMatrix->GetYaxis()->GetTitleOffset());

    normalizedMigrationMatrix->Draw("col.z");
    diagonalLine->DrawClone();

    canvasMigrationMatrixNormalized->SaveAs(Form("%s_%s_%s.png", canvasMigrationMatrixNormalized->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    canvasMigrationMatrixNormalized->SaveAs(Form("%s_%s_%s.pdf", canvasMigrationMatrixNormalized->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(canvasMigrationMatrixNormalized);

    // Plot: Normalized migration matrix (with cut-off weight)
    TCanvas* canvasMigrationMatrixNormalizedWithCutoffWeight = new TCanvas("canvasMigrationMatrixNormalizedWithCutoffWeight", Form("%s) Normalized migration matrix with cut-off weight", speciesSymbol.c_str()));
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.18);
    gPad->SetTopMargin(0.03);
    gPad->SetBottomMargin(0.14);
    gPad->SetLogx();
    gPad->SetLogy();
    gPad->SetLogz();

    auto* normalizedMigrationMatrixWithCutoffWeight = static_cast<TH2*>(unfolding.MigrationMatrixAfterCutoff()->Clone());
    Utilities::NormalizeHistogramYSlices(normalizedMigrationMatrixWithCutoffWeight, "", true);

    normalizedMigrationMatrixWithCutoffWeight->GetXaxis()->SetNoExponent(kTRUE);
    normalizedMigrationMatrixWithCutoffWeight->GetXaxis()->SetLabelOffset(normalizedMigrationMatrixWithCutoffWeight->GetXaxis()->GetLabelOffset() + 0.01);

    normalizedMigrationMatrixWithCutoffWeight->GetYaxis()->SetNoExponent(kTRUE);

    normalizedMigrationMatrixWithCutoffWeight->GetZaxis()->SetTitle("Normalized entries");
    normalizedMigrationMatrixWithCutoffWeight->GetZaxis()->SetRangeUser(1e-3, 1);
    normalizedMigrationMatrixWithCutoffWeight->GetZaxis()->SetLabelSize(normalizedMigrationMatrixWithCutoffWeight->GetYaxis()->GetLabelSize());
    normalizedMigrationMatrixWithCutoffWeight->GetZaxis()->SetLabelOffset(normalizedMigrationMatrixWithCutoffWeight->GetYaxis()->GetLabelOffset());
    normalizedMigrationMatrixWithCutoffWeight->GetZaxis()->SetTitleSize(normalizedMigrationMatrixWithCutoffWeight->GetYaxis()->GetTitleSize());
    normalizedMigrationMatrixWithCutoffWeight->GetZaxis()->SetTitleOffset(normalizedMigrationMatrixWithCutoffWeight->GetYaxis()->GetTitleOffset());

    normalizedMigrationMatrixWithCutoffWeight->Draw("col.z");
    diagonalLine->DrawClone();

    canvasMigrationMatrixNormalizedWithCutoffWeight->SaveAs(Form("%s_%s_%s.png", canvasMigrationMatrixNormalizedWithCutoffWeight->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    canvasMigrationMatrixNormalizedWithCutoffWeight->SaveAs(Form("%s_%s_%s.pdf", canvasMigrationMatrixNormalizedWithCutoffWeight->GetName(), speciesSymbol.c_str(), inputFileSuffix.c_str()));
    gROOT->GetListOfCanvases()->Remove(canvasMigrationMatrixNormalizedWithCutoffWeight);
  }

  TGraphAsymmErrors* unfoldedEventCounts = static_cast<TGraphAsymmErrors*>(eventCounts->Clone(Form("%sUnfolded", eventCounts->GetName())));

  for (int i = 0; i < unfoldedEventCounts->GetN(); ++i) {
    double x, y;
    unfoldedEventCounts->GetPoint(i, x, y);

    // Keep empty bins.
    if (y == 0) {
      assert(unfoldedEventCounts->GetErrorY(i) == 0);
      continue;
    }

    if (i < int(binningOffset) || i > int(reducedBinning.NumberOfBins() + binningOffset - 1)) {
      unfoldedEventCounts->SetPoint(i, x, 0.0);
      unfoldedEventCounts->SetPointError(i, 0.0, 0.0, 0.0, 0.0);
      continue;
    }

    double xTrue = Modelling::LaffertyWyatt(binning.LowEdge(i + 1),
                                            binning.UpEdge(i + 1),
                                            gAMSLaffertyWyattSpectralIndex);
    assert(std::abs(x - xTrue) < 1e-5);
    assert(std::abs(hUnfoldedEventCounts->GetBinLowEdge(i + 1 - binningOffset) - binning.LowEdge(i + 1)) < 1e-5);

    double ey = hUnfoldedEventCounts->GetBinError(i + 1 - binningOffset);
    unfoldedEventCounts->SetPoint(i, x, hUnfoldedEventCounts->GetBinContent(i + 1 - binningOffset));
    unfoldedEventCounts->SetPointError(i, 0.0, 0.0, ey, ey);
  }

  return unfoldedEventCounts;
}
