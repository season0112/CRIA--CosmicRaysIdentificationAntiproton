#include "TemplateFitterforIntermediateEnergy.hh"

#include "BinningTools.hh"
#include "Utilities.hh"

#include <TCanvas.h>
#include <TClass.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TMarker.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TObject.h>
#include <TStyle.h>

#include <cassert>
#include <cmath>
#include <vector>

#include <QString>
#include <QVector>

#define INFO_OUT_TAG "TemplateFitterforIntermediateEnergy"
#include "debugging.hh"

namespace MYUtilities {

TemplateFitter::TemplateFitter(int printlevel)
  : fHasBeenFit(false)
  , fInputData(0)
  , fInputTemplates()
  , fStartValues()
  , fFixParameter()
  , fFitResult(0)
  , fResultHistos()
  , fMinimizer(0)
  , fPrintLevel(printlevel)
  , fUnboundedFit(false)
  , fFitRangeLow(0)
  , fFitRangeHigh(0)
  , fNDF(-1)
  , fChi2(-1.0) {
}

TemplateFitter::~TemplateFitter() {

  if (fFitResult)
    delete fFitResult;

  for (unsigned int k = 0; k < fResultHistos.size(); ++k)
    delete fResultHistos[k];

  if (fMinimizer)
    delete fMinimizer;
}

bool TemplateFitter::AutoStartValues() {

  if (fInputTemplates.size() != 2) {
    WARN_OUT << "Algorithm only implemented for the case of exactly two templates!" << std::endl;
    return false;
  }
  if (fInputData->GetDimension() != 1) {
    WARN_OUT << "Algorithm only implemented for the case of 1D histograms!" << std::endl;
    return false;
  }

  double cut[2], nEvtOutside[2], nEvtGuess[2], startValue[2];
  for (int i = 0; i < 2; ++i)
    cut[i] = fInputTemplates[i]->GetMean();
  bool isInverted = (cut[1] < cut[0]);
  if (!isInverted) {
    nEvtOutside[0] = fInputData->Integral(1, fInputData->GetXaxis()->FindFixBin(cut[0]));
    nEvtOutside[1] = fInputData->Integral(fInputData->GetXaxis()->FindFixBin(cut[1]), fInputData->GetNbinsX());
  }
  else {
    nEvtOutside[0] = fInputData->Integral(fInputData->GetXaxis()->FindFixBin(cut[0]), fInputData->GetNbinsX());
    nEvtOutside[1] = fInputData->Integral(1, fInputData->GetXaxis()->FindFixBin(cut[1]));
  }
  for (int i = 0; i < 2; ++i)
    nEvtGuess[i] = 2.0 * nEvtOutside[i];

  int nData = fInputData->GetEntries();
  double scale = double(nData) / (nEvtGuess[0] + nEvtGuess[1]);

  for (int i = 0; i < 2; ++i) {
    startValue[i] = nEvtGuess[i] / double(nData) * scale;
    SetStartValue(i, startValue[i]);
  }

  DEBUG_OUT << "cuts: " << cut[0] << " " << cut[1] << ", inv: " << isInverted << ", nEvtOutside: " << nEvtOutside[0] << " " << nEvtOutside[1]
            << ", nData: " << nData << ", scale: " << scale << ", start values: " << startValue[0] << " " << startValue[1]
            << std::endl;
  return true;
}


void TemplateFitter::AddTemplateHistogram(const TH1* hist, double start) {
  assert(hist->GetDimension() == 1);
  assert(hist->IsA()->GetName() == std::string("TH1F") || hist->IsA()->GetName() == std::string("TH1D"));
  fInputTemplates.push_back(hist);
  fStartValues.push_back(start);
  fFixParameter.push_back(false);
}

void TemplateFitter::SetStartValue(unsigned int i, double start) {
  assert(i < fStartValues.size());
  fStartValues.at(i) = start;
}

void TemplateFitter::FixParameter(unsigned int i, double start) {
  assert(i < fStartValues.size());
  assert(i < fFixParameter.size());
  fStartValues.at(i) = start;
  fFixParameter.at(i) = true;
}

void TemplateFitter::FreeParameter(unsigned int i) {
  assert(i < fFixParameter.size());
  fFixParameter.at(i) = false;
}

void TemplateFitter::ResetStartValue(unsigned int i) {
  assert(i < fStartValues.size());
  fStartValues.at(i) = -1.0;
}

int RootColor(int i) {
  static std::vector<int> basicColors;
  if (!basicColors.size()) {
    basicColors.push_back(kAzure);
    basicColors.push_back(kViolet);
    basicColors.push_back(kPink);
    basicColors.push_back(kOrange);
    basicColors.push_back(kGreen);
    basicColors.push_back(kBlue);
    basicColors.push_back(kMagenta);
    basicColors.push_back(kRed);
    basicColors.push_back(kCyan);
  }
  static std::vector<int> subColors;
  if (!subColors.size()) {
    for (int j = 0; j < 14; ++j) {
      if (j < 10) {
        if (j > 7) // skip the light colors
          continue;
        subColors.push_back(-1 * j);
      }
      else
        subColors.push_back(j - 9);
    }
  }

  i *= 2;

  const unsigned int basic = i % basicColors.size();
  unsigned int sub = floor(i / basicColors.size());
  if (sub >= subColors.size())
    sub %= subColors.size();

  assert(basic < basicColors.size());
  assert(sub < subColors.size());

  return basicColors.at(basic) + subColors.at(sub);
}

int TemplateFitter::TemplateColor(int i) const {
  QVector<int> firstColors;
  firstColors << kAzure - 3 << kPink << kGreen << kYellow << kOrange;
  if (i < firstColors.size())
    return firstColors.at(i);
  return RootColor(i);
}

int TemplateFitter::TemplateFillStyle(int i) const {
  QVector<int> fillStyles;
  fillStyles << 3017 << 3018 << 3004 << 3005 << 3006 << 3007 << 3021 << 3022;
  return fillStyles.at(i % fillStyles.size());
}

TCanvas* TemplateFitter::CreateResultDrawing(std::string canvasname, int width, int height, bool defaultStyle) const {

  // FIXME take width and height from current style

  if (!CheckFit())
    return 0;

  TCanvas* graphicsOutput = new TCanvas(canvasname.c_str(), canvasname.c_str(), width, height);
  graphicsOutput->cd();
  //gPad->SetLogy();
  //gPad->SetGrid();

  TLegend* legend = new TLegend(0.20, 0.78, 0.46, 0.88, NULL, "brNDC");
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kBlack);

  TH1* dataCopy = fInputData->DrawCopy("PE");
  dataCopy->SetStats(0); // disable statistic boxes in template canvas
  const int dataColor = kBlack;
  const int dataMarker = 20;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerColor(dataColor);
  dataCopy->SetMarkerStyle(dataMarker);
  dataCopy->SetMarkerSize(1.7);
  //legend->AddEntry(dataCopy, Form("Data (%i)", (int) fInputData->GetEntries()), "L");
  legend->AddEntry(dataCopy, Form("Data (%i)", (int) fInputData->Integral(1, fInputData->GetXaxis()->GetNbins())), "LP"); // To avoid overfow events.

  TH1* fitResultCopy = fFitResult->DrawCopy("hist.same");
  fitResultCopy->SetStats(0);
  const int fitColor = kBlack;
  const int fitStyle = 2;
  fitResultCopy->SetLineWidth(4);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);
  legend->AddEntry(fitResultCopy, Form("Fit Result (#chi^{2}/ndf=%5.2f)", fChi2 / fNDF), "L");

  for (unsigned int i = 0; i < fResultHistos.size(); ++i) {
    TH1* resultCopy = fResultHistos.at(i)->DrawCopy("hist.same");
    resultCopy->SetStats(0);
    if (defaultStyle) {
      resultCopy->SetLineWidth(2);
      const int color = TemplateColor(i);
      const int fillStyle = TemplateFillStyle(i);
      resultCopy->SetFillStyle(fillStyle);
      resultCopy->SetFillColor(color);
      resultCopy->SetLineColor(color);
      resultCopy->SetMarkerColor(color);
    }
    /*
    // split to long titles in two lines
    QString title = resultCopy->GetTitle();
    //QString fittedValue = QString(" #color[%1]{Fitted: %2 #pm %3}")
    QString fittedValue = QString("(%2 #pm %3)")
    //                        .arg(resultCopy->GetLineColor())
                            .arg(GetAbsoluteResult().at(i), 0, 'f', 1, '0')
                            .arg(GetAbsoluteResultError().at(i), 0, 'f', 1, '0');
    title.append(fittedValue);
    if (title.contains(",")) {
      title.replace(",", "}{");
      title.prepend("#splitline{");
      title.append("}");
    }
    legend->AddEntry(resultCopy, qPrintable(title), "F");
    */
  }

  // Redraw  fitResult and data to be on top of the templates
  gStyle->SetErrorX(0);
  fitResultCopy->Draw("hist.same");
  dataCopy->Draw("P E1 same");
  dataCopy->SetLabelFont(62,"xy");
  dataCopy->SetLabelSize(0.05,"xy");
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  dataCopy->GetXaxis()->SetTitle("#Lambda_{TRD}");
  dataCopy->GetYaxis()->SetTitle("N");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.05);
  dataCopy->GetYaxis()->SetTitleSize(0.04);

  TLine line;
  line.SetLineStyle(9);
  line.SetLineWidth(2);
  double lineheight = fInputData->GetMaximum();
  if (fFitRangeLow != 1 && fFitRangeHigh != fInputData->GetNbinsX()) { //fit range set
    line.DrawLine(fInputData->GetBinCenter(fFitRangeLow),  0, fInputData->GetBinCenter(fFitRangeLow), lineheight);
    line.DrawLine(fInputData->GetBinCenter(fFitRangeHigh), 0, fInputData->GetBinCenter(fFitRangeHigh), lineheight);
  }

  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->Draw();

  TText *Text_Pbar = new TText(0.69, 0.40, "Antiprotons");
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(30);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();

  TText *Text_Electron = new TText(0.82, 0.8, "Electrons");
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kRed);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(30);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();

  return graphicsOutput;
}

std::tuple<TH1F*> ReverseXaxis(TH1 *InputHist){
//void ReverseXaxis(TH1 *InputHist){ 
    double binnumber = InputHist->GetNbinsX();
    double low       = -(InputHist->GetBinLowEdge(binnumber)+InputHist->GetBinWidth(binnumber));
    double high      = -(InputHist->GetBinLowEdge(1));
  
    TH1F *OutputHist = new TH1F("", "", binnumber, low, high );
    for (int i=0; i<binnumber; i++){
        OutputHist->SetBinContent(i+1, InputHist->GetBinContent(binnumber-i));
        OutputHist->SetBinError  (i+1, InputHist->GetBinError(binnumber-i));
    }

    OutputHist->SetTitle(InputHist->GetTitle());

    InputHist = OutputHist;
    return {OutputHist};
}

TCanvas* TemplateFitter::CreateResultDrawing_LogY(std::string canvasname, int width, int height, double fitrange, bool defaultStyle) const {

  // FIXME take width and height from current style

  if (!CheckFit())
    return 0;

  TCanvas* graphicsOutput = new TCanvas(canvasname.c_str(), canvasname.c_str(), width, height);
  graphicsOutput->cd();
  gPad->SetLogy();
  //gPad->SetGrid();

  //TLegend* legend = new TLegend(0.20, 0.65, 0.46, 0.88, NULL, "brNDC");
  TLegend* legend = new TLegend(0.60, 0.65, 0.86, 0.88, NULL, "brNDC");
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kBlack);

  TH1* dataCopy = fInputData->DrawCopy("PE");
  std::tie(dataCopy) = ReverseXaxis(dataCopy);
  dataCopy->Draw("PE");
  dataCopy->SetTitle("");
  dataCopy->SetStats(0); // disable statistic boxes in template canvas
  const int dataColor = kBlack;
  const int dataMarker = 15;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerColor(dataColor);
  dataCopy->SetMarkerStyle(dataMarker);
  dataCopy->SetMarkerSize(1.7);
  //legend->AddEntry(dataCopy, Form("Data (%i)", (int) fInputData->GetEntries()), "L");
  //legend->AddEntry(dataCopy, Form("Data (%i)", (int) fInputData->Integral(1, fInputData->GetXaxis()->GetNbins())), "LP"); // To avoid overfow events.
  legend->AddEntry(dataCopy, Form("Data (%i)", (int) fInputData->Integral(1, fInputData->GetXaxis()->FindFixBin(fitrange))), "LP"); 

  TH1* fitResultCopy = fFitResult->DrawCopy("hist.same");
  std::tie(fitResultCopy) = ReverseXaxis(fitResultCopy);
  fitResultCopy->Draw("hist.same");
  fitResultCopy->SetStats(0);
  const int fitColor = kBlack;
  const int fitStyle = 2;
  fitResultCopy->SetLineWidth(4);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);

  for (unsigned int i = 0; i < fResultHistos.size(); ++i) {
    TH1* resultCopy = fResultHistos.at(i)->DrawCopy("hist.same");
    std::tie(resultCopy) = ReverseXaxis(resultCopy);
    resultCopy->Draw("hist.same");
    resultCopy->SetStats(0);
    if (defaultStyle) {
      resultCopy->SetLineWidth(2);
      const int color = TemplateColor(i);
      const int fillStyle = TemplateFillStyle(i);
      resultCopy->SetFillStyle(fillStyle);
      resultCopy->SetFillColor(color);
      resultCopy->SetLineColor(color);
      resultCopy->SetMarkerColor(color);
    }
    QString title = resultCopy->GetTitle();
    QString fittedValue = QString("(%2 #pm %3)")
                            .arg(GetAbsoluteResult().at(i), 0, 'f', 1, '0')
                            .arg(GetAbsoluteResultError().at(i), 0, 'f', 1, '0');
    title.append(fittedValue);
    if (title.contains(",")) {
      title.replace(",", "}{");
      title.prepend("#splitline{");
      title.append("}");
    }
    legend->AddEntry(resultCopy, qPrintable(title), "F");
  }

  // Redraw  fitResult and data to be on top of the templates
  gStyle->SetErrorX(0);
  fitResultCopy->Draw("hist.same");
  dataCopy     ->Draw("P E1 same");
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  //dataCopy->GetXaxis()->SetTitle("Sign(R) #bullet #Lambda_{TRD}");
  dataCopy->GetXaxis()->SetTitle("#Lambda_{TRD}");
  dataCopy->GetYaxis()->SetTitle("N");
  dataCopy->SetLabelFont(62,"xy");
  dataCopy->SetLabelSize(0.045,"xy");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.045);
  dataCopy->GetYaxis()->SetTitleSize(0.045);

  // Fit Range Line (Not working? FIXME)
  TLine line;
  line.SetLineStyle(9);
  line.SetLineWidth(2);
  double lineheight = fInputData->GetMaximum();
  if (fFitRangeLow != 1 && fFitRangeHigh != fInputData->GetNbinsX()) { //fit range set
    line.DrawLine(fInputData->GetBinCenter(fFitRangeLow), 0, fInputData->GetBinCenter(fFitRangeLow), lineheight);
    line.DrawLine(fInputData->GetBinCenter(fFitRangeHigh), 0, fInputData->GetBinCenter(fFitRangeHigh), lineheight);
  }
  TLine line2;
  double lineheight2 = gPad->GetUymax();
  line2.SetLineStyle(9);
  line2.SetLineWidth(2);
  line2.DrawLine(-fitrange, 0, -fitrange, lineheight);

  //Option1: Calculate myself.
  double chi2 = 0;
  double tem  = 0;
  for (unsigned int i = 0; i < dataCopy->GetNbinsX(); i++) {
      tem  = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / fitResultCopy->GetBinContent(i+1);
      chi2 = chi2 + tem;
  }
  // Option2:Taken from ROOT.
  double chi2_v2 = dataCopy->Chi2Test(fitResultCopy, "UU CHI2");
  std::cout<< "chi2:" << chi2 <<std::endl;
  std::cout<< "chi2_v2:" << chi2_v2 <<std::endl;
  std::cout<< "Chi2():" << Chi2() <<std::endl;
  legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Chi2(), NDF(), Chi2() / NDF()));
    

  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->SetBorderSize(0);
  legend->Draw();

  /*
  //TText *Text_Pbar = new TText(0.48, 0.50, "Antiprotons"); // Time averaged setting
  TText *Text_Pbar = new TText(0.50, 0.45, "Antiprotons"); // Time dependent setting
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(30);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();

  TText *Text_Electron = new TText(0.75, 0.6, "Electrons"); // Time averaged and dependent setting
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kRed);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(30);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();
  */

  return graphicsOutput;
}



void TemplateFitter::Clear() {

  fHasBeenFit = false;
  if (fInputData)
    fInputData = 0;
  fInputTemplates.clear();

  fStartValues.clear();
  fFixParameter.clear();

  return;
}

TH1F* TemplateFitter::CloneCombinedResultHistogram(const char* hname = "") const {

  if (!CheckFit())
    return 0;

  return (TH1F*) fFitResult->Clone(hname);
}

TCanvas* TemplateFitter::CreateContourPlot(int par1, int par2, std::string canvasname, int width, int height) const {

  if (!CheckFit())
    return 0;

  gStyle->SetOptStat(0);
  TCanvas* graphicsOutput = new TCanvas(canvasname.c_str(), canvasname.c_str(), width, height);
  graphicsOutput->cd();

  const int Nbins = 800;
  char htitle[256];
  sprintf(htitle, ";Parameter %d;Parameter %d", par1, par2);
  TH2D* hist = new TH2D("contourplot", htitle, Nbins, fRelative.at(par1) - 4 * fRelativeErrorLower.at(par1), fRelative.at(par1) + 4 * fRelativeErrorUpper.at(par1), Nbins, fRelative.at(par2) - 4 * fRelativeErrorLower.at(par2), fRelative.at(par2) + 4 * fRelativeErrorUpper.at(par2));

  if (fInputTemplates.size() < 2) {
    WARN_OUT << "Less than two parameters!" << std::endl;
    return 0;
  }
  if (par1 == par2) {
    WARN_OUT << "Cannot correlate parameter with itself!" << std::endl;
    return 0;
  }
  if (par1 < 0 || par2 < 0) {
    WARN_OUT << "Some parameter index is smaller than zero!" << std::endl;
    return 0;
  }
  if (par1 >= (int) fInputTemplates.size() || par2 >= (int) fInputTemplates.size()) {
    WARN_OUT << "Some parameter index too large!" << std::endl;
    return 0;
  }

  std::vector<double> pars = fRelative;
  for (int x = 1; x <= Nbins; x++)
    for (int y = 1; y <= Nbins; y++) {

      pars.at(par1) = hist->GetXaxis()->GetBinCenter(x);
      pars.at(par2) = hist->GetYaxis()->GetBinCenter(y);
      int NDF;
      hist->SetBinContent(x, y, ManualChiSqr(NDF, pars));
    }

  hist->DrawCopy("colz");

  double levels[3] = { fChi2 + 1, fChi2 + 4, fChi2 + 9 };
  hist->SetContour(3, levels);
  hist->SetLineColor(1);
  hist->SetLineWidth(3);
  hist->DrawCopy("cont2 same");

  TMarker mark;
  mark.SetMarkerStyle(20);
  mark.SetMarkerColor(1);
  mark.SetMarkerSize(2.0);
  mark.DrawMarker(fRelative.at(par1), fRelative.at(par2));

  char textline[256];
  sprintf(textline, "#chi^{2}/NDF = %.2f/%d ", fChi2, fNDF);
  TLatex text;
  text.SetNDC();
  text.SetTextSize(0.03);
  text.DrawLatex(0.15, 0.85, textline);

  return graphicsOutput;
}


double TemplateFitter::Fit(int i, double lowBound, double highBound) {

  //reset things in case there was a previous fit with the same templates.
  fHasBeenFit = false;

  if (fInputTemplates.size() == 0) {
    WARN_OUT << "No templates added! Abort Fit!" << std::endl;
    return -4.0;
  }

  if (!fInputData) {
    WARN_OUT << "No data added! Abort Fit!" << std::endl;
    return -4.0;
  }

  if (lowBound != 0 || highBound != 0) { //fit range specified
    fFitRangeLow = std::max(fInputData->FindFixBin(lowBound), 1);
    fFitRangeHigh = std::min(fInputData->FindFixBin(highBound), fInputData->GetNbinsX());
  }
  else {
    fFitRangeLow = 1;
    fFitRangeHigh = fInputData->GetNbinsX();
  }

  fRelative.clear();
  fAbsolute.clear();
  fAbsoluteAllRange.clear();
  fRelativeError.clear();
  fAbsoluteError.clear();
  fAbsoluteAllRangeError.clear();
  fRelativeErrorUpper.clear();
  fAbsoluteErrorUpper.clear();
  fAbsoluteAllRangeErrorUpper.clear();
  fRelativeErrorLower.clear();
  fAbsoluteErrorLower.clear();
  fAbsoluteAllRangeErrorLower.clear();
  fGlobalCorrelationFactor.clear();

  fNentries.clear();
  fBinWidth.clear();
  fDataSet.clear();
  fRelError.clear();

  if (fFitResult)
    delete fFitResult;
  fFitResult = 0;
  for (unsigned int k = 0; k < fResultHistos.size(); ++k)
    delete fResultHistos[k];
  fResultHistos.clear();

  if (fMinimizer) {
    delete fMinimizer;
    fMinimizer = 0;
  }

  if (!TestHistogramCompatibility()) {
    WARN_OUT << "Histogram compatibility failed... Fit aborted!" << std::endl;
    return -2.00;
  }

  //Fill fDataSet and other vectors needed for the fit.
  //Also the fit range is taken care of in here
  PrepareFitData();

  const double default_start_value = 1.0 / (double) fInputTemplates.size();
  for (unsigned int k = 0; k < fStartValues.size(); k++) {
    if (fStartValues.at(k) < 0.0 || fStartValues.at(k) > 1.0) {
      fStartValues.at(k) = default_start_value;
    }
  }

  double ReducedChi;
  if (i == 0)
    ReducedChi = ExecuteFit(chisqr_fit_func);
  else if (i == 1) {
    ReducedChi = ExecuteFit(llh_fit_func);
  }
  else {
    WARN_OUT << "No valid fit strategy choosen!" << std::endl;
    ReducedChi = -5.0;
  }

  //assemble fFitResult
  for (unsigned int k = 0; k < fInputTemplates.size(); k++) {

    fAbsolute.push_back(fRelative.at(k) * fNentries.at(0));
    fAbsoluteError.push_back(fRelativeError.at(k) * fNentries.at(0));
    fAbsoluteErrorUpper.push_back(fRelativeErrorUpper.at(k) * fNentries.at(0));
    fAbsoluteErrorLower.push_back(fRelativeErrorLower.at(k) * fNentries.at(0));

    DEBUG_OUT << "k=" << k << ", rel: " << fRelative[k] << " entries " << fNentries[0] << std::endl;

    fAbsoluteAllRange.push_back(fAbsolute.back() * fInputTemplates.at(k)->Integral() / fNentries.at(k + 1));
    fAbsoluteAllRangeError.push_back(fAbsoluteError.back() * fInputTemplates.at(k)->Integral() / fNentries.at(k + 1));
    fAbsoluteAllRangeErrorUpper.push_back(fAbsoluteErrorUpper.back() * fInputTemplates.at(k)->Integral() / fNentries.at(k + 1));
    fAbsoluteAllRangeErrorLower.push_back(fAbsoluteErrorLower.back() * fInputTemplates.at(k)->Integral() / fNentries.at(k + 1));

    DEBUG_OUT << "absolute " << fAbsolute.back() << " * " << fInputTemplates.at(k)->Integral() << " / " << fNentries.at(k + 1)
              << " = " << fAbsolute.back() * fInputTemplates.at(k)->Integral() / fNentries.at(k + 1) << std::endl;

    fResultHistos.push_back((TH1*) fInputTemplates.at(k)->Clone());

    double scaleFactor = fAbsoluteAllRange.back() / fResultHistos.back()->Integral();
    DEBUG_OUT << "scale " << fAbsoluteAllRange.back() << " / " << fResultHistos.back()->Integral()
              << " = " << scaleFactor << std::endl;

    DEBUG_OUT << "integral 1: " << fResultHistos.back()->Integral() << std::endl;

    fResultHistos.back()->Scale(scaleFactor);

    DEBUG_OUT << "integral 2: " << fResultHistos.back()->Integral() << std::endl;

    if (k == 0) {
      fFitResult = (TH1*) fResultHistos.back()->Clone("hFitResult");
      fFitResult->SetTitle("Fit result");
    }
    else
      fFitResult->Add(fResultHistos.back());
  }

  return ReducedChi;
}

void TemplateFitter::PrepareFitData() {
  //fill data to fit and normalize
  //respect fit range (only done for 1D fit)

  const int mNbinsX = fInputData->GetNbinsX();
  const int mNbinsY = fInputData->GetNbinsY();

  std::vector<double> BinWidthX;
  for (int i = 1; i <= mNbinsX; i++)
    BinWidthX.push_back(fInputData->GetXaxis()->GetBinWidth(i));

  std::vector<double> BinWidthY;
  for (int i = 1; i <= mNbinsY; i++)
    BinWidthY.push_back(fInputData->GetYaxis()->GetBinWidth(i));

  const int dimension = fInputData->GetDimension();

  std::vector<double> kDataSet;
  std::vector<double> kRelError;
  if (dimension == 1) {

    fNentries.push_back(fInputData->Integral(fFitRangeLow, fFitRangeHigh));
    for (unsigned int i = 0; i < fInputTemplates.size(); i++)
      fNentries.push_back(fInputTemplates.at(i)->Integral(fFitRangeLow, fFitRangeHigh));

    //fill from data
    for (int i = fFitRangeLow; i <= fFitRangeHigh; ++i) {
      kDataSet.push_back((double) (fInputData->GetBinContent(i) / BinWidthX.at(i - 1) / fNentries.at(0)));
      if ((double) fInputData->GetBinContent(i) != 0)
        kRelError.push_back((double) fInputData->GetBinError(i) / (double) fInputData->GetBinContent(i));
      else
        kRelError.push_back(0);
      fBinWidth.push_back(BinWidthX.at(i - 1));
    }
    fDataSet.push_back(kDataSet);
    fRelError.push_back(kRelError);
    kDataSet.clear();
    kRelError.clear();

    //Fill from templates
    for (unsigned int k = 0; k < fInputTemplates.size(); ++k) {
      for (int i = fFitRangeLow; i <= fFitRangeHigh; ++i) {
        kDataSet.push_back((double) (fInputTemplates.at(k)->GetBinContent(i) / BinWidthX.at(i - 1) / fNentries.at(k + 1)));
        if ((double) fInputTemplates.at(k)->GetBinContent(i) != 0)
          kRelError.push_back((double) fInputTemplates.at(k)->GetBinError(i) / (double) fInputTemplates.at(k)->GetBinContent(i));
        else
          kRelError.push_back(0);
      }
      fDataSet.push_back(kDataSet);
      fRelError.push_back(kRelError);
      kDataSet.clear();
      kRelError.clear();
    }
  }
  else if (dimension == 2) {

    fNentries.push_back(fInputData->Integral());
    for (unsigned int i = 0; i < fInputTemplates.size(); i++)
      fNentries.push_back(fInputTemplates.at(i)->Integral());

    //fill from data
    for (int i = 1; i <= mNbinsX; ++i)
      for (int j = 1; j <= mNbinsY; ++j) {
        kDataSet.push_back((double) (fInputData->GetBinContent(i, j) / BinWidthX.at(i - 1) / BinWidthY.at(j - 1) / fNentries.at(0)));

        if ((double) fInputData->GetBinContent(i, j) != 0)
          kRelError.push_back((double) fInputData->GetBinError(i, j) / (double) fInputData->GetBinContent(i, j));
        else
          kRelError.push_back(0);
        fBinWidth.push_back(BinWidthX.at(i - 1) * BinWidthY.at(j - 1));
      }
    fDataSet.push_back(kDataSet);
    fRelError.push_back(kRelError);
    kDataSet.clear();
    kRelError.clear();

    //Fill from templates
    for (unsigned int k = 0; k < fInputTemplates.size(); ++k) {
      for (int i = 1; i <= mNbinsX; ++i)
        for (int j = 1; j <= mNbinsY; ++j) {
          kDataSet.push_back((double) (fInputTemplates.at(k)->GetBinContent(i, j) / BinWidthX.at(i - 1) / BinWidthY.at(j - 1) / fNentries.at(k + 1)));
          if ((double) fInputTemplates.at(k)->GetBinContent(i, j) != 0)
            kRelError.push_back((double) fInputTemplates.at(k)->GetBinError(i, j) / (double) fInputTemplates.at(k)->GetBinContent(i, j));
          else
            kRelError.push_back(0);
        }
      fDataSet.push_back(kDataSet);
      fRelError.push_back(kRelError);
      kDataSet.clear();
      kRelError.clear();
    }
  }
  else {
    assert(0);
  }

  //for the FIT now the vectors fDataSet, fRelError, fBinWidth, and fNentries can be used
  return;
}

double TemplateFitter::ExecuteFit(void (*function)(int&, double*, double&, double*, int)) {

  /* In List: Entry:    Designation
   *          0:        fNentries
   *          1:        fBinWidth
   *          2 until including 2 + NumberOfTemplates: contents of fDataSet
   *          3 until including 3 + 2 * NumberOfTemplates: contents fRelError
   *  DataSet at entry k corresponds to RelError at entry k+NumberOfTemplates+1
   */
  TemplateFitObject fitData;
  fitData.vec.push_back(fNentries);
  fitData.vec.push_back(fBinWidth);
  for (unsigned int i = 0; i <= fInputTemplates.size(); ++i)
    fitData.vec.push_back(fDataSet[i]);

  for (unsigned int i = 0; i <= fInputTemplates.size(); ++i)
    fitData.vec.push_back(fRelError[i]);

  if (fMinimizer)
    delete fMinimizer;

  fMinimizer = new TMinuit(static_cast<int>(fInputTemplates.size()));
  fMinimizer->SetPrintLevel(fPrintLevel);
  fMinimizer->SetFCN(function);
  fMinimizer->SetObjectFit(&fitData);

  int fitStatus = 0;
  double arglist[2] = { 0 };

  arglist[0] = 1.0; // Likelihood-fit.
  gMinuit->mnexcm("SET ERR", arglist, 1, fitStatus);
  if (fitStatus != 0) {
    WARN_OUT << "Template fit initialization failed!" << std::endl;
    return -1;
  }

  // Define parameters.
  const double parameterRelativeStepSize = 0.01;
  const double parameterLowerLimit = 0.0;
  const double parameterUpperLimit = fUnboundedFit ? 0.0 : 1.2;
  for (unsigned int i = 0; i < fInputTemplates.size(); i++) {
    char parameterName[32];
    sprintf(parameterName, "Amplitude_%d", i);
    fMinimizer->mnparm(i, parameterName, fStartValues.at(i), parameterRelativeStepSize * fStartValues.at(i), parameterLowerLimit, parameterUpperLimit, fitStatus);
    if (fitStatus != 0) {
      WARN_OUT << "Defining parameters failed." << std::endl;
      return -1;
    }
  }

  // Fix individual parameters.
  for (unsigned int i = 0; i < fInputTemplates.size(); ++i) {
    if (!fFixParameter.at(i))
      continue;
    fMinimizer->mnfixp(i, fitStatus);
    if (fitStatus != 0) {
      WARN_OUT << "Fixing parameters failed." << std::endl;
      return -1;
    }
  }

  const unsigned int minimizationCalls = 10000;    // Max. number of minimization calls
  const double parameterMaxTolerance = 0.01;       // Stop after reaching estimated vertical distance to minimum (EDM).
  const double parameterMaxToleranceCoarse = 0.02; // Stop after reaching estimated vertical distance to minimum (EDM).
  arglist[0] = minimizationCalls;
  arglist[1] = parameterMaxTolerance;

  // Refine initial start values using SIMPLEX method.
  fMinimizer->mnexcm("SIMPLEX", arglist, 2, fitStatus); // Status is ignored for the initial start value refining.

  // Execute MIGRAD minimization.
  fMinimizer->mnexcm("MIGRAD", arglist, 2, fitStatus);
  if (fitStatus != 0)
    fMinimizer->mnexcm("MIGRAD", arglist, 2, fitStatus); // Retry fit.

  if (fitStatus != 0) {
    arglist[1] = parameterMaxToleranceCoarse;
    fMinimizer->mnexcm("MIGRAD", arglist, 2, fitStatus); // Retry fit.
    arglist[1] = parameterMaxTolerance;
  }

  // Refine error estimation using MINOS technique, but only if MIGRAD method worked flawlessly.
  if (fitStatus == 0)
    fMinimizer->mnexcm("MINOS", arglist, 2, fitStatus);

  // Retrieve MIGRAD error matrix.
  unsigned int nDim = fInputTemplates.size();
  fErrorMatrix.clear();
  fErrorMatrix.resize(nDim * nDim);
  fMinimizer->mnemat(fErrorMatrix.data(), nDim);

  // Retrieve fit results.
  double value, error, error_parabola, errorUp, errorLow, correlation;
  for (unsigned int k = 0; k < fInputTemplates.size(); ++k) {
    fMinimizer->GetParameter(k, value, error);
    fMinimizer->mnerrs(k, errorUp, errorLow, error_parabola, correlation);

    fRelative.push_back(value);
    fRelativeError.push_back(error);
    fRelativeErrorUpper.push_back(errorUp);
    fRelativeErrorLower.push_back(std::abs(errorLow));

    fGlobalCorrelationFactor.push_back(correlation);
  }

  // Manually calculate goodness-of-fit parameter.
  fChi2 = ManualChiSqr(fNDF, fRelative);

  if (fitStatus == 0) {
    fHasBeenFit = true;
    return fNDF <= 0 ? -3.0 : fChi2 / fNDF;
  }

  return -1.0;
}


void TemplateFitter::chisqr_fit_func(int&, double*, double& f, double* par, int) {

  /* In List: Entry:    Designation
   *          0:        fNentries
   *          1:        fBinWidth
   *          2 bis einschliesslich 2 + Anzahl Templates: fDataSet Inhalte
   *          3 + Anzahl Templates bis einschliesslich 3+2*AnzahlTemplates: fRelError Inhalte
   *  Zu DataSet[i] bei index k gehoert RelError[i] bei k+Number of Templates+1
   */
  TemplateFitObject* list = (TemplateFitObject*) gMinuit->GetObjectFit();
  int NumberOfTemplates = (list->vec.size() - 4) / 2;
  const int index_data_data = 2;
  const int index_data_error = 3 + NumberOfTemplates;

  double Chisq = 0;

  //Calculate Chisquare!
  for (unsigned int i = 0; i < list->vec[index_data_data].size(); i++) {

    double TempSum = 0;
    for (int k = 0; k < NumberOfTemplates; k++)
      TempSum += par[k] * list->vec[index_data_data + 1 + k].at(i);

    if (list->vec[index_data_data].at(i) == 0 || TempSum == 0)
      continue;

    double nominator = std::pow(TempSum - list->vec[index_data_data].at(i), 2);

    double sigsqr = std::pow(list->vec[index_data_error].at(i) * list->vec[index_data_data].at(i), 2);
    for (int k = 0; k < NumberOfTemplates; k++)
      sigsqr += std::pow(list->vec[index_data_error + 1 + k].at(i) * list->vec[index_data_data + 1 + k].at(i) * std::pow(par[k], 2), 2);

    Chisq += nominator / sigsqr;
  }

  f = Chisq;
}

void TemplateFitter::llh_fit_func(int&, double*, double& f, double* par, int) {

  /* In List: Entry:    Designation
   *          0:        fNentries
   *          1:        fBinWidth
   *          2 bis einschliesslich 2 + Anzahl Templates: fDataSet Inhalte
   *          3 + Anzahl Templates bis einschliesslich 3+2*AnzahlTemplates: fRelError Inhalte
   *  Zu DataSet[i] bei index k gehoert RelError[i] bei k+Number of Templates+1
   */
  TemplateFitObject* list = (TemplateFitObject*) gMinuit->GetObjectFit();
  int NumberOfTemplates = (list->vec.size() - 4) / 2;
  const int index_Nentries = 0;
  const int index_BinWidth = 1;
  const int index_data_data = 2;

  //Calculate Log Likelihood! (Binned extended version: Cowan p.89 Eq.: 6.46)
  // -log(Llh)= -(Sum(n_i*log(m_i))-m_tot), m=model, n=data
  //Relative parameters need to be converted to absolute event numbers!

  double m_tot = 0;
  for (int i = 0; i < NumberOfTemplates; i++)
    for (unsigned int j = 0; j < list->vec[index_data_data].size(); j++)
      m_tot += par[i] * list->vec[index_data_data + 1 + i].at(j) * list->vec[index_BinWidth].at(j) * list->vec[index_Nentries].at(0);

  double LogL = 0;
  for (unsigned int i = 0; i < list->vec[index_data_data].size(); i++) {

    double logi = 0;
    double TempSum = 0;
    for (int k = 0; k < NumberOfTemplates; k++)
      TempSum += list->vec[index_data_data + 1 + k].at(i);

    if (TempSum != 0) { //only use points for which template information exists
      for (int k = 0; k < NumberOfTemplates; k++)
        logi += par[k] * list->vec[index_data_data + 1 + k].at(i) * list->vec[index_Nentries].at(0) * list->vec[index_BinWidth].at(i);


      if (logi <= 0.0)
        logi = 1.e-10;

      logi = std::log(logi);
      LogL += list->vec[index_data_data].at(i) * list->vec[index_Nentries].at(0) * list->vec[index_BinWidth].at(i) * logi;
    }
  }
  f = -2 * (LogL - m_tot);
}

bool TemplateFitter::TestHistogramCompatibility() const {

  for (unsigned int i = 0; i < fInputTemplates.size(); ++i) {
    //check minimum, maximum, bin number and bin width
    if (fInputTemplates.at(i)->GetNbinsX() != fInputData->GetNbinsX())
      return false;

    if (Binning::FromTAxis(fInputTemplates.at(i)->GetXaxis()) != Binning::FromTAxis(fInputData->GetXaxis()))
      return false;
  }

  return true;
}

const std::vector<double>& TemplateFitter::GetRelativeResult() const {

  CheckFit();
  return fRelative;
}

const std::vector<double>& TemplateFitter::GetAbsoluteResult() const {

  CheckFit();
  return fAbsolute;
}

const std::vector<double>& TemplateFitter::GetAbsoluteFullRangeResult() const {

  CheckFit();
  return fAbsoluteAllRange;
}

const std::vector<double>& TemplateFitter::GetRelativeResultError() const {

  CheckFit();
  return fRelativeError;
}

const std::vector<double>& TemplateFitter::GetAbsoluteResultError() const {

  CheckFit();
  return fAbsoluteError;
}

const std::vector<double>& TemplateFitter::GetAbsoluteFullRangeResultError() const {

  CheckFit();
  return fAbsoluteAllRangeError;
}

const std::vector<double>& TemplateFitter::GetRelativeResultErrorUpper() const {

  CheckFit();
  return fRelativeErrorUpper;
}

const std::vector<double>& TemplateFitter::GetAbsoluteResultErrorUpper() const {

  CheckFit();
  return fAbsoluteErrorUpper;
}

const std::vector<double>& TemplateFitter::GetAbsoluteFullRangeResultErrorUpper() const {

  CheckFit();
  return fAbsoluteAllRangeErrorUpper;
}

const std::vector<double>& TemplateFitter::GetRelativeResultErrorLower() const {

  CheckFit();
  return fRelativeErrorLower;
}

const std::vector<double>& TemplateFitter::GetAbsoluteResultErrorLower() const {

  CheckFit();
  return fAbsoluteErrorLower;
}

const std::vector<double>& TemplateFitter::GetAbsoluteFullRangeResultErrorLower() const {

  CheckFit();
  return fAbsoluteAllRangeErrorLower;
}

const std::vector<double>& TemplateFitter::GetGlobalCorrelationFactors() const {

  CheckFit();
  return fGlobalCorrelationFactor;
}

const std::vector<double>& TemplateFitter::GetErrorMatrix() const {

  CheckFit();
  return fErrorMatrix;
}

bool TemplateFitter::CheckFit() const {
  if (!fHasBeenFit)
    DEBUG_OUT << "No fit was done or fit not successful!" << std::endl;
  return fHasBeenFit;
}

double TemplateFitter::ManualChiSqr(int& NDF, std::vector<double> par) const {

  double Chisq = 0;
  NDF = 0;

  //Calculate Chisquare!
  for (unsigned int i = 0; i < fDataSet[0].size(); i++) {

    double TempSum = 0;
    for (unsigned int k = 0; k < fInputTemplates.size(); k++)
      TempSum += par.at(k) * fDataSet[k + 1].at(i);

    if (TempSum == 0 || fDataSet[0].at(i) == 0)
      continue;

    double nominator = std::pow(TempSum - fDataSet[0].at(i), 2);
    double sigsqr = std::pow(fDataSet[0].at(i) * fRelError[0].at(i), 2);
    for (unsigned int k = 0; k < fInputTemplates.size(); k++)
      sigsqr += std::pow(fDataSet[k + 1].at(i) * std::pow(par.at(k), 2) * fRelError[k + 1].at(i), 2); //poisson error

    Chisq += nominator / sigsqr;
    NDF += 1;
  }

  NDF -= fInputTemplates.size();

  return Chisq;
}

} // namespace MYUtilities
