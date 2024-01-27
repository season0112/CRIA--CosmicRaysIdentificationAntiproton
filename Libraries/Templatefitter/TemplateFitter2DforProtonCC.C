
#include "Statistics.hh"
#include "TemplateFitter2DforProtonCC.hh"

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <TLatex.h>
#include <TLine.h>

#include <cassert>
#include <cmath>

#include <QString>

#define INFO_OUT_TAG "TemplateFitter2DforProtonCC"
#include "debugging.hh"

namespace MYUtilities {

TemplateFitter2D::TemplateFitter2D(int printlevel)
  : TemplateFitter(printlevel) {
}

TemplateFitter2D::~TemplateFitter2D() {
}

void TemplateFitter2D::AddTemplateHistogram(const TH1* hist, double start) {
  assert(hist->GetDimension() == 2);
  fInputTemplates.push_back(hist);
  fStartValues.push_back(start);
  fFixParameter.push_back(false);
}

TCanvas* TemplateFitter2D::CreateResultDrawing(std::string canvasname, int width, int height) const {

  if (!fHasBeenFit) {
    WARN_OUT << "No fit was done!" << std::endl;
    return 0;
  }

  gStyle->SetOptStat(0);
  TCanvas* graphicsOutput = new TCanvas(canvasname.c_str(), canvasname.c_str(), width, height);
  graphicsOutput->cd();
  graphicsOutput->Divide(2, 2);


  graphicsOutput->cd(1); //data
  gPad->SetLogz();
  //gPad->SetGrid();
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  TH1* dataCopy = fInputData->DrawCopy("COLZ");
  dataCopy->SetStats(0); // disable statistic boxes in template canvas
  dataCopy->GetXaxis()->SetTitle("#Lambda_{CC}");
  dataCopy->GetYaxis()->SetTitle("#Lambda_{TRD}");
  dataCopy->GetZaxis()->SetTitle("N");
  dataCopy->SetLabelFont(62,"xyz");
  dataCopy->SetLabelSize(0.05,"xyz");
  dataCopy->SetTitle("");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetZaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.055);
  dataCopy->GetYaxis()->SetTitleSize(0.055);
  dataCopy->GetZaxis()->SetTitleSize(0.055);
  dataCopy->GetZaxis()->SetTitleOffset(0.4);

  TText *Text_Pbar = new TText(0.78, 0.80, "Antiprotons"); 
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(17);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();
  TText *Text_Electron = new TText(0.78, 0.23, "Electrons"); 
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kGreen);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(17);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();
  TText *Text_CCP = new TText(0.40, 0.75, "Charge confused protons"); 
  Text_CCP->SetTextAlign(22);
  Text_CCP->SetTextColor(kRed);
  Text_CCP->SetTextFont(43);
  Text_CCP->SetTextSize(17);
  Text_CCP->SetNDC(1);
  Text_CCP->Draw();
  TLine *line_Pbar = new TLine(0.88, 0.58, 0.80, 0.750);
  line_Pbar->SetLineColor(kBlack);
  line_Pbar->SetLineWidth(2);
  line_Pbar->SetNDC(1);
  line_Pbar->Draw();
  TLine *line_Electron = new TLine(0.88, 0.37, 0.80, 0.25);
  line_Electron->SetLineColor(kBlack);
  line_Electron->SetLineWidth(2);
  line_Electron->SetNDC(1);
  line_Electron->Draw();
  TLine *line_CCP = new TLine(0.20, 0.56, 0.35, 0.70);
  line_CCP->SetLineColor(kBlack);
  line_CCP->SetLineWidth(2);
  line_CCP->SetNDC(1);
  line_CCP->Draw();


  graphicsOutput->cd(2); //fit result
  gPad->SetLogz();
  //gPad->SetGrid();
  fFitResult->SetStats(0);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  TH1* fitResultCopy = fFitResult->DrawCopy("COLZ");
  fitResultCopy->GetXaxis()->SetTitle("#Lambda_{CC}");
  fitResultCopy->GetYaxis()->SetTitle("#Lambda_{TRD}");
  fitResultCopy->GetZaxis()->SetTitle("N");
  fitResultCopy->SetLabelFont(62,"xyz");
  fitResultCopy->SetLabelSize(0.05,"xyz");
  fitResultCopy->SetTitle("");
  fitResultCopy->GetXaxis()->SetTitleFont(62);
  fitResultCopy->GetYaxis()->SetTitleFont(62);
  fitResultCopy->GetZaxis()->SetTitleFont(62);
  fitResultCopy->GetXaxis()->SetTitleSize(0.055);
  fitResultCopy->GetYaxis()->SetTitleSize(0.055);
  fitResultCopy->GetZaxis()->SetTitleSize(0.055);
  fitResultCopy->GetZaxis()->SetTitleOffset(0.4);

  TText *Text_Pbar_2 = new TText(0.78, 0.80, "Antiprotons");
  Text_Pbar_2->SetTextAlign(22);
  Text_Pbar_2->SetTextColor(kBlue);
  Text_Pbar_2->SetTextFont(43);
  Text_Pbar_2->SetTextSize(17);
  Text_Pbar_2->SetNDC(1);
  Text_Pbar_2->Draw();
  TText *Text_Electron_2 = new TText(0.78, 0.23, "Electrons");
  Text_Electron_2->SetTextAlign(22);
  Text_Electron_2->SetTextColor(kGreen);
  Text_Electron_2->SetTextFont(43);
  Text_Electron_2->SetTextSize(17);
  Text_Electron_2->SetNDC(1);
  Text_Electron_2->Draw();
  TText *Text_CCP_2 = new TText(0.40, 0.75, "Charge confused protons");
  Text_CCP_2->SetTextAlign(22);
  Text_CCP_2->SetTextColor(kRed);
  Text_CCP_2->SetTextFont(43);
  Text_CCP_2->SetTextSize(17);
  Text_CCP_2->SetNDC(1);
  Text_CCP_2->Draw();
  TLine *line_Pbar_2 = new TLine(0.88, 0.58, 0.80, 0.750);
  line_Pbar_2->SetLineColor(kBlack);
  line_Pbar_2->SetLineWidth(2);
  line_Pbar_2->SetNDC(1);
  line_Pbar_2->Draw();
  TLine *line_Electron_2 = new TLine(0.88, 0.37, 0.80, 0.25);
  line_Electron_2->SetLineColor(kBlack);
  line_Electron_2->SetLineWidth(2);
  line_Electron_2->SetNDC(1);
  line_Electron_2->Draw();
  TLine *line_CCP_2 = new TLine(0.20, 0.56, 0.35, 0.70);
  line_CCP_2->SetLineColor(kBlack);
  line_CCP_2->SetLineWidth(2);
  line_CCP_2->SetNDC(1);
  line_CCP_2->Draw();


  graphicsOutput->cd(3); //fit-data difference
  gPad->SetLogz();
  //gPad->SetGrid();
  TH2F* differenceCopy = (TH2F*) fInputData->Clone("difference_hist");
  differenceCopy->SetTitle("(data-fit)**2/sigma_data**2");
  differenceCopy->GetXaxis()->SetTitle("#bf{#Lambda_{CC}}");
  differenceCopy->GetYaxis()->SetTitle("#bf{#Lambda_{TRD}}");
  //---> This assumes the template has no error!!
  for (int k = 1; k <= differenceCopy->GetXaxis()->GetNbins(); k++)
    for (int l = 1; l <= differenceCopy->GetYaxis()->GetNbins(); l++) {
      float dat = differenceCopy->GetBinContent(k, l);
      float fitres = fFitResult->GetBinContent(k, l);
      double upE, lowE;
      PoissonUncertainty((int) dat, lowE, upE);
      float error = (dat > fitres ? lowE : upE);

      float reldiffabs = (dat + fitres > 0 ? std::pow((dat - fitres), 2) / std::pow(error, 2) : 0);
      differenceCopy->SetBinContent(k, l, reldiffabs);
    }
  differenceCopy->DrawCopy("COLZ");
  differenceCopy->SetStats(0);


  graphicsOutput->cd(4); //legend
  TLegend* legend = new TLegend(0.2, 0.2, 0.8, 0.8, NULL, "brNDC");
  legend->SetFillColor(kWhite);
  //legend->SetLineColor(kWhite);
  legend->AddEntry(dataCopy, Form("Data distribution (%i entries)", (int) fInputData->GetEntries()));
  legend->AddEntry(fitResultCopy, Form("Template fit result (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  legend->Draw();


  return graphicsOutput;
}

TH2F* TemplateFitter2D::CloneCombinedResultHistogram(const char* hname = "") const {

  if (!fHasBeenFit) {
    WARN_OUT << "No fit was done!" << std::endl;
    return 0;
  }

  return (TH2F*) fFitResult->Clone(hname);
}

double TemplateFitter2D::Fit(int i, double lowBound, double highBound) {

  if (lowBound != 0 || highBound != 0)
    WARN_OUT << "TemplateFitter2D does not support fit ranges." << std::endl;

  return TemplateFitter::Fit(i, 0, 0);
}

bool TemplateFitter2D::TestHistogramCompatibility() const {

  //INFO_OUT << "Compatibility..." << std::endl;

  for (int i = 0; i < (int) fInputTemplates.size(); i++) {
    //check minimum, maximum, bin number and bin width
    assert(fInputTemplates.at(i));
    assert(fInputData);
    if (fInputTemplates.at(i)->GetNbinsX() != fInputData->GetNbinsX())
      return false;

    if (fInputTemplates.at(i)->GetNbinsY() != fInputData->GetNbinsY())
      return false;

    if (fInputTemplates.at(i)->GetEntries() == 0)
      return false;

    for (int k = 1; k <= fInputData->GetNbinsX(); k++)
      if (fInputTemplates.at(i)->GetXaxis()->GetBinLowEdge(k) != fInputData->GetXaxis()->GetBinLowEdge(k)) {
        return false;
      }

    for (int k = 1; k <= fInputData->GetNbinsY(); k++)
      if (fInputTemplates.at(i)->GetYaxis()->GetBinLowEdge(k) != fInputData->GetYaxis()->GetBinLowEdge(k)) {
        return false;
      }
  }

  return true;
}

TH1F* TemplateFitter2D::GetResultHistogramXprojection(const char* hname = "", double Ylow, double Yhigh) const {

  if (!fHasBeenFit) {
    WARN_OUT << "No fit was done!" << std::endl;
    return 0;
  }

  if (Ylow == 0 && Yhigh == 0)
    return (TH1F*) Cast(fFitResult)->ProjectionX(hname, 0, -1, "e");
  else
    return (TH1F*) Cast(fFitResult)->ProjectionX(hname, fFitResult->GetYaxis()->FindFixBin(Ylow), fFitResult->GetYaxis()->FindFixBin(Yhigh), "e");
}

TH1F* TemplateFitter2D::GetResultHistogramYprojection(const char* hname = "", double Xlow, double Xhigh) const {

  if (!fHasBeenFit) {
    WARN_OUT << "No fit was done!" << std::endl;
    return 0;
  }

  if (Xlow == 0 && Xhigh == 0)
    return (TH1F*) Cast(fFitResult)->ProjectionY(hname, 0, -1, "e");
  else
    return (TH1F*) Cast(fFitResult)->ProjectionY(hname, fFitResult->GetXaxis()->FindFixBin(Xlow), fFitResult->GetXaxis()->FindFixBin(Xhigh), "e");
}

TCanvas* TemplateFitter2D::CreateResultDrawingXprojection(std::string canvasname, int width, int height, double Ylow, double Yhigh, int BinRemovedNumber_FromLeft, int RebinNumber_X) const {

  if (!fHasBeenFit) {
    WARN_OUT << "No fit was done!" << std::endl;
    return 0;
  }

  bool allrange = false;
  if (Ylow == 0 && Yhigh == 0)
    allrange = true;

  gStyle->SetOptStat(0);
  TCanvas* graphicsOutput = new TCanvas(canvasname.c_str(), canvasname.c_str(), width, height);
  graphicsOutput->cd();

  //TLegend* legend = new TLegend(0.6034913, 0.7009346, 0.8944306, 0.8913551, NULL, "brNDC");
  TLegend* legend = new TLegend(0.20, 0.55, 0.43, 0.85, NULL, "brNDC");
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kBlack);

  TH1* dataCopy_FullRange;
  if (allrange)
    dataCopy_FullRange = Cast(fInputData)->ProjectionX()->DrawCopy("PE");
  else
    dataCopy_FullRange = Cast(fInputData)->ProjectionX("", fInputData->GetYaxis()->FindFixBin(Ylow), fInputData->GetYaxis()->FindFixBin(Yhigh))->DrawCopy("PE");
  int newBinNumber = dataCopy_FullRange->GetNbinsX() - BinRemovedNumber_FromLeft;
  int oldBinNumber = dataCopy_FullRange->GetNbinsX();
  double TotalEntriesNumber = 0;
  TH1* dataCopy = new TH1F("", "", newBinNumber, dataCopy_FullRange->GetBinLowEdge(1+BinRemovedNumber_FromLeft), dataCopy_FullRange->GetBinLowEdge(oldBinNumber) + dataCopy_FullRange->GetBinWidth(oldBinNumber));
  for (int i = 1; i < newBinNumber+1; i++) {
      dataCopy->SetBinContent(i, dataCopy_FullRange->GetBinContent(i+BinRemovedNumber_FromLeft));
      dataCopy->SetBinError(  i, dataCopy_FullRange->GetBinError  (i+BinRemovedNumber_FromLeft));
      TotalEntriesNumber = TotalEntriesNumber + dataCopy_FullRange->GetBinContent(i+BinRemovedNumber_FromLeft);
  }
  dataCopy->Rebin(RebinNumber_X);
  graphicsOutput->GetListOfPrimitives()->Remove(dataCopy_FullRange);
  graphicsOutput->Modified();
  dataCopy->Draw("PE");
  dataCopy->SetStats(0); // disable statistic boxes in template canvas
  dataCopy->SetTitle("");
  const int dataColor = kBlack;
  const int dataMarker = 15;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerColor(dataColor);
  dataCopy->SetMarkerStyle(dataMarker);
  dataCopy->SetMarkerSize(2.0);
  if (allrange){
    legend->AddEntry(dataCopy, Form("Data (%.0f)", fInputData->GetEntries()));
  }
  else{
    legend->AddEntry(dataCopy, Form("Data (%.0f)", TotalEntriesNumber));
  }

  TH1* fitResultCopy;
  if (allrange)
    fitResultCopy = Cast(fFitResult)->ProjectionX()->DrawCopy("hist.same");
  else
    fitResultCopy = Cast(fFitResult)->ProjectionX("", fFitResult->GetYaxis()->FindFixBin(Ylow), fFitResult->GetYaxis()->FindFixBin(Yhigh))->DrawCopy("hist.same");
  fitResultCopy->SetStats(0);
  const int fitColor = kBlack;
  const int fitStyle = 2;
  fitResultCopy->SetLineWidth(4);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);
  /* Disable to show chi2 in total fit.
  if (allrange){
    legend->AddEntry(fitResultCopy, Form("Fit Result (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));}
  else{
    legend->AddEntry(fitResultCopy, Form("Fit Result  (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));}
  */
  fitResultCopy->Rebin(RebinNumber_X);

  for (unsigned int i = 0; i < fResultHistos.size(); ++i) {
    TH1* resultCopy;
    if (allrange){
      resultCopy = Cast(fResultHistos.at(i))->ProjectionX()->DrawCopy("hist.same");
    }
    else{
      resultCopy = Cast(fResultHistos.at(i))->ProjectionX("", fResultHistos.at(i)->GetYaxis()->FindFixBin(Ylow), fResultHistos.at(i)->GetYaxis()->FindFixBin(Yhigh))->DrawCopy("hist.same");
      resultCopy->SetTitle(Cast(fResultHistos.at(i))->GetTitle());
    }
    double ResultCopyNumber = 0;
    for (int k = 1; k < resultCopy->GetNbinsX()+1 ; k++) {
        ResultCopyNumber = ResultCopyNumber + resultCopy->GetBinContent(k);
    }
    std::cout<< "ResultCopyNumber:" << ResultCopyNumber << std::endl;
    resultCopy->SetStats(0);
    resultCopy->Rebin(RebinNumber_X);
    resultCopy->SetLineWidth(2);
    const int color = TemplateColor(i);
    const int fillStyle = TemplateFillStyle(i);
    resultCopy->SetFillStyle(fillStyle);
    resultCopy->SetFillColor(color);
    resultCopy->SetLineColor(color);
    resultCopy->SetMarkerColor(color);
    QString title = resultCopy->GetTitle();
    QString fittedValue = QString("(%2 #pm %3)")
                            //.arg(GetAbsoluteResult().at(i)     , 0, 'f', 1, '0')
                            //.arg(GetAbsoluteResultError().at(i), 0, 'f', 1, '0');
                            //.arg(resultCopy->GetEntries(), 0, 'f', 1, '0')
                            //.arg(GetAbsoluteResultError().at(i)/GetAbsoluteResult().at(i)*resultCopy->GetEntries(), 0, 'f', 1, '0');
                            .arg(ResultCopyNumber, 0, 'f', 1, '0')
                            .arg(GetAbsoluteResultError().at(i)/GetAbsoluteResult().at(i)*ResultCopyNumber, 0, 'f', 1, '0');
                            //.arg(GetRelativeResult().at(i)*TotalEntriesNumber   , 0, 'f', 1, '0')
                            //.arg(GetRelativeResultError().at(i)*TotalEntriesNumber, 0, 'f', 1, '0');
    /*
    std::cout<< "GetRelativeResult:"         << GetRelativeResult().at(i)*TotalEntriesNumber << std::endl;
    std::cout<< "fResultHistos.at(i):"       << fResultHistos.at(i)->GetEntries() << std::endl;
    std::cout<< "GetAbsoluteResult().at(i):" << GetAbsoluteResult().at(i)         << std::endl;
    std::cout<< "TotalEntriesNumber:"        << TotalEntriesNumber                << std::endl;
    std::cout<< "fInputData->GetEntries():"  << fInputData->GetEntries()          << std::endl;
    */

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
  fitResultCopy->SetMinimum(0);
  dataCopy     ->SetMinimum(0);
  fitResultCopy->Draw("hist.same");
  dataCopy->Draw("P E1 same");
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  dataCopy->SetLabelFont(62,"xy");
  dataCopy->SetLabelSize(0.045,"xy");
  dataCopy->GetXaxis()->SetTitle("#Lambda_{CC}");
  dataCopy->GetYaxis()->SetTitle("N");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.045);
  dataCopy->GetYaxis()->SetTitleSize(0.045);

  // Calculate Chi2
  // Option1: Calculate meself.
  double chi2 = 0;
  double tem  = 0;
  for (unsigned int i = 0; i < dataCopy->GetNbinsX(); i++) {
      tem  = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / fitResultCopy->GetBinContent(i+1);
      chi2 = chi2 + tem;
  }
  // Option2:Taken from ROOT.  
  double chi2_v2 = dataCopy->Chi2Test(fitResultCopy, "UU CHI2");
  int ndf = dataCopy->GetNbinsX() - fResultHistos.size();
  legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", chi2, ndf, chi2/ndf));

  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->SetBorderSize(0);
  legend->Draw();
  /*
  TText *Text_Pbar = new TText(0.75, 0.60, "Antiprotons");
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(30);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();

  TText *Text_CCProton = new TText(0.55, 0.35, "Charge confused protons");
  Text_CCProton->SetTextAlign(22);
  Text_CCProton->SetTextColor(kRed);
  Text_CCProton->SetTextFont(43);
  Text_CCProton->SetTextSize(30);
  Text_CCProton->SetNDC(1);
  Text_CCProton->Draw();

  TText *Text_Electron = new TText(0.75, 0.45, "Electrons");
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kGreen);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(30);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();
  */
  return graphicsOutput;
}

TCanvas* TemplateFitter2D::CreateResultDrawingYprojection(std::string canvasname, int width, int height, double Xlow, double Xhigh, int BinRemovedNumber_FromRight, int RebinNumber_Y) const {

  if (!fHasBeenFit) {
    WARN_OUT << "No fit was done!" << std::endl;
    return 0;
  }

  bool allrange = false;
  if (Xlow == 0 && Xhigh == 0)
    allrange = true;

  gStyle->SetOptStat(0);
  TCanvas* graphicsOutput = new TCanvas(canvasname.c_str(), canvasname.c_str(), width, height);
  graphicsOutput->cd();

  //TLegend* legend = new TLegend(0.6034913, 0.7009346, 0.8944306, 0.8913551, NULL, "brNDC");
  TLegend* legend = new TLegend(0.20, 0.55, 0.43, 0.85, NULL, "brNDC");
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kBlack);

  TH1* dataCopy_FullRange;
  if (allrange)
    dataCopy_FullRange = Cast(fInputData)->ProjectionY()->DrawCopy("PE");
  else
    dataCopy_FullRange = Cast(fInputData)->ProjectionY("", fInputData->GetXaxis()->FindFixBin(Xlow), fInputData->GetXaxis()->FindFixBin(Xhigh))->DrawCopy("PE");
  int newBinNumber = dataCopy_FullRange->GetNbinsX() - BinRemovedNumber_FromRight;
  double TotalEntriesNumber = 0;
  TH1* dataCopy = new TH1F("", "", newBinNumber, dataCopy_FullRange->GetBinLowEdge(1), dataCopy_FullRange->GetBinLowEdge(newBinNumber) + dataCopy_FullRange->GetBinWidth(newBinNumber));
  for (int i = 1; i < newBinNumber+1; i++) {
      dataCopy->SetBinContent(i, dataCopy_FullRange->GetBinContent(i));
      dataCopy->SetBinError(  i, dataCopy_FullRange->GetBinError(i));
      TotalEntriesNumber = TotalEntriesNumber + dataCopy_FullRange->GetBinContent(i);
  }
  dataCopy->Rebin(RebinNumber_Y);
  graphicsOutput->GetListOfPrimitives()->Remove(dataCopy_FullRange);
  graphicsOutput->Modified();
  dataCopy->Draw("PE");
  dataCopy->SetStats(0); // disable statistic boxes in template canvas
  dataCopy->SetTitle("");
  const int dataColor = kBlack;
  const int dataMarker = 15;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerColor(dataColor);
  dataCopy->SetMarkerStyle(dataMarker);
  dataCopy->SetMarkerSize(2.0);
  if (allrange){
    legend->AddEntry(dataCopy, Form("Data (%.0f)", fInputData->GetEntries()));
  }
  else{
    legend->AddEntry(dataCopy, Form("Data (%.0f)", TotalEntriesNumber));
  }

  TH1* fitResultCopy;
  if (allrange)
    fitResultCopy = Cast(fFitResult)->ProjectionY()->DrawCopy("hist.same");
  else
    fitResultCopy = Cast(fFitResult)->ProjectionY("", fFitResult->GetXaxis()->FindFixBin(Xlow), fFitResult->GetXaxis()->FindFixBin(Xhigh))->DrawCopy("hist.same");
  fitResultCopy->SetStats(0);
  const int fitColor = kBlack;
  const int fitStyle = 2;
  fitResultCopy->SetLineWidth(4);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);
  /* Disable to show chi2 in total fit.
  if (allrange){
    legend->AddEntry(fitResultCopy, Form("Fit Result (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));}
  else{
    legend->AddEntry(fitResultCopy, Form("Fit Result  (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));}
  */
  fitResultCopy->Rebin(RebinNumber_Y);

  for (unsigned int i = 0; i < fResultHistos.size(); ++i) {
    TH1* resultCopy;
    if (allrange){
      resultCopy = Cast(fResultHistos.at(i))->ProjectionY()->DrawCopy("hist.same");
    }
    else{
      resultCopy = Cast(fResultHistos.at(i))->ProjectionY("", fResultHistos.at(i)->GetXaxis()->FindFixBin(Xlow), fResultHistos.at(i)->GetXaxis()->FindFixBin(Xhigh))->DrawCopy("hist.same");
      resultCopy->SetTitle(Cast(fResultHistos.at(i))->GetTitle()); 
    }
    double ResultCopyNumber = 0;
    for (int k = 1; k < resultCopy->GetNbinsX()+1 ; k++) {
        ResultCopyNumber = ResultCopyNumber + resultCopy->GetBinContent(k);
    }
    resultCopy->Rebin(RebinNumber_Y);
    resultCopy->SetStats(0);
    resultCopy->SetLineWidth(2);
    const int color = TemplateColor(i);
    const int fillStyle = TemplateFillStyle(i);
    resultCopy->SetFillStyle(fillStyle);
    resultCopy->SetFillColor(color);
    resultCopy->SetLineColor(color);
    resultCopy->SetMarkerColor(color);
    resultCopy->SetLineWidth(3.0);
    QString title = resultCopy->GetTitle();
    QString fittedValue = QString("(%2 #pm %3)")
                            //.arg(GetAbsoluteResult().at(i), 0, 'f', 1, '0')
                            //.arg(GetAbsoluteResultError().at(i), 0, 'f', 1, '0');
                            //.arg(resultCopy->GetEntries(), 0, 'f', 1, '0')
                            //.arg(GetAbsoluteResultError().at(i)/GetAbsoluteResult().at(i)*resultCopy->GetEntries(), 0, 'f', 1, '0');
                            .arg(ResultCopyNumber, 0, 'f', 1, '0')
                            .arg(GetAbsoluteResultError().at(i)/GetAbsoluteResult().at(i)*ResultCopyNumber, 0, 'f', 1, '0');
                            //.arg(GetRelativeResult().at(i)*TotalEntriesNumber   , 0, 'f', 1, '0')
                            //.arg(GetRelativeResultError().at(i)*TotalEntriesNumber, 0, 'f', 1, '0');

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
  dataCopy->Draw("P E1 same");
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  dataCopy->SetLabelFont(62,"xy");
  dataCopy->SetLabelSize(0.045,"xy");
  dataCopy->GetXaxis()->SetTitle("#Lambda_{TRD}");
  dataCopy->GetYaxis()->SetTitle("N");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.045);
  dataCopy->GetYaxis()->SetTitleSize(0.045);

  // Calculate Chi2
  //Option1: Calculate meself.
  double chi2 = 0;
  double tem  = 0;
  for (unsigned int i = 0; i < dataCopy->GetNbinsX(); i++) {
      if (fitResultCopy->GetBinContent(i+1)==0) {
          tem = 0;
      }
      else{
          tem  = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / fitResultCopy->GetBinContent(i+1);
      }
      chi2 = chi2 + tem;
  }
  // Option2:Taken from ROOT.
  double chi2_v2 = dataCopy->Chi2Test(fitResultCopy, "UU CHI2");
  int ndf = dataCopy->GetNbinsX() - fResultHistos.size();
  legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", chi2, ndf, chi2/ndf));    

  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->SetBorderSize(0);
  legend->Draw();

  /*
  TText *Text_Pbar = new TText(0.8, 0.30, "Antiprotons");
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(30);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();

  TText *Text_CCProton = new TText(0.39, 0.70, "Charge confused protons");
  Text_CCProton->SetTextAlign(22);
  Text_CCProton->SetTextColor(kRed);
  Text_CCProton->SetTextFont(43);
  Text_CCProton->SetTextSize(30);
  Text_CCProton->SetNDC(1);
  Text_CCProton->Draw();

  TText *Text_Electron = new TText(0.37, 0.44, "Electrons");
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kGreen);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(30);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();
  */

  return graphicsOutput;
}

const TH2* TemplateFitter2D::Cast(const TH1* histogram) const {
  const TH2* castHistogram = dynamic_cast<const TH2*>(histogram);
  assert(castHistogram);
  return castHistogram;
}

} // namespace MYUtilities
