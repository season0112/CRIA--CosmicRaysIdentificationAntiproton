
#include "Statistics.hh"
#include "TemplateFitterforLowEnergy2D.hh"

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaletteAxis.h>

#include <cassert>
#include <cmath>

#include <QString>

#define INFO_OUT_TAG "TemplateFitterforLowEnergy2D"
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
  gPad->SetLeftMargin(0.18);
  TH1* dataCopy = fInputData->DrawCopy("COLZ");
  dataCopy->SetStats(0); // disable statistic boxes in template canvas
  dataCopy->GetXaxis()->SetTitle("#Lambda_{TRD}");
  dataCopy->GetYaxis()->SetTitle("#frac{1}{#beta_{TOF}} - #frac{1}{#beta(R,m_{p})}");
  dataCopy->GetZaxis()->SetTitle("N");
  dataCopy->SetTitle("");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetZaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.055);
  dataCopy->GetYaxis()->SetTitleSize(0.045);
  dataCopy->GetZaxis()->SetTitleSize(0.055);
  dataCopy->GetXaxis()->SetLabelFont(62);
  dataCopy->GetYaxis()->SetLabelFont(62);
  dataCopy->GetZaxis()->SetLabelFont(62);
  dataCopy->GetXaxis()->SetLabelSize(0.05);
  dataCopy->GetYaxis()->SetLabelSize(0.05);
  dataCopy->GetZaxis()->SetLabelSize(0.05);
  dataCopy->GetZaxis()->SetTitleOffset(0.4);  

  // Move the palette
  TPaletteAxis *palette_1 = (TPaletteAxis*)dataCopy->GetListOfFunctions()->FindObject("palette");
  palette_1->SetX1NDC(0.905); 
  palette_1->SetX2NDC(0.955); 
  palette_1->SetY1NDC(0.149); 
  palette_1->SetY2NDC(0.903);

  TText *Text_Pbar = new TText(0.52, 0.25, "Antiprotons");
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(17);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();
  
  TText *Text_Electron = new TText(0.30, 0.78, "Electrons");
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kRed);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(17);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();
  
  TText *Text_and = new TText(0.41, 0.78, "and");
  Text_and->SetTextAlign(22);
  Text_and->SetTextColor(kBlack);
  Text_and->SetTextFont(43);
  Text_and->SetTextSize(17);
  Text_and->SetNDC(1);
  Text_and->Draw();

  TText *Text_Sec = new TText(0.54, 0.78, " Secondaries");
  Text_Sec->SetTextAlign(22);
  Text_Sec->SetTextColor(kGreen);
  Text_Sec->SetTextFont(43);
  Text_Sec->SetTextSize(17);
  Text_Sec->SetNDC(1);
  Text_Sec->Draw();

  TLine *line_Pbar = new TLine(0.50, 0.28, 0.50, 0.54);
  line_Pbar->SetLineColor(kBlack);
  line_Pbar->SetLineWidth(2);
  line_Pbar->SetNDC(1);
  line_Pbar->Draw();
  /*
  TLine *line_Electron = new TLine(0.88, 0.37, 0.80, 0.25);
  line_Electron->SetLineColor(kBlack);
  line_Electron->SetLineWidth(2);
  line_Electron->SetNDC(1);
  line_Electron->Draw();
  */
  TLine *line_Sec = new TLine(0.20, 0.56, 0.35, 0.74);
  line_Sec->SetLineColor(kBlack);
  line_Sec->SetLineWidth(2);
  line_Sec->SetNDC(1);
  line_Sec->Draw();

  graphicsOutput->cd(2); //fit result
  gPad->SetLogz();
  //gPad->SetGrid();
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.18);
  fFitResult->SetStats(0);
  TH1* fitResultCopy = fFitResult->DrawCopy("COLZ");
  fitResultCopy->SetTitle("");
  fitResultCopy->GetXaxis()->SetTitle("#Lambda_{TRD}");
  fitResultCopy->GetYaxis()->SetTitle("#frac{1}{#beta_{TOF}} - #frac{1}{#beta(R,m_{p})}");
  fitResultCopy->GetZaxis()->SetTitle("N");
  fitResultCopy->GetXaxis()->SetTitleFont(62);
  fitResultCopy->GetYaxis()->SetTitleFont(62);
  fitResultCopy->GetZaxis()->SetTitleFont(62);
  fitResultCopy->GetXaxis()->SetTitleSize(0.055);
  fitResultCopy->GetYaxis()->SetTitleSize(0.045);
  fitResultCopy->GetZaxis()->SetTitleSize(0.055);
  fitResultCopy->GetXaxis()->SetLabelFont(62);
  fitResultCopy->GetYaxis()->SetLabelFont(62);
  fitResultCopy->GetZaxis()->SetLabelFont(62);
  fitResultCopy->GetXaxis()->SetLabelSize(0.05);
  fitResultCopy->GetYaxis()->SetLabelSize(0.05);
  fitResultCopy->GetZaxis()->SetLabelSize(0.05);
  fitResultCopy->GetZaxis()->SetTitleOffset(0.4);
  /*
  // Move the palette FIXME:palette_2 doesn't exist.
  TPaletteAxis *palette_2 = (TPaletteAxis*)fitResultCopy->GetListOfFunctions()->FindObject("palette");
  palette_2->SetX1NDC(0.903);
  palette_2->SetX2NDC(0.953);
  palette_2->SetY1NDC(0.15);
  palette_2->SetY2NDC(0.905);
  */

  TText *Text_Pbar_2 = new TText(0.52, 0.25, "Antiprotons");
  Text_Pbar_2->SetTextAlign(22);
  Text_Pbar_2->SetTextColor(kBlue);
  Text_Pbar_2->SetTextFont(43);
  Text_Pbar_2->SetTextSize(17);
  Text_Pbar_2->SetNDC(1);
  Text_Pbar_2->Draw();

  TText *Text_Electron_2 = new TText(0.30, 0.78, "Electrons");
  Text_Electron_2->SetTextAlign(22);
  Text_Electron_2->SetTextColor(kRed);
  Text_Electron_2->SetTextFont(43);
  Text_Electron_2->SetTextSize(17);
  Text_Electron_2->SetNDC(1);
  Text_Electron_2->Draw();

  TText *Text_and_2 = new TText(0.41, 0.78, "and");
  Text_and_2->SetTextAlign(22);
  Text_and_2->SetTextColor(kBlack);
  Text_and_2->SetTextFont(43);
  Text_and_2->SetTextSize(17);
  Text_and_2->SetNDC(1);
  Text_and_2->Draw();

  TText *Text_Sec_2 = new TText(0.54, 0.78, " Secondaries");
  Text_Sec_2->SetTextAlign(22);
  Text_Sec_2->SetTextColor(kGreen);
  Text_Sec_2->SetTextFont(43);
  Text_Sec_2->SetTextSize(17);
  Text_Sec_2->SetNDC(1);
  Text_Sec_2->Draw();

  TLine *line_Pbar_2 = new TLine(0.50, 0.28, 0.50, 0.54);
  line_Pbar_2->SetLineColor(kBlack);
  line_Pbar_2->SetLineWidth(2);
  line_Pbar_2->SetNDC(1);
  line_Pbar_2->Draw();

  TLine *line_Sec_2 = new TLine(0.20, 0.56, 0.35, 0.74);
  line_Sec_2->SetLineColor(kBlack);
  line_Sec_2->SetLineWidth(2);
  line_Sec_2->SetNDC(1);
  line_Sec_2->Draw();


  graphicsOutput->cd(3); //fit-data difference
  gPad->SetLogz();
  gPad->SetGrid();
  TH2F* differenceCopy = (TH2F*) fInputData->Clone("difference_hist");
  differenceCopy->SetTitle("(data-fit)**2/sigma_data**2");
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
  //gPad->SetLogy();
  //gPad->SetGrid();

  //TLegend* legend = new TLegend(0.6034913, 0.7009346, 0.8944306, 0.8913551, NULL, "brNDC");
  TLegend* legend = new TLegend(0.42, 0.76, 0.88, 0.87, NULL, "brNDC");
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
  const int dataColor = kBlack;
  const int dataMarker = 20;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerColor(dataColor);
  dataCopy->SetMarkerStyle(dataMarker);
  dataCopy->SetMarkerSize(2.0);
  dataCopy->SetTitle("");
  if (allrange){
    legend->AddEntry(dataCopy, Form("Data (%.0f entries)", fInputData->GetEntries()));
  }
  else{
    //legend->AddEntry(dataCopy, Form("Data (%.2f < Y < %.2f) (%.0f entries)", Ylow, Yhigh, dataCopy->GetEntries()));
    legend->AddEntry(dataCopy, Form("Data (%.0f)", TotalEntriesNumber)); //dataCopy->GetEntries()
  }

  TH1* fitResultCopy;
  if (allrange)
    fitResultCopy = Cast(fFitResult)->ProjectionX()->DrawCopy("hist.same");
  else
    fitResultCopy = Cast(fFitResult)->ProjectionX("", fFitResult->GetYaxis()->FindFixBin(Ylow), fFitResult->GetYaxis()->FindFixBin(Yhigh))->DrawCopy("hist.same");
  fitResultCopy->SetStats(0);
  const int fitColor = kBlack;
  const int fitStyle = 2;
  fitResultCopy->SetLineWidth(2);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);
  /* Disable to show chi2 in total fit.
  if (allrange)
    //legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Chi2(), NDF(), Chi2() / NDF()));
    legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  else
    //legend->AddEntry(fitResultCopy, Form("Fit (%.2f < Y < %.2f) (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Ylow, Yhigh, Chi2(), NDF(), Chi2() / NDF()));
    legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  */
  fitResultCopy->Rebin(RebinNumber_X);
  
  for (unsigned int i = 0; i < fResultHistos.size(); ++i) {
    TH1* resultCopy;
    if (allrange){
      resultCopy = Cast(fResultHistos.at(i))->ProjectionX()->DrawCopy("hist.same");}
    else{
      resultCopy = Cast(fResultHistos.at(i))->ProjectionX("", fResultHistos.at(i)->GetYaxis()->FindFixBin(Ylow), fResultHistos.at(i)->GetYaxis()->FindFixBin(Yhigh))->DrawCopy("hist.same");
      resultCopy->SetTitle(Cast(fResultHistos.at(i))->GetTitle());}
    resultCopy->SetStats(0);
    resultCopy->Rebin(RebinNumber_X);
    resultCopy->SetLineWidth(2);
    const int color = TemplateColor(i);
    const int fillStyle = TemplateFillStyle(i);
    resultCopy->SetFillStyle(fillStyle);
    resultCopy->SetFillColor(color);
    resultCopy->SetLineColor(color);
    resultCopy->SetMarkerColor(color);
    /*
    // split to long titles in two lines
    QString title = resultCopy->GetTitle();
    //QString fittedValue = QString(" #color[%1]{(%2 #pm %3)}")
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
  //dataCopy->GetXaxis()->SetTitleOffset(1.7);
  dataCopy->GetXaxis()->SetTitle("#Lambda_{TRD}");
  dataCopy->GetYaxis()->SetTitle("N");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.05);
  dataCopy->GetYaxis()->SetTitleSize(0.04);


  // Calculate Chi2
  /*
  //Option1: Calculate mesulf.
  double chi2 = 0;
  double tem  = 0;
  for (unsigned int i = 0; i < dataCopy->GetNbinsX(); i++) {
      tem  = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / fitResultCopy->GetBinContent(i+1);
      //tem = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / pow(dataCopy->GetBinError(i+1), 2);  
      //tem = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / pow( pow(dataCopy->GetBinContent(i+1), 0.5), 2);
      //tem = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / pow( fitResultCopy->GetBinError(i+1) * pow(fitResultCopy->GetBinCenter(i+1),2), 2);
      chi2 = chi2 + tem;
  }
  std::cout<< chi2 << std::endl;
  std::cout<< ndf << std::endl;
  std::cout<< chi2/ndf << std::endl;
  std::cout<< dataCopy->Chi2Test(fitResultCopy, "UU P") << std::endl;
  std::cout<< dataCopy->Chi2Test(fitResultCopy, "UU CHI2") << std::endl; 
  */
  // Option2:Taken from ROOT.
  double chi2 = dataCopy->Chi2Test(fitResultCopy, "UU CHI2");
  int ndf = dataCopy->GetNbinsX() - fResultHistos.size();
  legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", chi2, ndf, chi2/ndf));

  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->Draw();


  TText *Text_Pbar = new TText(0.51, 0.30, "Antiprotons");
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(30);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();

  TText *Text_Electron = new TText(0.3, 0.6, "Electrons");
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kRed);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(30);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();

  TText *Text_Secondary = new TText(0.25, 0.23, "Secondaries");
  Text_Secondary->SetTextAlign(22);
  Text_Secondary->SetTextColor(kGreen);
  Text_Secondary->SetTextFont(43);
  Text_Secondary->SetTextSize(30);
  Text_Secondary->SetNDC(1);
  Text_Secondary->Draw();

  return graphicsOutput;
}

TCanvas* TemplateFitter2D::CreateResultDrawingXprojection_LogY(std::string canvasname, int width, int height, double Ylow, double Yhigh, int BinRemovedNumber_FromLeft, int RebinNumber_X) const {

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
  gPad->SetLogy();

  TLegend* legend = new TLegend(0.58, 0.67, 0.88, 0.87, NULL, "brNDC");
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
  const int dataColor = kBlack;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerStyle(15);
  dataCopy->SetMarkerSize(2.0);
  dataCopy->SetMarkerColor(1);
  dataCopy->SetTitle("");
  if (allrange){
    legend->AddEntry(dataCopy, Form("Data (%.0f entries)", fInputData->GetEntries()));}
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
  fitResultCopy->SetLineWidth(2);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);
  /* Disable to show chi2 in total fit.
  if (allrange){
      //legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Chi2(), NDF(), Chi2() / NDF()));
      legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  }
  else{
      legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  }
  */
  fitResultCopy->Rebin(RebinNumber_X);

  for (unsigned int i = 0; i < fResultHistos.size(); ++i) {
    TH1* resultCopy;
    if (allrange){
      resultCopy = Cast(fResultHistos.at(i))->ProjectionX()->DrawCopy("hist.same");}
    else{
      resultCopy = Cast(fResultHistos.at(i))->ProjectionX("", fResultHistos.at(i)->GetYaxis()->FindFixBin(Ylow), fResultHistos.at(i)->GetYaxis()->FindFixBin(Yhigh))->DrawCopy("hist.same");
      resultCopy->SetTitle(Cast(fResultHistos.at(i))->GetTitle());}
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
                            .arg(resultCopy->GetEntries(), 0, 'f', 1, '0')
                            .arg(GetAbsoluteResultError().at(i)/GetAbsoluteResult().at(i)*resultCopy->GetEntries(), 0, 'f', 1, '0');
    title.append(fittedValue);
    if (title.contains(",")) {
      title.replace(",", "}{");
      title.prepend("#splitline{");
      title.append("}");
    }
    legend->AddEntry(resultCopy, qPrintable(title), "F");
  }

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
  /*
  //Option1: Calculate meself.
  double chi2 = 0;
  double tem  = 0;
  for (unsigned int i = 0; i < dataCopy->GetNbinsX(); i++) {
      tem  = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / fitResultCopy->GetBinContent(i+1);
      chi2 = chi2 + tem;
  }
  */
  // Option2:Taken from ROOT.
  double chi2 = dataCopy->Chi2Test(fitResultCopy, "UU CHI2");
  int ndf = dataCopy->GetNbinsX() - fResultHistos.size();
  legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", chi2, ndf, chi2/ndf));
 
  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->SetBorderSize(0);
  legend->Draw();
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
  //gPad->SetGrid();

  //TLegend* legend = new TLegend(0.1034913, 0.7009346, 0.3944306, 0.8913551, NULL, "brNDC");
  TLegend* legend = new TLegend(0.42, 0.76, 0.88, 0.87, NULL, "brNDC");
  //TLegend* legend = new TLegend(0.18, 0.76, 0.44, 0.87, NULL, "brNDC");
  legend->SetFillColor(kWhite);
  //legend->SetLineColor(kWhite);
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
  const int dataColor = kBlack;
  const int dataMarker = 20;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerColor(dataColor);
  dataCopy->SetMarkerStyle(dataMarker);
  dataCopy->SetMarkerSize(2.0);
  dataCopy->SetTitle("");
  if (allrange){
    //legend->AddEntry(dataCopy, Form("Data (%.0f)", fInputData->GetEntries()));
    legend->AddEntry(dataCopy, Form("Data (%.0f entries)", fInputData->GetEntries()));
  }
  else{
    //legend->AddEntry(dataCopy, Form("Data (Y projection, %.2f < X < %.2f) (%.0f entries)", Xlow, Xhigh, fInputData->GetEntries()));
    legend->AddEntry(dataCopy, Form("Data (%.0f)", TotalEntriesNumber ));  ////dataCopy->GetEntries()
  } 
  
  TH1* fitResultCopy;
  if (allrange)
    fitResultCopy = Cast(fFitResult)->ProjectionY()->DrawCopy("hist.same");
  else
    fitResultCopy = Cast(fFitResult)->ProjectionY("", fFitResult->GetXaxis()->FindFixBin(Xlow), fFitResult->GetXaxis()->FindFixBin(Xhigh))->DrawCopy("hist.same");
  fitResultCopy->SetStats(0);
  const int fitColor = kBlack;
  const int fitStyle = 2;
  fitResultCopy->SetLineWidth(2);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);
  /* Disable to show chi2 in total fit.
  if (allrange){
    //legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Chi2(), NDF(), Chi2() / NDF()));
    legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  }
  else{
    //legend->AddEntry(fitResultCopy, Form("Fit (%.2f < X < %.2f) (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Xlow, Xhigh, Chi2(), NDF(), Chi2() / NDF()));
    //legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Chi2(), NDF(), Chi2() / NDF()));
    legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  }
  */
  fitResultCopy->Rebin(RebinNumber_Y);  

  
  for (unsigned int i = 0; i < fResultHistos.size(); ++i) {
    TH1* resultCopy;
    if (allrange){
      resultCopy = Cast(fResultHistos.at(i))->ProjectionY()->DrawCopy("hist.same");}
    else{
      resultCopy = Cast(fResultHistos.at(i))->ProjectionY("", fResultHistos.at(i)->GetXaxis()->FindFixBin(Xlow), fResultHistos.at(i)->GetXaxis()->FindFixBin(Xhigh))->DrawCopy("hist.same");
      resultCopy->SetTitle(Cast(fResultHistos.at(i))->GetTitle());}
    resultCopy->SetStats(0);
    resultCopy->Rebin(RebinNumber_Y);
    resultCopy->SetLineWidth(2);
    const int color = TemplateColor(i);
    const int fillStyle = TemplateFillStyle(i);
    resultCopy->SetFillStyle(fillStyle);
    resultCopy->SetFillColor(color);
    resultCopy->SetLineColor(color);
    resultCopy->SetMarkerColor(color);
    resultCopy->SetLineWidth(3.0);
    /*
    // split to long titles in two lines
    QString title = resultCopy->GetTitle();
    //QString fittedValue = QString(" #color[%1]{(%2 #pm %3)}")
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
  dataCopy->GetXaxis()->SetTitleOffset(1.7);
  dataCopy->GetXaxis()->SetTitle("#frac{1}{#bf{#beta}_{TOF}} - #frac{1}{#bf{#beta}(R,m_{p})}");
  dataCopy->GetYaxis()->SetTitle("N");
  dataCopy->GetXaxis()->SetTitleFont(62);  
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.035);
  dataCopy->GetYaxis()->SetTitleSize(0.04);

  // Calculate Chi2
  /*
  //Option1: Calculate mesulf.
  double chi2 = 0;
  double tem  = 0;
  for (unsigned int i = 0; i < dataCopy->GetNbinsX(); i++) {
      tem  = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / fitResultCopy->GetBinContent(i+1);
      chi2 = chi2 + tem;
  }
  */
  // Option2:Taken from ROOT.
  double chi2 = dataCopy->Chi2Test(fitResultCopy, "UU CHI2");
  int ndf = dataCopy->GetNbinsX() - fResultHistos.size();
  legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", chi2, ndf, chi2/ndf));  
 

  //gStyle->SetLegendFont(62);
  //gStyle->SetLegendTextSize();
  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->Draw();

  TText *Text_Pbar = new TText(0.32, 0.27, "Antiprotons");
  Text_Pbar->SetTextAlign(22);
  Text_Pbar->SetTextColor(kBlue);
  Text_Pbar->SetTextFont(43);
  Text_Pbar->SetTextSize(30);
  Text_Pbar->SetNDC(1);
  Text_Pbar->Draw();

  TText *Text_Electron = new TText(0.28, 0.4, "Electrons");
  Text_Electron->SetTextAlign(22);
  Text_Electron->SetTextColor(kRed);
  Text_Electron->SetTextFont(43);
  Text_Electron->SetTextSize(30);
  Text_Electron->SetNDC(1);
  Text_Electron->Draw();

  TText *Text_Secondary = new TText(0.32, 0.185, "Secondaries");
  Text_Secondary->SetTextAlign(22);
  Text_Secondary->SetTextColor(kGreen);
  Text_Secondary->SetTextFont(43);
  Text_Secondary->SetTextSize(30);
  Text_Secondary->SetNDC(1);
  Text_Secondary->Draw();

  return graphicsOutput;
}

TCanvas* TemplateFitter2D::CreateResultDrawingYprojection_LogY(std::string canvasname, int width, int height, double Xlow, double Xhigh, int BinRemovedNumber_FromRight, int RebinNumber_Y) const {

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
  gPad->SetLogy();

  TLegend* legend = new TLegend(0.2, 0.67, 0.45, 0.87, NULL, "brNDC");
  //TLegend* legend = new TLegend(0.18, 0.76, 0.44, 0.87, NULL, "brNDC");
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
  const int dataColor = kBlack;
  const int dataMarker = 20;
  dataCopy->SetLineWidth(2);
  dataCopy->SetLineColor(dataColor);
  dataCopy->SetMarkerColor(dataColor);
  dataCopy->SetMarkerStyle(dataMarker);
  dataCopy->SetMarkerSize(2.0);
  dataCopy->SetTitle("");
  if (allrange){
    //legend->AddEntry(dataCopy, Form("Data (%.0f)", fInputData->GetEntries()));
    legend->AddEntry(dataCopy, Form("Data (%.0f entries)", fInputData->GetEntries()));
  }
  else{
    //legend->AddEntry(dataCopy, Form("Data (%.0f)", dataCopy->GetEntries()));
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
  fitResultCopy->SetLineWidth(2);
  fitResultCopy->SetLineColor(fitColor);
  fitResultCopy->SetLineStyle(fitStyle);
  /* Disable to show chi2 in total fit.
  if (allrange){
    //legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", Chi2(), NDF(), Chi2() / NDF()));
    legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  }
  else{
    legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f)", Chi2() / NDF()));
  }
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
                            .arg(resultCopy->GetEntries(), 0, 'f', 1, '0')
                            .arg(GetAbsoluteResultError().at(i)/GetAbsoluteResult().at(i)*resultCopy->GetEntries(), 0, 'f', 1, '0');
    title.append(fittedValue);
    if (title.contains(",")) {
      title.replace(",", "}{");
      title.prepend("#splitline{");
      title.append("}");
    }
    legend->AddEntry(resultCopy, qPrintable(title), "F");

  }

  gStyle->SetErrorX(0);
  fitResultCopy->Draw("hist.same");
  dataCopy->Draw("P E1 same");

  dataCopy->SetLabelFont(62,"xy");
  dataCopy->SetLabelSize(0.045,"xy");
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  dataCopy->GetXaxis()->SetTitleOffset(1.7);
  dataCopy->GetXaxis()->SetTitle("#frac{1}{#bf{#beta}_{TOF}} - #frac{1}{#bf{#beta}(R,m_{p})}");
  dataCopy->GetYaxis()->SetTitle("N");
  dataCopy->GetXaxis()->SetTitleFont(62);
  dataCopy->GetYaxis()->SetTitleFont(62);
  dataCopy->GetXaxis()->SetTitleSize(0.035);
  dataCopy->GetYaxis()->SetTitleSize(0.04);
  // To have space for legend
  dataCopy->GetYaxis()->SetRangeUser(1, dataCopy->GetMaximum()*10); 

  // Calculate Chi2
  /*
  //Option1: Calculate myself.
  double chi2 = 0;
  double tem  = 0;
  for (unsigned int i = 0; i < dataCopy->GetNbinsX(); i++) {
      tem  = pow(dataCopy->GetBinContent(i+1) - fitResultCopy->GetBinContent(i+1), 2) / fitResultCopy->GetBinContent(i+1);
      chi2 = chi2 + tem;
  }
  */
  // Option2:Taken from ROOT.
  double chi2 = dataCopy->Chi2Test(fitResultCopy, "UU CHI2");
  int ndf = dataCopy->GetNbinsX() - fResultHistos.size();
  legend->AddEntry(fitResultCopy, Form("Fit (#chi^{2}/ndf=%5.2f/%d=%5.2f)", chi2, ndf, chi2/ndf));    

  legend->SetTextFont(62);
  legend->SetTextSize(0.02);
  legend->SetBorderSize(0);
  legend->Draw();
  return graphicsOutput;
}

const TH2* TemplateFitter2D::Cast(const TH1* histogram) const {
  const TH2* castHistogram = dynamic_cast<const TH2*>(histogram);
  assert(castHistogram);
  return castHistogram;
}

} // namespace MYUtilities
