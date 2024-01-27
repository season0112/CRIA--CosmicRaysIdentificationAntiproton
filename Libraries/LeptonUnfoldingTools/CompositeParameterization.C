#include "CompositeParameterization.hh"

#include "CompositeParameterizationFunction.hh"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>

#include <Utilities.hh>

#include <iostream>
#include <fstream>

#define INFO_OUT_TAG "CompositeParameterization"
#include <debugging.hh>

CompositeParameterization::CompositeParameterization(const std::string& paramFilename) :
  MigrationParameterization() {

  grLandauMPV = new TGraph;
  grLandauSigma = new TGraph;
  grCbAlpha = new TGraph;
  grCbN = new TGraph;
  grCbMu = new TGraph;
  grCbRelSigma = new TGraph;
  grLandauFraction = new TGraph;
  grProbUnderflow = new TGraph;
  grProbOverflow = new TGraph;

  std::ifstream inputfile(paramFilename);
  if (!inputfile.good())
    FATAL_OUT << "Error reading parameter file " << paramFilename << std::endl;

  double E;
  double landauMpv;
  double landauSigma;
  double cbAlpha;
  double cbN;
  double cbMu;
  double cbRelSigma;
  double landauFraction;
  double probUnderflow;
  double probOverflow;

  while (inputfile.good()) {

    inputfile >> E >> landauMpv >> landauSigma >> cbAlpha >> cbN >> cbMu >> cbRelSigma >> landauFraction >> probUnderflow >> probOverflow;

    if (inputfile.eof())
      break;

    int ip = grCbAlpha->GetN();
    grLandauMPV->SetPoint(ip, E, landauMpv);
    grLandauSigma->SetPoint(ip, E, landauSigma);
    grCbAlpha->SetPoint(ip, E, cbAlpha);
    grCbN->SetPoint(ip, E, cbN);
    grCbMu->SetPoint(ip, E, cbMu);
    grCbRelSigma->SetPoint(ip, E, cbRelSigma);
    grLandauFraction->SetPoint(ip, E, landauFraction);
    grProbUnderflow->SetPoint(ip, E, probUnderflow);
    grProbOverflow->SetPoint(ip, E, probOverflow);
  }
}

CompositeParameterization::~CompositeParameterization() {

  delete grLandauMPV;
  delete grLandauSigma;
  delete grCbAlpha;
  delete grCbN;
  delete grCbMu;
  delete grCbRelSigma;
  delete grLandauFraction;
  delete grProbUnderflow;
  delete grProbOverflow;
}

TF1* CompositeParameterization::MakeFunction(double E) const {

  CompositeParameterizationFunction* cpf = new CompositeParameterizationFunction(fEmin, fEmax);
  cpf->SetUnderOverflow(grProbUnderflow->Eval(E), grProbOverflow->Eval(E));

  TF1* fc = new TF1("CompositeParameterizationFunc", cpf, &CompositeParameterizationFunction::Evaluate, fEmin, fEmax, 7, "CompositeParameterizationFunction", "Evaluate");
  fc->SetParNames("landauMPV", "landauSigma", "#alpha", "n", "#mu", "#sigma", "landauFraction");
  fc->SetNpx(1000);

  double cbN = grCbN->Eval(E);
  if (cbN < 1.001)
    cbN = 1.001;
  fc->SetParameters(grLandauMPV->Eval(E), grLandauSigma->Eval(E),
                    grCbAlpha->Eval(E), cbN,
                    grCbMu->Eval(E), grCbRelSigma->Eval(E) * E,
                    grLandauFraction->Eval(E));

  for (int ip = 0; ip < fc->GetNpar(); ++ip) {
    DEBUG_OUT << "Par " << ip << " " << fc->GetParName(ip) << " = " << fc->GetParameter(ip) << std::endl;
  }

  return fc;
}

double CompositeParameterization::UnderflowProbability(double E) const {

  return grProbUnderflow->Eval(E);
}

double CompositeParameterization::OverflowProbability(double E) const {

  return grProbOverflow->Eval(E);
}

TGraph* CompositeParameterization::GraphLandauMPV() const {
  return grLandauMPV;
}

TGraph* CompositeParameterization::GraphLandauSigma() const {
  return grLandauSigma;
}

TGraph* CompositeParameterization::GraphCbAlpha() const {
  return grCbAlpha;
}

TGraph* CompositeParameterization::GraphCbN() const {
  return grCbN;
}

TGraph* CompositeParameterization::GraphCbMu() const {
  return grCbMu;
}

TGraph* CompositeParameterization::GraphCbRelSigma() const {
  return grCbRelSigma;
}

TGraph* CompositeParameterization::GraphLandauFraction() const {
  return grLandauFraction;
}

TCanvas* CompositeParameterization::DrawParameterCanvas() {

  TCanvas* compositeParamCanvas = new TCanvas("compositeParamCanvas", "composite", gStyle->GetCanvasDefW(), gStyle->GetCanvasDefH());
  compositeParamCanvas->Divide(4,2);

  double rightMargin = 0.08;

  compositeParamCanvas->cd(1);
  gPad->SetRightMargin(rightMargin);
  TH1* bg1 = Utilities::DrawEmptyPlot("", "E (GeV)", "Landau MPV", fEmin, fEmax, 0.0, 1.5, true);
  bg1->GetXaxis()->SetMoreLogLabels(false);
  bg1->GetYaxis()->SetTitleOffset(1.55);
  grLandauMPV->Draw("PL");
  gPad->SetLogx();

  compositeParamCanvas->cd(2);
  gPad->SetRightMargin(rightMargin);
  TH1* bg2 = Utilities::DrawEmptyPlot("", "E (GeV)", "Landau #sigma", fEmin, fEmax, 0.0, 0.4, true);
  bg2->GetXaxis()->SetMoreLogLabels(false);
  bg2->GetYaxis()->SetTitleOffset(1.55);
  grLandauSigma->Draw("PL");
  gPad->SetLogx();

  compositeParamCanvas->cd(3);
  gPad->SetRightMargin(rightMargin);
  TH1* bg3 = Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball #alpha", fEmin, fEmax, 0.0, 4.0, true);
  bg3->GetXaxis()->SetMoreLogLabels(false);
  bg3->GetYaxis()->SetTitleOffset(1.55);
  grCbAlpha->Draw("PL");
  gPad->SetLogx();

  compositeParamCanvas->cd(4);
  gPad->SetRightMargin(rightMargin);
  TH1* bg4 = Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball N", fEmin, fEmax, 0.0, 8.0, true);
  bg4->GetXaxis()->SetMoreLogLabels(false);
  bg4->GetYaxis()->SetTitleOffset(1.55);
  grCbN->Draw("PL");
  gPad->SetLogx();

  compositeParamCanvas->cd(5);
  gPad->SetRightMargin(rightMargin);
  TH1* bg5 = Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball #mu", fEmin, fEmax, fEmin, fEmax, true);
  bg5->GetXaxis()->SetMoreLogLabels(false);
  bg5->GetYaxis()->SetTitleOffset(1.55);
  grCbMu->Draw("PL");
  gPad->SetLogx();
  gPad->SetLogy();

  compositeParamCanvas->cd(6);
  gPad->SetRightMargin(rightMargin);
  TH1* bg6 = Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball #sigma_{E}/E", fEmin, fEmax, 0.0, 0.2, true);
  bg6->GetXaxis()->SetMoreLogLabels(false);
  bg6->GetYaxis()->SetTitleOffset(1.55);
  grCbRelSigma->Draw("PL");
  gPad->SetLogx();

  compositeParamCanvas->cd(7);
  gPad->SetRightMargin(rightMargin);
  TH1* bg7 = Utilities::DrawEmptyPlot("", "E (GeV)", "Landau fraction", fEmin, fEmax, -0.001, 0.02, true);
  bg7->GetXaxis()->SetMoreLogLabels(false);
  bg7->GetYaxis()->SetTitleOffset(1.55);
  grLandauFraction->Draw("PL");
  gPad->SetLogx();

  compositeParamCanvas->cd(8);
  gPad->SetRightMargin(rightMargin);
  TH1* bg8 = Utilities::DrawEmptyPlot("", "E (GeV)", "prob (under-/overflow)", fEmin, fEmax, 0.0001, 0.3, true);
  bg8->GetXaxis()->SetMoreLogLabels(false);
  bg8->GetYaxis()->SetTitleOffset(1.55);
  grProbOverflow->SetMarkerStyle(kOpenCircle);
  grProbOverflow->Draw("P");
  grProbUnderflow->Draw("PL");
  gPad->SetLogx();
  gPad->SetLogy();

  return compositeParamCanvas;
}





