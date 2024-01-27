#include "CrystalBallParameterization.hh"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TStyle.h>

#include <Utilities.hh>

#define INFO_OUT_TAG "CrystalBallParameterization"
#include <debugging.hh>

CrystalBallParameterization::CrystalBallParameterization(const std::string& paramFilename) :
  MigrationParameterization() {

  // input file format:
  // E/GeV   alpha   N   amplitude
  fCrystalBallAlpha = new TGraph(paramFilename.c_str(), "%lg %lg %*lg %*lg");
  fCrystalBallN = new TGraph(paramFilename.c_str(), "%lg %*lg %lg %*lg");
  fCrystalBallAmp = new TGraph(paramFilename.c_str(), "%lg %*lg %*lg %lg");

}

CrystalBallParameterization::~CrystalBallParameterization() {

  delete fCrystalBallAlpha;
  delete fCrystalBallAmp;
  delete fCrystalBallN;
}

TF1* CrystalBallParameterization::MakeFunction(double E) const {

  // careful: these projections are not pdfs, because they are not divided by bin width, yet we use a "pdf" function to parameterize them.

  TF1* fc = new TF1(Form("fc_%g", E), "[4]*ROOT::Math::crystalball_pdf(x, [0], [1], [2], [3])", 0.5, 1.e3);
  fc->SetParNames("#alpha", "n", "sigma", "#mu", "amp");
  double sigma = GaussianEnergyResolution(E);
  fc->SetNpx(1000);
  fc->SetParameter(0, Alpha(E));
  fc->SetParameter(1, N(E));
  fc->SetParameter(4, Amp(E));
  fc->SetParLimits(0, 1.2, 3.0);
  fc->SetParLimits(1, 1.0, 15.0);
  fc->SetParLimits(4, 0.01, 100.);
  fc->FixParameter(2, sigma);
  fc->FixParameter(3, E);

  return fc;
}

TCanvas* CrystalBallParameterization::DrawParameterCanvas() {

  TCanvas* crystalBallParamCanvas = new TCanvas("crystalBallParamCanvas", "crystal ball", gStyle->GetCanvasDefW(), gStyle->GetCanvasDefH());
  crystalBallParamCanvas->Divide(2,2);

  double xmin = 0.5;
  double xmax = 1.e3;

  crystalBallParamCanvas->cd(1);
  Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball #alpha", xmin, xmax, 1.2, 3.2, true);
  fCrystalBallAlpha->Draw("P");
  gPad->SetLogx();

  crystalBallParamCanvas->cd(2);
  Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball N", xmin, xmax, 0.0, 13.0, true);
  fCrystalBallN->Draw("P");
  gPad->SetLogx();

  crystalBallParamCanvas->cd(3);
  Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball amp", xmin, xmax, 0.1, 150.0, true);
  fCrystalBallAmp->Draw("P");
  gPad->SetLogx();
  gPad->SetLogy();

  crystalBallParamCanvas->cd(4);
  Utilities::DrawEmptyPlot("", "E (GeV)", "crystal ball #sigma_{E}/E", xmin, xmax, 0.0, 0.15, true);
  TGraph* grSigma = new TGraph;
  for (int ip=0; ip < fCrystalBallAmp->GetN(); ++ip) {
    double E = fCrystalBallAlpha->GetX()[ip];
    grSigma->SetPoint(ip, E, GaussianEnergyResolution(E)/E);
  }
  grSigma->Draw("P");
  gPad->SetLogx();

  return crystalBallParamCanvas;
}

double CrystalBallParameterization::GaussianEnergyResolution(double E) const {

  double constantTerm = fGaussianConstantTerm;
  if (E>100.0)
    constantTerm = fGaussianConstantTerm + (0.1-fGaussianConstantTerm) * std::log(E/100.0) / std::log(430./100.);

  return E * Utilities::QuadraticSum(fGaussianStatisticTerm/std::sqrt(E), constantTerm);
}

double CrystalBallParameterization::Alpha(double E) const {

  return fCrystalBallAlpha->Eval(E);
}

double CrystalBallParameterization::N(double E) const {

  return fCrystalBallN->Eval(E);
}

double CrystalBallParameterization::Amp(double E) const {

  return Utilities::LogarithmicEval(fCrystalBallAmp, E);
}




