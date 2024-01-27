#include "CompositeParameterizationFunction.hh"

#include <TF1.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/PdfFuncMathCore.h>

CompositeParameterizationFunction::CompositeParameterizationFunction(double Emin, double Emax) :
  fEmin(Emin),
  fEmax(Emax) {

  fLandauPart = new TF1("fLandauPart", "landau", Emin, Emax);
}

CompositeParameterizationFunction::~CompositeParameterizationFunction() {

  delete fLandauPart;
}

void CompositeParameterizationFunction::SetParameters(double landauMpv, double landauSigma, double cbAlpha, double cbN, double cbMu, double cbSigma, double landauFraction) {

  fLandauPart->SetParameters(1.0, landauMpv, landauSigma);

  fLandauMpv = landauMpv;
  fLandauSigma = landauSigma;
  fCBAlpha = cbAlpha;
  fCBN = cbN;
  fCBMu = cbMu;
  fCBSigma = cbSigma;
  fLandauFraction = landauFraction;
}

void CompositeParameterizationFunction::SetUnderOverflow(double underflowProb, double overflowProb) {

  fUnderflowProb = underflowProb;
  fOverflowProb = overflowProb;
}

double CompositeParameterizationFunction::Evaluate(double* x, double* p) {

  SetParameters(p[0], p[1], p[2], p[3], p[4], p[5], p[6]);

  double E = x[0];

  double landauPart = 0.0;
  if (fLandauFraction > 0.0 && E < fCBMu) {
    double landauIntegral = fLandauPart->Integral(fEmin, fCBMu);
    double landauAmp = fLandauFraction / landauIntegral;
    landauPart = landauAmp * TMath::Landau(E, fLandauMpv, fLandauSigma, kFALSE);
  }

  double cbAmp = 1.0 - fUnderflowProb - fOverflowProb - fLandauFraction;
  if (cbAmp < 0.0)
    cbAmp = 0.0;
  if (cbAmp > 1.0)
    cbAmp = 1.0;

  if (fUnderflowProb > 0.01) {
    double cbIntegral = ROOT::Math::crystalball_cdf(fEmax, fCBAlpha, fCBN, fCBSigma, fCBMu) - ROOT::Math::crystalball_cdf(fEmin, fCBAlpha, fCBN, fCBSigma, fCBMu);
    cbAmp /= cbIntegral;
  }

  return landauPart + cbAmp * ROOT::Math::crystalball_pdf(E, fCBAlpha, fCBN, fCBSigma, fCBMu);
}
