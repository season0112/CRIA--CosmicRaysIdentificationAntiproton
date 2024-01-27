#include "RooLangausPdf.hh"

#include <TMath.h>

RooLangausPdf::RooLangausPdf(const char *name, const char *title,
                               RooAbsReal& _x, RooAbsReal& _widthLandau, RooAbsReal& _mpvLandau, RooAbsReal& _totalArea, RooAbsReal& _sigmaGauss)
  : RooAbsPdf(name, title)
  , x("x", "x", this, _x)
  , widthLandau("widthLandau", "widthLandau", this, _widthLandau)
  , mpvLandau("mpvLandau", "mpvLandau", this, _mpvLandau)
  , totalArea("totalArea", "totalArea", this, _totalArea)
  , sigmaGauss("sigmaGauss", "sigmaGauss", this, _sigmaGauss) {

}

RooLangausPdf::RooLangausPdf(const RooLangausPdf& other, const char* name)
  : RooAbsPdf(other, name)
  , x("x", this, other.x)
  , widthLandau("widthLandau", this, other.widthLandau)
  , mpvLandau("mpvLandau", this, other.mpvLandau)
  , totalArea("totalArea", this, other.totalArea)
  , sigmaGauss("sigmaGauss", this, other.sigmaGauss) {

}

Double_t RooLangausPdf::evaluate() const {

  static double sOneOverSqrtTwoPi = 0.3989422804014;
  static double sMPVShift = -0.22278298;
  static double sNumberOfConvolutionSteps = 100.0;
  static double sConvolutionExtensionInGaussianSigmas = 5.0;

  // Convolution integral stepping / range
  double xLowerBound = x - sConvolutionExtensionInGaussianSigmas * sigmaGauss;
  double xUpperBound = x + sConvolutionExtensionInGaussianSigmas * sigmaGauss;
  double step = (xUpperBound - xLowerBound) / sNumberOfConvolutionSteps;
  double mpvShiftCorrected = mpvLandau - sMPVShift * widthLandau;

  // Compute convolution integral via summation.
  double sum = 0.0;
  double xPosition = 0.0;
  for (double i = 1.0; i <= sNumberOfConvolutionSteps / 2.0; ++i) {
    xPosition = xLowerBound + (i - 0.5) * step;
    sum += TMath::Landau(xPosition, mpvShiftCorrected, widthLandau) / widthLandau * TMath::Gaus(x, xPosition, sigmaGauss);

    xPosition = xUpperBound - (i - 0.5) * step;
    sum += TMath::Landau(xPosition, mpvShiftCorrected, widthLandau) / widthLandau * TMath::Gaus(x, xPosition, sigmaGauss);
  }

  return (totalArea * step * sum * sOneOverSqrtTwoPi / sigmaGauss);
}
