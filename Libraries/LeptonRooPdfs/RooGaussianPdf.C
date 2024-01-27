#include "RooGaussianPdf.hh"

#include <RooMath.h>

#include <cmath>

RooGaussianPdf::RooGaussianPdf(const char* name, const char* title, RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _sigma)
  : RooAbsPdf(name, title)
  , x("x", "Observable", this, _x)
  , mean("mean", "Mean", this, _mean)
  , sigma("sigma", "Width", this, _sigma) {

}

RooGaussianPdf::RooGaussianPdf(const RooGaussianPdf& other, const char* name)
  : RooAbsPdf(other, name)
  , x("x", this, other.x)
  , mean("mean", this, other.mean)
  , sigma("sigma", this, other.sigma) {

}


Double_t RooGaussianPdf::evaluate() const {

  return std::exp(-0.5 * std::pow(x - mean, 2) / std::pow(sigma, 2));
}

Int_t RooGaussianPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const {

  if (matchArgs(allVars, analVars, x))
    return 1;
  if (matchArgs(allVars, analVars, mean))
    return 2;

  return 0;
}

Double_t RooGaussianPdf::analyticalIntegral(Int_t code, const char* rangeName) const {

  assert(code == 1 || code == 2);

  static const double sSqrt2 = std::sqrt(2.0) ;
  static const double sSqrtPiBy2 = std::sqrt(std::atan2(0.0, -1.0) / 2.0);

  double xScaled = sSqrt2 * sigma;

  if (code == 1)
    return sSqrtPiBy2 * sigma * (RooMath::erf((x.max(rangeName) - mean) / xScaled) - RooMath::erf((x.min(rangeName) - mean) / xScaled));

  if (code == 2)
    return sSqrtPiBy2 * sigma * (RooMath::erf((mean.max(rangeName) - x) / xScaled) - RooMath::erf((mean.min(rangeName) - x) / xScaled));
  
  return 0.0;
}
