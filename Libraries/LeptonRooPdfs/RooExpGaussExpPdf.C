#include "RooExpGaussExpPdf.hh"

#include <cmath>
#include <cassert>

RooExpGaussExpPdf::RooExpGaussExpPdf(const char *name, const char *title,
			             RooAbsReal& _x,
				     RooAbsReal& _mu, RooAbsReal& _sigma,
				     RooAbsReal& _kLow, RooAbsReal& _kHigh)
  : RooAbsPdf(name, title)
  , x("x", "x", this, _x)
  , mu("mu", "mu", this, _mu)
  , sigma("sigma", "sigma", this, _sigma)
  , kLow("kLow", "kLow", this, _kLow)
  , kHigh("kHigh", "kHigh", this, _kHigh) {

}

RooExpGaussExpPdf::RooExpGaussExpPdf(const RooExpGaussExpPdf& other, const char* name)
  : RooAbsPdf(other, name)
  , x("x", this, other.x)
  , mu("mu", this, other.mu)
  , sigma("sigma", this, other.sigma)
  , kLow("kLow", this, other.kLow)
  , kHigh("kHigh", this, other.kHigh) {

}

Double_t RooExpGaussExpPdf::evaluate() const {

  double xbar = (x - mu) / sigma;
  double normalization = std::exp(3.0 * std::pow(kLow,  2) / 2.0) / kLow
	               + std::exp(3.0 * std::pow(kHigh, 2) / 2.0) / kHigh
		       + std::sqrt(M_PI / 2.0) * (std::erf(kHigh / std::sqrt(2)) - std::erf(kLow / std::sqrt(2)));
  if (xbar <= -kLow)
    return std::exp(std::pow(kLow, 2) / 2.0 + kLow * xbar) / normalization;
  if (xbar > -kLow && xbar <= kHigh) 
    return std::exp(-std::pow(xbar, 2) / 2.0) / normalization;
  if (xbar > kHigh)
    return std::exp((kHigh * kHigh / 2) - kHigh * xbar) / normalization;
  assert(false);
  return 0.0;
}
