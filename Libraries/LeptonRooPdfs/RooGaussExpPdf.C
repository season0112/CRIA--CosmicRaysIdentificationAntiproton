#include "RooGaussExpPdf.hh"

#include <cmath>
#include <cassert>

RooGaussExpPdf::RooGaussExpPdf(const char *name, const char *title,
                               RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _k)
  : RooAbsPdf(name, title)
  , x("x", "x", this, _x)
  , mu("mu", "mu", this, _mu)
  , sigma("sigma", "sigma", this, _sigma)
  , k("k", "k", this, _k) {

}

RooGaussExpPdf::RooGaussExpPdf(const RooGaussExpPdf& other, const char* name)
  : RooAbsPdf(other, name)
  , x("x", this, other.x)
  , mu("mu", this, other.mu)
  , sigma("sigma", this, other.sigma)
  , k("k", this, other.k) {

}

Double_t RooGaussExpPdf::evaluate() const {

  double xbar = (x - mu) / sigma;
  double normalization = std::exp(3.0 * std::pow(k, 2) / 2.0) / k + std::sqrt(M_PI / 2.0) * (1.0 - std::erf(k / std::sqrt(2)));
  assert(normalization > 0.0);
  if (xbar >= -k)
    return std::exp(-std::pow(xbar, 2) / 2.0) / normalization;
  return std::exp((std::pow(k, 2) / 2.0) + k * xbar) / normalization;
}
