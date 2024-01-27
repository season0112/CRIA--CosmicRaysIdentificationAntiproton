#include "RooAsymmLandauLikePdf.hh"

#include <cmath>

#include <TMath.h>

RooAsymmLandauLikePdf::RooAsymmLandauLikePdf(const char *name, const char *title,
                                             RooAbsReal& _x, RooAbsReal& _peak, RooAbsReal& _width, RooAbsReal& _asymm)
  : RooAbsPdf(name, title)
  , x("x", "x", this, _x)
  , peak("peak", "peak", this, _peak)
  , width("width", "width", this, _width)
  , asymm("asymm", "asymm", this, _asymm) {

}

RooAsymmLandauLikePdf::RooAsymmLandauLikePdf(const RooAsymmLandauLikePdf& other, const char* name)
  : RooAbsPdf(other, name)
  , x("x", this, other.x)
  , peak("peak", this, other.peak)
  , width("width", this, other.width)
  , asymm("asymm", this, other.asymm) {

}

Double_t RooAsymmLandauLikePdf::evaluate() const {

  double var = (x - peak) / width / asymm;
  double val = std::pow(asymm, asymm - 1) / TMath::Gamma(asymm) * std::exp(-asymm * (var + std::exp(-var))) / width;
  return val > 0.0 ? val : 1e-30;
}
