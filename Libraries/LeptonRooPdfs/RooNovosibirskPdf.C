#include "RooNovosibirskPdf.hh"

#include <cmath>

RooNovosibirskPdf::RooNovosibirskPdf(const char *name, const char *title, RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _tau)
  : RooAbsPdf(name, title)
  , x("x", "x", this, _x)
  , mu("mu", "mu", this, _mu)
  , sigma("sigma", "sigma", this, _sigma)
  , tau("tau", "tau", this, _tau) {

}

RooNovosibirskPdf::RooNovosibirskPdf(const RooNovosibirskPdf& other, const char* name)
  : RooAbsPdf(other, name)
  , x("x", this, other.x)
  , mu("mu", this, other.mu)
  , sigma("sigma", this, other.sigma)
  , tau("tau", this, other.tau) {

}

Double_t RooNovosibirskPdf::evaluate() const {

  static const double sSqrtLog4 = std::sqrt(std::log(4.0));

  double y = 0.0;
  if (tau == 0.0)
    y = std::exp(-0.5 * (std::pow(x - mu, 2) / std::pow(sigma, 2)));
  else {
    double xLambda = 1.0 + std::sinh(tau * sSqrtLog4) / (sigma * tau * sSqrtLog4) * tau * (x - mu);
    y = std::exp(-0.5 * (std::pow(std::log(xLambda), 2) / std::pow(tau, 2) + std::pow(tau, 2)));
  }

  if (y <= 0.0 || std::isinf(y) || std::isnan(y))
    return 1e-30;
  return y;
}
