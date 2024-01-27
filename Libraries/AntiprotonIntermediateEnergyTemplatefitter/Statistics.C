#include "Statistics.hh"

#include <Math/ProbFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>
#include <TMath.h>
#include <TRandom3.h>

#include <cassert>

#define INFO_OUT_TAG "Statistics"
#include "debugging.hh"

namespace MYUtilities {

void PoissonUncertainty(int val, double& lowerError, double& upperError) {

  if (val > 20) { //sqrt
    lowerError = std::sqrt(val);
    upperError = std::sqrt(val);
  }
  else {
    switch (val) {
    case 0:
      lowerError = 0.0;
      upperError = 1.29;
      break;
    case 1:
      lowerError = 0.63;
      upperError = 1.75;
      break;
    case 2:
      lowerError = 1.26;
      upperError = 2.25;
      break;
    case 3:
      lowerError = 1.90;
      upperError = 2.30;
      break;
    case 4:
      lowerError = 1.66;
      upperError = 2.78;
      break;
    case 5:
      lowerError = 2.25;
      upperError = 2.81;
      break;
    case 6:
      lowerError = 2.18;
      upperError = 3.28;
      break;
    case 7:
      lowerError = 2.75;
      upperError = 3.30;
      break;
    case 8:
      lowerError = 2.70;
      upperError = 3.32;
      break;
    case 9:
      lowerError = 2.67;
      upperError = 3.79;
      break;
    case 10:
      lowerError = 3.22;
      upperError = 3.81;
      break;
    case 11:
      lowerError = 3.19;
      upperError = 3.82;
      break;
    case 12:
      lowerError = 3.17;
      upperError = 4.29;
      break;
    case 13:
      lowerError = 3.73;
      upperError = 4.30;
      break;
    case 14:
      lowerError = 3.70;
      upperError = 4.32;
      break;
    case 15:
      lowerError = 3.68;
      upperError = 4.32;
      break;
    case 16:
      lowerError = 3.67;
      upperError = 4.80;
      break;
    case 17:
      lowerError = 4.21;
      upperError = 4.81;
      break;
    case 18:
      lowerError = 4.19;
      upperError = 4.82;
      break;
    case 19:
      lowerError = 4.18;
      upperError = 4.82;
      break;
    case 20:
      lowerError = 4.17;
      upperError = 5.30;
      break;
    default:
      lowerError = 0;
      upperError = 0;
      std::cout << "WARNING! No valid input for PoissonError" << std::endl;
    } //end switch
  }   // end else

  return;
}


double PoissonPvalue(unsigned int nEvents, double expectedBg) {

  if (!nEvents)
    return 1.0;

  return 1.0 - ROOT::Math::poisson_cdf(nEvents - 1, expectedBg);
}


double PoissonSignificance(unsigned int nEvents, double expectedBg) {

  if (double(nEvents) <= expectedBg)
    return 0.0;

  double k = double(nEvents);
  double r = -2.0 * (k * std::log(expectedBg / k) + k - expectedBg);
  return (std::sqrt(r));
}

double LiMaSignificance(int nOn, int nOff, double alpha) {

  assert(alpha >= 0.0);
  if (alpha < 0.0001)
    return std::sqrt(1.0 * nOn);

  double excess = nOn - nOff * alpha;
  double n1, n2, sign, a;

  if (excess > 0.0) {
    n1 = nOn;
    n2 = nOff;
    sign = 1.0;
    a = alpha;
  }
  else {
    n2 = nOn;
    n1 = nOff;
    sign = -1.0;
    a = 1.0 / alpha;
  }

  // nOn = 0 or nOff = 0
  if (n2 == 0)
    return sign * std::sqrt(2 * n1 * std::log((1 + a) / a));
  if (n1 == 0)
    return sign * std::sqrt(2 * n2 * std::log(1 + a));

  // the standard Li & Ma formula:
  double nt = n1 + n2;
  double pa = 1 + a;

  double t1 = n1 * std::log((pa / a) * (n1 / nt));
  double t2 = n2 * std::log(pa * (n2 / nt));
  double sig = std::sqrt(2) * std::sqrt(std::abs(t1 + t2));

  return sign * sig;
}


double PowerLawRandomNumber(double xlow, double xup, double spectralIndex, TRandom* rand) {

  double r = 0.0;
  if (rand)
    r = rand->Uniform();
  else {
    static TRandom3* sRandomGenerator = new TRandom3(0);
    r = sRandomGenerator->Uniform();
  }

  if (std::abs(spectralIndex + 1.0) > 0.001) {

    // standard case
    double a = std::pow(xlow, spectralIndex + 1.0);
    double b = std::pow(xup, spectralIndex + 1.0);
    return std::pow((b - a) * r + a, 1.0 / (spectralIndex + 1.0));
  }
  else {

    // E^(-1) special case
    return xlow * std::pow(xup / xlow, r);
  }
}


} // namespace MYUtilities
