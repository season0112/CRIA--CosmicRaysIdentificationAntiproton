#ifndef Statistics_hh
#define Statistics_hh

#include "Quantity.hh"
#include "StackVector.hh"
#include <cmath>
#include <limits>

class TRandom;

namespace MYUtilities {

struct Moments {
  Moments()
    : mean(0)
    , meanError(-1)
    , rms(0)
    , rmsError(-1) {
  }

  double mean;
  double meanError;
  double rms;
  double rmsError;
};

enum MomentsMode {
  AllEntries,
  AbsOfEntries,
  OnlyPositiveEntries
};

template <typename T>
static Moments CalculateMoments(const T& input, MomentsMode mode = AllEntries) {

  Moments moments;
  int inputSize = input.size();
  if (!inputSize)
    return moments;

  double n = 0;
  double sum = 0;
  double sum2 = 0;

  for (int i = 0; i < inputSize; i++) {
    if ((mode == AllEntries)
        || (mode == AbsOfEntries)
        || (mode == OnlyPositiveEntries && input[i] > 0)) {

      n += 1.0;

      if (mode == AbsOfEntries)
        sum += std::abs(input[i]);
      else
        sum += input[i];

      sum2 += input[i] * input[i];
    }
  }

  moments.mean = n ? sum / n : 0;
  if (n > 1) {
    moments.rms = std::sqrt((sum2 / n - moments.mean * moments.mean));
    moments.meanError = moments.rms / std::sqrt(n - 1.0);
    moments.rmsError = moments.rms / std::sqrt(2.0 * n); // ??
  }

  return moments;
}

template <typename T, typename DataType>
static DataType FindMaximum(const T& input, MomentsMode mode = AllEntries) {

  DataType max = -std::numeric_limits<DataType>::max();
  for (unsigned int i = 0; i < input.size(); ++i) {
    if (mode == OnlyPositiveEntries && input.at(i) <= 0)
      continue;
    if (input.at(i) > max)
      max = input.at(i);
  }
  return max;
}

template <typename T, typename DataType>
static DataType FindMinimum(const T& input, MomentsMode mode = AllEntries) {

  DataType min = std::numeric_limits<DataType>::max();
  for (unsigned int i = 0; i < input.size(); ++i) {
    if (mode == OnlyPositiveEntries && input.at(i) <= 0)
      continue;
    if (input.at(i) < min)
      min = input.at(i);
  }
  return min;
}

/** This routine calculates the correct 1-sigma Poisson distributed uncertainty of the integer value val.
  *
  * Above 20 the Square-root is calculated. Below twenty the values given in (Feldman & Cousins, Phys. Rev. D, 1998, Table 2, no background) are used.
  * So all uncertaintites, also for small numbers, are correct.
  */
void PoissonUncertainty(int val, double& lowerError, double& upperError);

double PoissonPvalue(unsigned int nEvents, double expectedBg);
double PoissonSignificance(unsigned int nEvents, double expectedBg);


/** Li&Ma significance.
  *
  * See Li, T.-P. & Ma, Y.-Q., Analysis methods for results in gamma-ray astronomy ApJ, 1983, 272, 317-324
  */
double LiMaSignificance(int nOn, int nOff, double alpha);

/** Generate a random number distributed according to a power law within a given interval.
  *
  * \param[in] xlow Lower bound of interval.
  * \param[in] xup Upper bound of interval.
  * \param[in] spectralIndex Spectral index \f$ \gamma \f$ of power law  (\f$ E^gamma \f$).
  * \param[in] rand Random generator to use (or will create a static one for this function if none is provided).
  */
double PowerLawRandomNumber(double xlow, double xup, double spectralIndex, TRandom* rand = nullptr);

template <typename T>
Quantity WeightedMean(const T& values, const T& errors) {

  assert(values.size() == errors.size());

  double wSum = 0.0;
  double wValSum = 0.0;

  for (unsigned int i = 0; i < values.size(); ++i) {

    double w = 1.0 / std::pow(errors[i], 2);
    wSum += w;
    wValSum += w * values[i];
  }

  double wmean = wSum > 0.0 ? wValSum / wSum : 0.0;
  double wmeanErr = wSum > 0.0 ? std::sqrt(1.0 / wSum) : 0.0;

  return Quantity(wmean, wmeanErr);
}

} // namespace MYUtilities

#endif
