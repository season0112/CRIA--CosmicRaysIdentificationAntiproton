#ifndef MigrationRandomGenerator_hh
#define MigrationRandomGenerator_hh

#include <BinningDefinition.hh>

#include <map>
#include <vector>

class MigrationParameterization;

class TRandom;

class MigrationRandomGenerator {

public:
  MigrationRandomGenerator(const Binning::Definition& binning, const MigrationParameterization* param);

  double GetRandomSmearedEnergy(double Etruth) const;

  void SetRandomGenerator(TRandom* rand);

private:

  mutable TRandom* fRandom = nullptr;

  Binning::Definition fEnergyBinning;

  const MigrationParameterization* fParam = nullptr;

  std::map<unsigned int,std::vector<double>> fIntegralMap; // map bin number to vector of cdf values
};

#endif

