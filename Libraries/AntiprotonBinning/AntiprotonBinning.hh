#ifndef AntiprotonBinning_hh
#define AntiprotonBinning_hh
#include "BinningFunctions.hh"
#include "BinningTools.hh"
#include <vector>
#include <math.h>


namespace AntiprotonNewBinning {

int SplitTotal_3B = 135; 
int SplitTotal_6B = 138;
int SplitTotal_6M = 21;
int mergestep_3B = 3;
int mergestep_6B = 6;
int mergestep_6M = 1;

class Definition;

class NewBinning {
public:
  static const Binning::Definition& AntiprotonBinning450();
  static const Binning::Definition& AntiprotonBinning525_zhili();
  static const Binning::Definition& AntiprotonBinCenter450();
  static const Binning::Definition& AntiprotonBinCenter525_zhili();
};


class AntiprotonAllBinning {
public:
  static const std::vector<double> AntiprotonBinning_450;
  static const std::vector<double> AntiprotonBinning_525;
  static const std::vector<double> AntiprotonBinning_zhili525;
  static const std::vector<double> AntiprotonBinningCenter_450;
  static const std::vector<double> AntiprotonBinningCenter_525;
  static const std::vector<double> AntiprotonBinningCenter_zhili525;
};


class AntiprotonResults {
public:
  static const std::vector<double> PublishedRatioPRL;
  static const std::vector<double> PublishedRatioErrorPRL;
  static const std::vector<double> PublishedRatioStatisticErrorPRL;
  static const std::vector<double> PublishedRatioSystematicErrorPRL;
  static const std::vector<double> PublishedRatioRelativeErrorPRL;
  static const std::vector<double> PublishedRatioStatisticRelativeErrorPRL;
  static const std::vector<double> PublishedRatioSystematicRelativeErrorPRL;
  static const std::vector<double> PublishedRatioPRL_pbarNumber;

  static const std::vector<double> PhysicsReportRatio;
  static const std::vector<double> PhysicsReportRatioError;
  static const std::vector<double> PhysicsReportRatioStatisticError;
  static const std::vector<double> PhysicsReportRatioSystematicError;
  static const std::vector<double> PhysicsReportRatioRelativeError;
  static const std::vector<double> PhysicsReportStatisticRatioRelativeError;
  static const std::vector<double> PhysicsReportSystematicRatioRelativeError;
  static const std::vector<double> PhysicsReport_pbarNumber;

  static const std::vector<double> PublishedFluxPRL;
  static const std::vector<double> PublishedFluxErrorPRL;
  static const std::vector<double> PublishedFluxStatisticErrorPRL;
  static const std::vector<double> PublishedFluxSystematicErrorPRL;
  static const std::vector<double> PublishedFluxStatisticRelativeErrorPRL;
  static const std::vector<double> PublishedFluxSystematicRelativeErrorPRL;

  static const std::vector<double> PublishedProtonFluxPRL;
  static const std::vector<double> PublishedProtonFluxErrorPRL;
  static const std::vector<double> PublishedProtonFluxStatisticErrorPRL;
  static const std::vector<double> PublishedProtonFluxSystematicErrorPRL;
};

} //namespace AntiprotonNewBinning


#endif
