#ifndef ExampleAnalysisTree_hh
#define ExampleAnalysisTree_hh

#include "TreeInterface.hh"

namespace Mva{
    class MvaImplementation;
}

class ExampleAnalysisTree : public IO::TreeInterface {
public:
  ExampleAnalysisTree();
  //IO::TreeBranch<Float_t> TrackerChiSquareX           { "TrackerChiSquareX"    ,  0   };
  //IO::TreeBranch<Float_t> TrackerChiSquareY           { "TrackerChiSquareY"    ,  0   };
  IO::TreeBranch<Float_t> InnerTrackerCharge          { "InnerTrackerCharge"   ,  IO::ValueLimitMode::HighestValue };
  IO::TreeBranch<std::vector<Float_t>> TrackerCharges { "TrackerCharges"       ,  IO::TreeVectorSize(9) };
  IO::TreeBranch<Float_t> Rigidity                    { "Rigidity"             ,  0   };
  IO::TreeBranch<Float_t> RichBeta                    { "RichBeta",             0.0f,  };

  //IO::TreeBranch<std::vector<Float_t>> HitXCoordinate { "HitXCoordinate"       ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> HitYCoordinate { "HitYCoordinate"       ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> HitZLayer      { "HitZLayer"       ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> ChargeX        { "ChargeX"              ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> ChargeY        { "ChargeY"              ,  IO::TreeVectorSize(9) };

  IO::TreeBranch<UInt_t> TimeStamp                    { "TimeStamp"                 ,  0   };
  IO::TreeBranch<Int_t> Run                           { "Run"                  , -1   };
  IO::TreeBranch<Int_t> EventNumber                   { "EventNumber"          , -1   };
  IO::TreeBranch<Double_t> Weight                     { "Weight"               , -1   };
  IO::TreeBranch<bool> IsEventMC                      { "IsEventMC", false };
  IO::TreeBranch<short> MCParticleID                  { "MCParticleID", 0 };
  IO::TreeBranch<double> MCPrimaryMomentum            { "MCPrimaryMomentum", 0.0 };
  IO::TreeBranch<UChar_t>  TriggerFlags               { "TriggerFlags"         , 0    };

  //IO::TreeBranch<std::vector<Float_t>> ResidualYInMicrons   { "ResidualYInMicrons"     ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> ClusterWidthX        { "ClusterWidthX"              ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> ClusterWidthY        { "ClusterWidthY"              ,  IO::TreeVectorSize(9) };

  IO::TreeBranch<short> Pattern                             { "Pattern"       , -1   };

  //IO::TreeBranch<Float_t> ZTrackerLayer1                    { "ZTrackerLayer1"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer2                    { "ZTrackerLayer2"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer3                    { "ZTrackerLayer3"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer4                    { "ZTrackerLayer4"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer5                    { "ZTrackerLayer5"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer6                    { "ZTrackerLayer6"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer7                    { "ZTrackerLayer7"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer8                    { "ZTrackerLayer8"             ,  0   };
  //IO::TreeBranch<Float_t> ZTrackerLayer9                    { "ZTrackerLayer9"             ,  0   };
 
  //IO::TreeBranch<std::vector<Float_t>> UnassociatedHitZLayer      { "UnassociatedHitZLayer"       ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> UnassociatedHitXCoordinate { "UnassociatedHitXCoordinate"       ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> UnassociatedHitYCoordinate { "UnassociatedHitYCoordinate"       ,  IO::TreeVectorSize(9) }; 
  //IO::TreeBranch<std::vector<Float_t>> UnassociatedChargeX        { "UnassociatedChargeX"              ,  IO::TreeVectorSize(9) };  
  //IO::TreeBranch<std::vector<Float_t>> UnassociatedChargeY        { "UnassociatedChargeY"              ,  IO::TreeVectorSize(9) };  
  //IO::TreeBranch<std::vector<Float_t>> RigidityWithoutThisHit     { "RigidityWithoutThisHit"           ,  IO::TreeVectorSize(9) }; 

  // Charge Confusion Ecal BDT
  IO::TreeBranch<double> EcalEnergyElectron                       { "EcalEnergyElectron"   ,  0.0 };
  IO::TreeBranch<Float_t> ElectronCCMVABDT                        { "ElectronCCMVABDT",               -2.0 };
  // Ecal BDT for Proton and Electron Seperation
  IO::TreeBranch<Float_t> EcalBDT_EnergyD                         { "EcalBDT_EnergyD"                , -2.0 };
  IO::TreeBranch<Float_t> EcalBDT_EnergyD_Smoothed                { "EcalBDT_EnergyD_Smoothed"       , -2.0 };


  //Charge Confusion
  //IO::TreeBranch<Float_t> CCBDT                                   { "CCBDT"                ,  0.0 };
  //IO::TreeBranch<Float_t> CCBDTLapp                               { "CCBDTLapp"            ,  0.0 };
  //IO::TreeBranch<Float_t> CCLikelihood                            { "CCLikelihood"         ,  0.0 };
  IO::TreeBranch<Float_t> ProtonCCMVABDT                          { "ProtonCCMVABDT"       , -2.0 };

  // MVA
  IO::TreeBranch<Float_t> RigidityInnerL1                { "RigidityInnerL1",             0.0f,  };
  IO::TreeBranch<Float_t> RigidityInnerL9                { "RigidityInnerL9",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerXInnerL1            { "Chi2TrackerXInnerL1",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerYInnerL1            { "Chi2TrackerYInnerL1",             0.0f,  };

  IO::TreeBranch<Float_t> RigidityInner                  { "RigidityInner",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerXInnerL9            { "Chi2TrackerXInnerL9",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerYInnerL9            { "Chi2TrackerYInnerL9",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerXInner              { "Chi2TrackerXInner",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerYInner              { "Chi2TrackerYInner",             0.0f,  };

  IO::TreeBranch<Float_t> Chi2TrackerX                   { "Chi2TrackerX",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerY                   { "Chi2TrackerY",             0.0f,  };
  IO::TreeBranch<Float_t> InverseRigidityError           { "InverseRigidityError",             0.0f,  };
  IO::TreeBranch<Float_t> YResidualJLayer1               { "YResidualJLayer1",             0.0f,  };
  IO::TreeBranch<Float_t> YResidualJLayer9               { "YResidualJLayer9",             0.0f,  };

  IO::TreeBranch<Float_t> RigidityInnerUpperHalf               { "RigidityInnerUpperHalf",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerXInnerUpperHalf           { "Chi2TrackerXInnerUpperHalf",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerYInnerUpperHalf           { "Chi2TrackerYInnerUpperHalf",             0.0f,  };
  IO::TreeBranch<Float_t> RigidityInnerLowerHalf               { "RigidityInnerLowerHalf",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerXInnerLowerHalf           { "Chi2TrackerXInnerLowerHalf",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerYInnerLowerHalf           { "Chi2TrackerYInnerLowerHalf",             0.0f,  };

  IO::TreeBranch<Float_t> RigidityAsymmetry                    { "RigidityAsymmetry",             0.0f,  };
  IO::TreeBranch<Float_t> RigidityAsymmetryL9                  { "RigidityAsymmetryL9",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerYAsymmetry                { "Chi2TrackerYAsymmetry",             0.0f,  };
  IO::TreeBranch<Float_t> InnerMaxSpanRigidityMatching         { "InnerMaxSpanRigidityMatching",             0.0f,  };
  IO::TreeBranch<Float_t> InnerMaxSpanRigidityMatchingv2       { "InnerMaxSpanRigidityMatchingv2",             0.0f,  };
  IO::TreeBranch<Float_t> L1L9RigidityMatching                 { "L1L9RigidityMatching",             0.0f,  };
  IO::TreeBranch<Float_t> L1L9RigidityMatchingv2               { "L1L9RigidityMatchingv2",             0.0f,  };

  IO::TreeBranch<Float_t> L24L58RigidityMatching               { "L24L58RigidityMatching",             0.0f,  };
  IO::TreeBranch<Float_t> L24L58RigidityMatchingv2             { "L24L58RigidityMatchingv2",             0.0f,  };

  IO::TreeBranch<Float_t> Log10Chi2TrackerYInner               { "Log10Chi2TrackerYInner",             0.0f,  };
  IO::TreeBranch<Float_t> Log10Chi2TrackerX                    { "Log10Chi2TrackerX",             0.0f,  };
  IO::TreeBranch<Float_t> Log10Chi2TrackerY                    { "Log10Chi2TrackerY",             0.0f,  };
  IO::TreeBranch<Float_t> TrackerL58L24ChargeAsymmetry         { "TrackerL58L24ChargeAsymmetry",             0.0f,  };

  IO::TreeBranch<Float_t> TrackerL9Charge                      { "TrackerL9Charge",             0.0f,  };
  IO::TreeBranch<Float_t> TrackerL78Charge                     { "TrackerL78Charge",             0.0f,  };
  IO::TreeBranch<Float_t> UpperTofCharge                       { "UpperTofCharge",             0.0f,  };
  IO::TreeBranch<Float_t> LowerTofCharge                       { "LowerTofCharge",             0.0f,  };

  IO::TreeBranch<Float_t> TrdLogLikelihoodRatioElectronProtonTracker           { "TrdLogLikelihoodRatioElectronProtonTracker",             0.0f,  };
  IO::TreeBranch<Float_t> TrdLogLikelihoodRatioProtonHeliumTracker             { "TrdLogLikelihoodRatioProtonHeliumTracker",             0.0f,  };
  
  IO::TreeBranch<Float_t> InverseRigidityErrorInner           { "InverseRigidityErrorInner",             0.0f,  };
  IO::TreeBranch<Float_t> RigidityAlcaraz                     { "RigidityAlcaraz",             0.0f,  };
  IO::TreeBranch<Float_t> InverseRigidityErrorAlcaraz         { "InverseRigidityErrorAlcaraz",             0.0f,  };
  IO::TreeBranch<Float_t> RigidityChikanian                   { "RigidityChikanian",             0.0f,  };
  IO::TreeBranch<Float_t> Chi2TrackerYChikanian               { "Chi2TrackerYChikanian",             0.0f,  };
  IO::TreeBranch<Float_t> InverseRigidityErrorChikanian       { "InverseRigidityErrorChikanian",             0.0f,  };
  IO::TreeBranch<Float_t> Log10Chi2TrackerXInner              { "Log10Chi2TrackerXInner",             0.0f,  };
  
  /*
  IO::TreeBranch<Float_t> XLayer1               { "XLayer1",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer1               { "YLayer1",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer2               { "XLayer2",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer2               { "YLayer2",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer3               { "XLayer3",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer3               { "YLayer3",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer4               { "XLayer4",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer4               { "YLayer4",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer5               { "XLayer5",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer5               { "YLayer5",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer6               { "XLayer6",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer6               { "YLayer6",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer7               { "XLayer7",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer7               { "YLayer7",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer8               { "XLayer8",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer8               { "YLayer8",             0.0f,  };
  IO::TreeBranch<Float_t> XLayer9               { "XLayer9",             0.0f,  };
  IO::TreeBranch<Float_t> YLayer9               { "YLayer9",             0.0f,  };
  */

  //IO::TreeBranch<std::vector<Float_t>> XLayer        { "XLayer"              ,  IO::TreeVectorSize(9) };
  //IO::TreeBranch<std::vector<Float_t>> YLayer        { "YLayer"              ,  IO::TreeVectorSize(9) };

  // ACC
  IO::TreeBranch<std::vector<Float_t>> ACCTime                   { "ACCTime",                  IO::TreeVectorSize(9)  };
  IO::TreeBranch<std::vector<Float_t>> ACCZ                      { "ACCZ",                     IO::TreeVectorSize(9)  };
  IO::TreeBranch<std::vector<Float_t>> ACCRawCharge              { "ACCRawCharge",             IO::TreeVectorSize(9)  };

private:
  virtual void Fill(const Analysis::Event&);
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new ExampleAnalysisTree; }
  virtual const IO::TreeBranch<UChar_t>* CurrentTriggerFlags() const { return nullptr; }
  Mva::MvaImplementation* fProtonChargeConfusionMva;
  Mva::MvaImplementation* fElectronChargeConfusionMva;
};
#endif
