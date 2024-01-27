#ifndef AntiprotonIntermediateEnergyTree_hh
#define AntiprotonIntermediateEnergyTree_hh

#include "TreeInterface.hh"
#include <bitset>

#define TRACK_FIT_Algorithm 8

namespace Mva{
    class MvaImplementation;
}

class AntiprotonIntermediateEnergyTree : public IO::TreeInterface {
public:
  AntiprotonIntermediateEnergyTree();
  IO::TreeBranch<Float_t> RichBeta                                             { "RichBeta",            0.0f,  };
  IO::TreeBranch<UChar_t> RichNumberOfRings                                    { "RichNumberOfRings",      0,  };
  IO::TreeBranch<Bool_t> RichIsNaF                                             { "RichIsNaF",           false  };


  IO::TreeBranch<Float_t> TofBeta                                              { "TofBeta",             0.0f,  };
  IO::TreeBranch<UInt_t> TimeStamp                                             { "TimeStamp"            ,  0   };
  IO::TreeBranch<Int_t> Run                                                    { "Run"                  , -1   };
  IO::TreeBranch<Int_t> EventNumber                                            { "EventNumber"          , -1   };
  IO::TreeBranch<Double_t> Weight                                              { "Weight"               , -1   };
  IO::TreeBranch<bool> IsEventMC                                               { "IsEventMC"            ,false };
  IO::TreeBranch<short> MCParticleID                                           { "MCParticleID"         , 0    };
  IO::TreeBranch<double> MCPrimaryMomentum                                     { "MCPrimaryMomentum"    , 0.0  };
  IO::TreeBranch<short> Pattern                                                { "Pattern"              , -1   };
  //  IO::TreeBranch<unsigned char> TriggerFlags                               { "TriggerFlags"         , 9    };
  IO::TreeBranch<UChar_t>  TriggerFlags                                        { "TriggerFlags"         , 0    };
  //  IO::TreeBranch<std::bitset<8>>  TriggerFlagsBitset                       { "TriggerFlagsBitset"   , ：00000000};
  IO::TreeBranch<double> EcalEnergyElectron                                    { "EcalEnergyElectron"   ,  0.0 };
  IO::TreeBranch<char> NumberOfEcalShower                                      { "NumberOfEcalShower"   , 99   };
  IO::TreeBranch<unsigned char> TofNumberOfLayers                              { "TofNumberOfLayers"    , 9    };


// dE/dX
  IO::TreeBranch<double> dEdX                                                  { "dEdX"                 ,-999.0}; 

// MVA
  IO::TreeBranch<Float_t> TrdLogLikelihoodRatioElectronProtonTracker           { "TrdLogLikelihoodRatioElectronProtonTracker",           0.0f,  };
  IO::TreeBranch<Float_t> TrdLogLikelihoodRatioProtonHeliumTracker             { "TrdLogLikelihoodRatioProtonHeliumTracker",             0.0f,  };
  IO::TreeBranch<Float_t> ProtonCCMVABDT                                       { "ProtonCCMVABDT" ,                                        -2.0 };
  IO::TreeBranch<Float_t> ElectronCCMVABDT                                     { "ElectronCCMVABDT",               -2.0 };

// Ecal BDT for Proton and Electron Seperation
  IO::TreeBranch<Float_t> EcalBDT_EnergyD                                      { "EcalBDT_EnergyD"                , -2.0 };
  IO::TreeBranch<Float_t> EcalBDT_EnergyD_Smoothed                             { "EcalBDT_EnergyD_Smoothed"       , -2.0 };


// Different Rigidity
  IO::TreeBranch<std::vector<float>> RigidityValue                             { "RigidityValue", IO::TreeVectorSize(TRACK_FIT_Algorithm) };


// MVA (simplied version)
  IO::TreeBranch<Float_t> Rigidity                               { "Rigidity",                               0.0f   };
  IO::TreeBranch<std::vector<Float_t>> TrackerCharges            { "TrackerCharges",          IO::TreeVectorSize(9) };
//  IO::TreeBranch<Float_t> InnerTrackerCharge                   { "InnerTrackerCharge"   ,  IO::ValueLimitMode::HighestValue };
//  IO::TreeBranch<Float_t> RigidityInnerL1                      { "RigidityInnerL1",                        0.0f,  };
//  IO::TreeBranch<Float_t> RigidityInnerL9                      { "RigidityInnerL9",                        0.0f,  };
//  IO::TreeBranch<Float_t> RigidityInner                        { "RigidityInner",                          0.0f,  };
//  IO::TreeBranch<Float_t> Chi2TrackerXInner                    { "Chi2TrackerXInner",                      0.0f,  };
//  IO::TreeBranch<Float_t> Chi2TrackerYInner                    { "Chi2TrackerYInner",                      0.0f,  };
//  IO::TreeBranch<Float_t> Chi2TrackerX                         { "Chi2TrackerX",                           0.0f,  };
//  IO::TreeBranch<Float_t> Chi2TrackerY                         { "Chi2TrackerY",                           0.0f,  };
//  IO::TreeBranch<Float_t> RigidityInnerUpperHalf               { "RigidityInnerUpperHalf",                 0.0f,  };
//  IO::TreeBranch<Float_t> RigidityInnerLowerHalf               { "RigidityInnerLowerHalf",                 0.0f,  };
//  IO::TreeBranch<Float_t> RigidityAsymmetry                    { "RigidityAsymmetry",                      0.0f,  };
//  IO::TreeBranch<Float_t> RigidityAsymmetryL9                  { "RigidityAsymmetryL9",                    0.0f,  };
//  IO::TreeBranch<Float_t> Chi2TrackerYAsymmetry                { "Chi2TrackerYAsymmetry",                  0.0f,  };
//  IO::TreeBranch<Float_t> InnerMaxSpanRigidityMatching         { "InnerMaxSpanRigidityMatching",           0.0f,  };
//  IO::TreeBranch<Float_t> L1L9RigidityMatching                 { "L1L9RigidityMatching",                   0.0f,  };
//  IO::TreeBranch<Float_t> L24L58RigidityMatching               { "L24L58RigidityMatching",                 0.0f,  };
//  IO::TreeBranch<Float_t> Log10Chi2TrackerYInner               { "Log10Chi2TrackerYInner",                 0.0f,  };
//  IO::TreeBranch<Float_t> Log10Chi2TrackerX                    { "Log10Chi2TrackerX",                      0.0f,  };
//  IO::TreeBranch<Float_t> Log10Chi2TrackerY                    { "Log10Chi2TrackerY",                      0.0f,  };
//  IO::TreeBranch<Float_t> TrackerL58L24ChargeAsymmetry         { "TrackerL58L24ChargeAsymmetry",           0.0f,  };
  IO::TreeBranch<Float_t> TrackerL9Charge                        { "TrackerL9Charge",                        0.0f,  };
  IO::TreeBranch<Float_t> TrackerL78Charge                       { "TrackerL78Charge",                       0.0f,  };
//  IO::TreeBranch<Float_t> UpperTofCharge                       { "UpperTofCharge",                         0.0f,  };
  IO::TreeBranch<Float_t> LowerTofCharge                         { "LowerTofCharge",                         0.0f,  };
//  IO::TreeBranch<Float_t> Log10Chi2TrackerXInner               { "Log10Chi2TrackerXInner",                 0.0f,  };






private:
  virtual void Fill(const Analysis::Event&);
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new AntiprotonIntermediateEnergyTree; }
  virtual const IO::TreeBranch<UChar_t>* CurrentTriggerFlags() const { return nullptr; }
  Mva::MvaImplementation* fProtonChargeConfusionMva;
  Mva::MvaImplementation* fElectronChargeConfusionMva;
};
#endif
