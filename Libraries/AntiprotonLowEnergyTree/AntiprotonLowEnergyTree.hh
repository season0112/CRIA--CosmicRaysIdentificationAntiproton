#ifndef AntiprotonLowEnergyTree_hh
#define AntiprotonLowEnergyTree_hh

#include "TreeInterface.hh"
#include <bitset>

#define TRACK_FIT_Algorithm 8

namespace Mva{
    class MvaImplementation;
}

class AntiprotonLowEnergyTree : public IO::TreeInterface {
public:
  AntiprotonLowEnergyTree();
  //Event
  IO::TreeBranch<UInt_t> TimeStamp                                             { "TimeStamp",                                               0   };
  IO::TreeBranch<Int_t> Run                                                    { "Run",                                                    -1   };
  IO::TreeBranch<Int_t> EventNumber                                            { "EventNumber",                                            -1   };
  IO::TreeBranch<Double_t> Weight                                              { "Weight",                                                 -1   };
  IO::TreeBranch<bool> IsEventMC                                               { "IsEventMC",                                             false };
  IO::TreeBranch<short> MCParticleID                                           { "MCParticleID",                                           0    };
  IO::TreeBranch<double> MCPrimaryMomentum                                     { "MCPrimaryMomentum",                                      0.0  };
  IO::TreeBranch<short> Pattern                                                { "Pattern",                                                -1   };

  //Trigger
  IO::TreeBranch<UChar_t>  TriggerFlags                                        { "TriggerFlags",                                           0    };

  // RICH
  IO::TreeBranch<Float_t> RICHNumberOfHits                                     { "RICHNumberOfHits",                                     0.0f,  };
  IO::TreeBranch<UChar_t> RichNumberOfRings                                    { "RichNumberOfRings",                                       0,  };
  IO::TreeBranch<Float_t> RichBeta                                             { "RichBeta",                                             0.0f,  };
  IO::TreeBranch<Float_t> RichCharge                                           { "RichCharge",                                           0.0f,  };
  IO::TreeBranch<Float_t> RichIsNaF                                            { "RichIsNaF",                                           -999.0  };
  IO::TreeBranch<Float_t> RICHRINGNumberOfHits                                 { "RICHRINGNumberOfHits",                                -999.0  };
  IO::TreeBranch<Float_t> BetaConsistency                                      { "BetaConsistency",                                     -999.0  };
  IO::TreeBranch<unsigned char> TileIndex                                      { "TileIndex",                                              254  };
  IO::TreeBranch<Float_t> NExpectedPhotoElectrons                              { "NExpectedPhotoElectrons",                             -999.0  };
  IO::TreeBranch<Float_t> NPhotoElectrons                                      { "NPhotoElectrons",                                     -999.0  };
  IO::TreeBranch<Float_t> NCollectedPhotoElectrons                             { "NCollectedPhotoElectrons",                            -999.0  };
  IO::TreeBranch<Float_t> IsGood                                               { "IsGood",                                              -999.0  };
  IO::TreeBranch<Float_t> NumberOfReflectedHits                                { "NumberOfReflectedHits",                               -999.0  };

  //TOF
  IO::TreeBranch<Float_t> TofBeta                                              { "TofBeta",                                           -999.0f,  };
  IO::TreeBranch<Float_t> TofBetaSize                                          { "TofBetaSize",                                       -999.0f,  };
  IO::TreeBranch<unsigned char> TofNumberOfLayers                              { "TofNumberOfLayers",                                       9,  };
  IO::TreeBranch<Float_t> BetaConverted                                        { "BetaConverted",                                     -999.0f,  };
  IO::TreeBranch<Float_t> UpperTofBeta                                         { "UpperTofBeta",                                      -999.0f,  };
  IO::TreeBranch<Float_t> LowerTofBeta                                         { "LowerTofBeta",                                      -999.0f,  };
  IO::TreeBranch<Float_t> TofMassonecharge                                     { "TofMassonecharge",                                   -999.0f, };
  IO::TreeBranch<Float_t> UpperTofCharge                                       { "UpperTofCharge",                                       0.0f,  };
  IO::TreeBranch<Float_t> LowerTofCharge                                       { "LowerTofCharge",                                       0.0f,  };
  IO::TreeBranch<std::vector<Float_t>> BetaFromDeDx                            { "BetaFromDeDx"   ,                      IO::TreeVectorSize(8)  };
  IO::TreeBranch<std::vector<Float_t>> TOFBETACharges                          { "TOFBETACharges",                       IO::TreeVectorSize(4)  };
  IO::TreeBranch<std::vector<Float_t>> TOFClusterEnergy                        { "TOFClusterEnergy",                     IO::TreeVectorSize(4)  };
  IO::TreeBranch<std::vector<Float_t>> TOFClusterCharge                        { "TOFClusterCharge",                     IO::TreeVectorSize(4)  };


  //ECAL
  IO::TreeBranch<double> EcalEnergyElectron                                    { "EcalEnergyElectron",                                     0.0  };
  IO::TreeBranch<char> NumberOfEcalShower                                      { "NumberOfEcalShower",                                     99   };
  IO::TreeBranch<Float_t> EcalBDT_EnergyD                                      { "EcalBDT_EnergyD",                                        -2.0 };
  IO::TreeBranch<Float_t> EcalBDT_EnergyD_Smoothed                             { "EcalBDT_EnergyD_Smoothed",                               -2.0 };

  // TRD
  //IO::TreeBranch<std::vector<Float_t>> dEdX_TRD                                { "dEdX_TRD" ,                            IO::TreeVectorSize(20) };
  IO::TreeBranch<Float_t> TrdLogLikelihoodRatioElectronProtonTracker           { "TrdLogLikelihoodRatioElectronProtonTracker",           0.0f,  };
  IO::TreeBranch<Float_t> TrdLogLikelihoodRatioProtonHeliumTracker             { "TrdLogLikelihoodRatioProtonHeliumTracker",             0.0f,  };
  IO::TreeBranch<Float_t> TRDVTracksSize                                       { "TRDVTracksSize",                                    -999.0f,  };
  IO::TreeBranch<Float_t> TRDHTracksSize                                       { "TRDHTracksSize",                                    -999.0f,  };       
  IO::TreeBranch<Float_t> NumberOfHitsOnTrack                                  { "NumberOfHitsOnTrack",                               -999.0f,  };
  IO::TreeBranch<Float_t> NumberOfLayersWithHit                                { "NumberOfLayersWithHit",                             -999.0f,  };
  IO::TreeBranch<Float_t> TrdSegmentsXZNumber                                  { "TrdSegmentsXZNumber",                               -999.0f,  };
  IO::TreeBranch<Float_t> TrdSegmentsYZNumber                                  { "TrdSegmentsYZNumber",                               -999.0f,  };
  IO::TreeBranch<Float_t> TrdVerticesXZNumber                                  { "TrdVerticesXZNumber",                               -999.0f,  };
  IO::TreeBranch<Float_t> TrdVerticesYZNumber                                  { "TrdVerticesYZNumber",                               -999.0f,  };
  IO::TreeBranch<UShort_t> TrdNumberOfHits                                     { "TrdNumberOfHits",                                           0 };
  IO::TreeBranch<UChar_t> TrdTrackNumberOfSubLayersXZ                          { "TrdTrackNumberOfSubLayersXZ",                               0 };
  IO::TreeBranch<UChar_t> TrdTrackNumberOfSubLayersYZ                          { "TrdTrackNumberOfSubLayersYZ",                               0 };
  IO::TreeBranch<Float_t> TrdKCharge                                           { "TrdKCharge",                                        -999.0f,  };
  //IO::TreeBranch<Float_t> TrdKLikelihood                                       { "TrdKLikelihood",                                    -999.0f,  };
  IO::TreeBranch<UChar_t> TrdActiveLayersTracker                               { "TrdActiveLayersTracker",                                    0 };
  IO::TreeBranch<UChar_t> TrdActiveLayersHybrid                                { "TrdActiveLayersHybrid",                                     0 };
  IO::TreeBranch<UChar_t> TrdActiveLayersStandalone                            { "TrdActiveLayersStandalone",                                 0 };


  //TRACKER
  IO::TreeBranch<Float_t> Rigidity                                             { "Rigidity",                                             0.0f   };
  IO::TreeBranch<std::vector<float>> RigidityValue                             { "RigidityValue",       IO::TreeVectorSize(TRACK_FIT_Algorithm) };
  IO::TreeBranch<unsigned char> ExtrapolatedRichTileIndex                      { "ExtrapolatedRichTileIndex",                             253   };
  IO::TreeBranch<float> ExtrapolatedRichExpectedPhotoElectronsProton           { "ExtrapolatedRichExpectedPhotoElectronsProton",      -999.0f   };
  IO::TreeBranch<float> ExtrapolatedRichExpectedPhotoElectronsPion             { "ExtrapolatedRichExpectedPhotoElectronsPion",        -999.0f   };
  IO::TreeBranch<float> ExtrapolatedRichExpectedPhotoElectronsElectron         { "ExtrapolatedRichExpectedPhotoElectronsElectron",    -999.0f   };
  IO::TreeBranch<float> CCLikelihood                                           { "CCLikelihood",                                      -999.0f   };
  IO::TreeBranch<float> CCBDT                                                  { "CCBDT",                                             -999.0f   };
  IO::TreeBranch<float> CCBDTLapp                                              { "CCBDTLapp",                                         -999.0f   };
  IO::TreeBranch<std::vector<Float_t>> ReconstructedHitLayer                   { "ReconstructedHitLayer",               IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> DepositedEnergyX                        { "DepositedEnergyX",                    IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> DepositedEnergyY                        { "DepositedEnergyY",                    IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> ClusterWidthX                           { "ClusterWidthX",                       IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> ClusterWidthY                           { "ClusterWidthY",                       IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> UnbiasedQX                              { "UnbiasedQX",                          IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> UnbiasedQY                              { "UnbiasedQY",                          IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> ChargeYiJiaXY                           { "ChargeYiJiaXY",                       IO::TreeVectorSize(9)   };
  IO::TreeBranch<std::vector<Float_t>> TrackerCharges                          { "TrackerCharges",                      IO::TreeVectorSize(9)   };

  // ACC 
  IO::TreeBranch<Float_t> ACCHits                                              { "ACCHits",                                           -999.0f,  };
  IO::TreeBranch<Float_t> NumberOfClustersMIT                                  { "NumberOfClustersMIT",                               -999.0f,  };  

  // MVA
  IO::TreeBranch<Float_t> ProtonCCMVABDT                                       { "ProtonCCMVABDT" ,                                        -2.0 };
  IO::TreeBranch<Float_t> ElectronCCMVABDT                                     { "ElectronCCMVABDT",                                       -2.0 };

private:
  virtual void Fill(const Analysis::Event&);
  virtual void UpdateInMemoryBranches();
  virtual IO::TreeInterface* Create() const { return new AntiprotonLowEnergyTree; }
  virtual const IO::TreeBranch<UChar_t>* CurrentTriggerFlags() const { return nullptr; }
  Mva::MvaImplementation* fProtonChargeConfusionMva;
  Mva::MvaImplementation* fElectronChargeConfusionMva;
};
#endif
