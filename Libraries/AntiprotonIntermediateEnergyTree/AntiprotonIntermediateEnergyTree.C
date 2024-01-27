#include "AnalysisEvent.hh"
#include "AntiprotonIntermediateEnergyTree.hh"

#include "AMSGeometry.h"
#include "TrackerTrack.h"
#include "MvaImplementation.hh"
#include "MvaInterface.hh"
#include "Tracker.h"
#include "TrackerTrackFit.h"
#include "TOFBeta.h"
#include "TOF.h"
#include "Trigger.h"
#include "Event.h"
#include <TrdLikelihoodCalculation.hh>
#include <cmath>
#include "Utilities.hh"
#include "AnalysisParticle.hh"
#include "EcalLongitudinalShowerFit.hh"
#include "Environment.hh"
#include "ParticleId.hh"
#include "ResolutionModels.hh"
#include "SlowControlLookup.hh"
#include <TrackerTrackCoordinates.h>
#include "ECALShower.h"
#include "EcalLongitudinalShowerFit.hh"
#include "ECAL.h"
#include "MC.h"
#define INFO_OUT_TAG "AntiprotonIntermediateEnergyTree"
#include "debugging.hh"

AntiprotonIntermediateEnergyTree::AntiprotonIntermediateEnergyTree()
  : IO::TreeInterface("AntiprotonIntermediateEnergyTree", "Antiproton Intermediate Energy Tree") {
  RegisterBranches();
 // Load proton charge-confusion MVA
   const auto* mvaInterfaceProtonCC = Mva::MvaInterface::CreateMvaInterfaceByName("Mva::ProtonChargeConfusionMvaInterface");
      fProtonChargeConfusionMva = mvaInterfaceProtonCC->CreateMvaImplementation();
 // Load electron charge-confusion MVA
    const auto* mvaInterfaceElectronCC = Mva::MvaInterface::CreateMvaInterfaceByName("Mva::ElectronChargeConfusionMvaInterface");
      fElectronChargeConfusionMva = mvaInterfaceElectronCC->CreateMvaImplementation();
}

template <typename T>
double CalculateMeanChargeForLayers(unsigned int first, unsigned int last, const std::vector<T>& charges) {
  double charge = 0;
  double numberOfLayers = 0;
  assert(first > 0);
  assert(last <= charges.size());
  for (unsigned int i = first; i <= last; ++i) {
    const double layerCharge = charges[i - 1];
    if (layerCharge > 0) {
      charge += layerCharge;
      ++numberOfLayers;
    }
  }
  return charge / numberOfLayers;
}

void AntiprotonIntermediateEnergyTree::Fill(const Analysis::Event& event) {
  TimeStamp = event.TimeStamp().GetSec();
  Run = event.Run();
  EventNumber = event.EventNumber();
  Weight = event.Weight();
  IsEventMC = event.IsMC();
  TriggerFlags = event.TriggerFlags();
  NumberOfEcalShower = event.NumberOfEcalShower();
  TofNumberOfLayers = event.TofNumberOfLayers();

  RichNumberOfRings = event.NumberOfRichRings();

  if (event.IsMC()) {
    MCPrimaryMomentum = event.McMomentum();
    MCParticleID = event.McParticleId();
  }
  if (const Analysis::Particle* primaryParticle = event.PrimaryParticle()) {  // { event.PrimaryParticle

      if (const AC::TrackerTrack* trackertrack=primaryParticle->TrackerTrack()){  // {primaryParticle->TrackerTrack
        Pattern=trackertrack->TrackerLayerPatternClassification();
        const int refitPattern = AC::PGMA + AC::RebuildFromTDV;
        const AC::ParticleHypothesis particleHypothesis = AC::DefaultMass;

    const size_t fitAlgorithms = 8;
    RigidityValue().resize(fitAlgorithms);
         for (auto fitAlgorithm : {
            AC::Default, AC::Choutko, AC::Alcaraz, AC::ChikanianF, AC::ChikanianC, AC::Kalman, AC::VertexFit}) {
          int fitIndex = trackertrack->GetFitFuzzy(
           fitAlgorithm , AC::All, refitPattern, particleHypothesis);
          if (fitIndex >= 0) {
            const auto& trackFit = trackertrack->TrackFits().at(fitIndex);
            RigidityValue().at(fitAlgorithm) = trackFit.Rigidity();
          }
        }

  ////    MVA  begin
  int iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::InnerPlusL1, refitPattern, particleHypothesis);
//  if (iFit >= 0) {
//    const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
//    RigidityInnerL1 = trackFit.Rigidity();
//  }

//  iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::InnerPlusL9, refitPattern, particleHypothesis);
//  if (iFit >= 0) {
//    const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
//    RigidityInnerL9 = trackFit.Rigidity();
//  }

  iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::Inner, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
    Rigidity = trackFit.Rigidity();
//    RigidityInner = trackFit.Rigidity();
//    Chi2TrackerXInner = trackFit.ChiSquareNormalizedX(); 
//    Chi2TrackerYInner = trackFit.ChiSquareNormalizedY();
//    Chi2TrackerX = trackFit.ChiSquareNormalizedX();
//    Chi2TrackerY = trackFit.ChiSquareNormalizedY();
  }

//  iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::All, refitPattern, particleHypothesis);
//  if (iFit >= 0) {
//    const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
//    Chi2TrackerY = trackFit.ChiSquareNormalizedY();
//  }

//  iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::UpperHalf, refitPattern, particleHypothesis);
//  if (iFit >= 0) {
//    const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
//    RigidityInnerUpperHalf = trackFit.Rigidity();
//  }

//  iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::LowerHalf, refitPattern, particleHypothesis);
//  if (iFit >= 0) {
//    const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
//    RigidityInnerLowerHalf = trackFit.Rigidity();
//  }

//  iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::UpperHalf, refitPattern, particleHypothesis);
//  if (iFit >= 0) {
//    const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
//    RigidityInnerUpperHalf = trackFit.Rigidity();
//  }

//  if (RigidityInnerL1() + RigidityInnerL9() != 0)
//    {RigidityAsymmetry = (RigidityInnerL1() - RigidityInnerL9()) / (RigidityInnerL1() + RigidityInnerL9());}

//  if (RigidityInnerL9() + RigidityInner() != 0)
//    {RigidityAsymmetryL9 = (RigidityInnerL9() - RigidityInner()) / (RigidityInnerL9() + RigidityInner());}

//  if (Chi2TrackerYInner() + Chi2TrackerY() != 0)
//    {Chi2TrackerYAsymmetry = (Chi2TrackerYInner() - Chi2TrackerY()) / (Chi2TrackerYInner() + Chi2TrackerY());}

//  if (RigidityInner() != 0 && Rigidity() != 0)
//    {InnerMaxSpanRigidityMatching = 100.0 * ((1.0 / RigidityInner()) - (1.0 / Rigidity()));}

//  if (RigidityInnerL1() != 0 && RigidityInnerL9() != 0)
//    {L1L9RigidityMatching = 100.0 * ((1.0 / RigidityInnerL1()) - (1.0 / RigidityInnerL9()));}

//  if (RigidityInnerUpperHalf() != 0 && RigidityInnerLowerHalf() != 0)
//    {L24L58RigidityMatching = 100.0 * ((1.0 / RigidityInnerUpperHalf()) - (1.0 / RigidityInnerLowerHalf()));}

//  if (Chi2TrackerXInner() > 0)
//    {Log10Chi2TrackerXInner = std::log10(Chi2TrackerXInner());}

//  if (Chi2TrackerYInner() > 0)
//    {Log10Chi2TrackerYInner = std::log10(Chi2TrackerYInner());}

//  if (Chi2TrackerX() > 0)
//    {Log10Chi2TrackerX = std::log10(Chi2TrackerX());}

//  if (Chi2TrackerY() > 0)
//    {Log10Chi2TrackerY = std::log10(Chi2TrackerY());}


  assert(TrackerCharges().empty());
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay1));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay2));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay3));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay4));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay5));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay6));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay7));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay8));
  TrackerCharges().emplace_back(primaryParticle->TrackerChargeFor(Analysis::Particle::TrkLay9));

//  const double charge58 = CalculateMeanChargeForLayers(5, 8, TrackerCharges());
//  const double charge24 = CalculateMeanChargeForLayers(2, 4, TrackerCharges());
//  if (charge58 > 0 && charge24 > 0 && InnerTrackerCharge() > 0)
//    {TrackerL58L24ChargeAsymmetry = (charge58 - charge24) / InnerTrackerCharge();}

  if (TrackerCharges().at(8) > 0)
    {TrackerL9Charge = TrackerCharges().at(8);}

  const double charge = CalculateMeanChargeForLayers(7, 8, TrackerCharges());
  if (charge > 0) 
    {TrackerL78Charge = charge;}

//  UpperTofCharge = primaryParticle->UpperTofCharge();
  LowerTofCharge = primaryParticle->LowerTofCharge();
  //// MVA end

    } //  primaryParticle->TrackerTrack }

  TofBeta = primaryParticle->BetaTof();
  EcalEnergyElectron = primaryParticle->EcalEnergyElectron();
  EcalBDT_EnergyD = primaryParticle->EcalEstimator(AC::ECALShower::BDTv7_EnergyD);
  EcalBDT_EnergyD_Smoothed = primaryParticle->EcalEstimator(AC::ECALShower::BDTv7_EnergyD_Smoothed);

  if (primaryParticle->RichRing()) {
    const auto* richRing = primaryParticle->RichRing();
    assert(richRing);
    RichIsNaF = richRing->IsNaF();
    RichBeta = richRing->Beta();
  }

  bool pXeOk = false;
  const double pXe = Utilities::SlowControlLookup::Self()->QueryXenonPressure(primaryParticle->AnaEvent()->TimeStamp(), pXeOk);
  Analysis::TrdLikelihoodCalculation likelihoodCalculation(primaryParticle->TrdHitsFromTrackerTrack(), pXe);
  TrdLogLikelihoodRatioElectronProtonTracker = likelihoodCalculation.ComputeLogLikelihoodRatio(ParticleId::Electron, ParticleId::Proton, primaryParticle->Rigidity());
  TrdLogLikelihoodRatioProtonHeliumTracker = likelihoodCalculation.ComputeLogLikelihoodRatio(ParticleId::Proton, ParticleId::Helium, primaryParticle->Rigidity());

  } //end of PrimaryParticle loop 

 // Proton CC MVA
   fProtonChargeConfusionMva->ProcessEvent(event);
  if (fProtonChargeConfusionMva->IsApplicable(std::abs(Rigidity())))
      ProtonCCMVABDT = fProtonChargeConfusionMva->EvaluateClassifier(fProtonChargeConfusionMva->CategoryForEvent(std::abs(Rigidity())));

 // Electron CC MVA
   fElectronChargeConfusionMva->ProcessEvent(event);
  if (fElectronChargeConfusionMva->IsApplicable(EcalEnergyElectron()))
      ElectronCCMVABDT = fElectronChargeConfusionMva->EvaluateClassifier(fElectronChargeConfusionMva->CategoryForEvent(EcalEnergyElectron()));

} //end of event loop

void AntiprotonIntermediateEnergyTree::UpdateInMemoryBranches() {
}
