#include "AnalysisEvent.hh"
#include "RigidityResolutionTree.hh"
#include "AMSGeometry.h"
#include "TrackerTrack.h"
#include "MvaImplementation.hh"
#include "MvaInterface.hh"
#include "Tracker.h"
#include "TrackerTrackFit.h"
#include "Event.h"
#include <TrdLikelihoodCalculation.hh>
#include <cmath>
#include <AnalysisEvent.hh>
#define INFO_OUT_TAG "RigidityResolutionTree"
#include "debugging.hh"
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

RigidityResolutionTree::RigidityResolutionTree()
  : IO::TreeInterface("RigidityResolutionTree", "Example analysis tree") {

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


void RigidityResolutionTree::Fill(const Analysis::Event& event) { /////////// { Analysis::Event
//    ZTrackerLayer=ZTrackerLayer();
        ZTrackerLayer1=AC::AMSGeometry::ZTrackerLayer1;
        ZTrackerLayer2=AC::AMSGeometry::ZTrackerLayer2;
        ZTrackerLayer3=AC::AMSGeometry::ZTrackerLayer3;
        ZTrackerLayer4=AC::AMSGeometry::ZTrackerLayer4;
        ZTrackerLayer5=AC::AMSGeometry::ZTrackerLayer5;
        ZTrackerLayer6=AC::AMSGeometry::ZTrackerLayer6;
        ZTrackerLayer7=AC::AMSGeometry::ZTrackerLayer7;
        ZTrackerLayer8=AC::AMSGeometry::ZTrackerLayer8;
        ZTrackerLayer9=AC::AMSGeometry::ZTrackerLayer9;
        assert(UnassociatedHitZLayer().empty());
        assert(UnassociatedHitXCoordinate().empty());
        assert(UnassociatedHitYCoordinate().empty());
        assert(UnassociatedChargeX().empty());
        assert(UnassociatedChargeY().empty());
  for(auto & unassociatedhit: event.RawEvent()->Tracker().UnassociatedReconstructedHits()){
    UnassociatedHitZLayer().emplace_back(unassociatedhit.Layer());
    UnassociatedHitXCoordinate().emplace_back(unassociatedhit.X());
    UnassociatedHitYCoordinate().emplace_back(unassociatedhit.Y());
    UnassociatedChargeX().emplace_back(unassociatedhit.QX());
    UnassociatedChargeY().emplace_back(unassociatedhit.QY());
    }

  TimeStamp = event.TimeStamp().GetSec();
  Run = event.Run();
  EventNumber = event.EventNumber();
  Weight = event.Weight();
  IsEventMC = event.IsMC();
  if (event.IsMC()) {
    MCPrimaryMomentum = event.McMomentum();
    MCParticleID = event.McParticleId();
  }
  if (const Analysis::Particle* primaryParticle = event.PrimaryParticle()) {  /////////// { event.PrimaryParticle

      EcalEnergyElectron = primaryParticle->EcalEnergyElectron();
/*      
  if (const AC::ECALShower* shower = primaryParticle->EcalShower()){
          EcalBDT_EnergyD = shower->EcalBDTv7_EnergyD();
          EcalBDT_EnergyD_Smoothed = shower->EcalBDTv7_EnergyD_Smoothed();
      }
*/
  if (primaryParticle && primaryParticle->HasEcalShower()) { /////////// { event.HasEcalShower
    const AC::ECALShower* shower = primaryParticle->EcalShower();
    assert(shower);
    EcalBDT_EnergyD = primaryParticle->EcalEstimator(AC::ECALShower::BDTv7_EnergyD);
    EcalBDT_EnergyD_Smoothed = primaryParticle->EcalEstimator(AC::ECALShower::BDTv7_EnergyD_Smoothed);
  }    ///////////  event.HasEcalShower }


    InnerTrackerCharge = primaryParticle->TrackerCharge();
    // Example on how to fill a vector.
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

      if (const AC::TrackerTrack* trackertrack=primaryParticle->TrackerTrack()){  /////////// { primaryParticle->TrackerTrack
        assert(HitXCoordinate().empty());
        assert(HitYCoordinate().empty());
        assert(HitZLayer().empty());
        assert(ChargeX().empty());
        assert(ChargeY().empty());
        assert(RigidityWithoutThisHit().empty());
        assert(ClusterWidthX().empty());
        assert(ClusterWidthY().empty());
        assert(ResidualYInMicrons().empty());
          for(auto& hit: trackertrack->ReconstructedHits()){
    //fill XYZ and ChargeXY Coordinate.
            HitZLayer().emplace_back(hit.Layer());
            HitXCoordinate().emplace_back(hit.X());
            HitYCoordinate().emplace_back(hit.Y());
            ChargeX().emplace_back(hit.QX());
            ChargeY().emplace_back(hit.QY());
            ClusterWidthX().emplace_back(hit.ClusterWidthX());
            ClusterWidthY().emplace_back(hit.ClusterWidthY());
            RigidityWithoutThisHit().emplace_back(hit.RigidityWithoutThisHit());
            ResidualYInMicrons().emplace_back(hit.ResidualYInMicrons());
          }
  //TrackerLayerPatternClassification
          CCBDT=trackertrack->CCBDT();
          CCBDTLapp=trackertrack->CCBDTLapp();
          CCLikelihood=trackertrack->CCLikelihood();
          Pattern=trackertrack->TrackerLayerPatternClassification();


    const auto refitPattern2 = AC::PGMA + AC::RebuildFromTDV;
    const auto particleHypothesis2 = AC::DefaultMass;

//    const size_t fitPatterns = 8;
    const size_t fitAlgorithms = 8;
/*
    RigidityValue().resize(fitPatterns);
    RigidityInverseError().resize(fitPatterns);
    RigidityChi2X().resize(fitPatterns);
    RigidityChi2Y().resize(fitPatterns);
    RigidityResidualYL1().resize(fitPatterns);
    RigidityResidualYL9().resize(fitPatterns);
*/

    RigidityValue().resize(fitAlgorithms);
    RigidityInverseError().resize(fitAlgorithms);
    RigidityChi2X().resize(fitAlgorithms);
    RigidityChi2Y().resize(fitAlgorithms);
    RigidityResidualYL1().resize(fitAlgorithms);
    RigidityResidualYL9().resize(fitAlgorithms);

/*
    for (auto fitPattern : {
        AC::All, AC::UpperHalf, AC::LowerHalf,
        AC::Inner, AC::TwoInnerPlusExt,
        AC::InnerPlusL1, AC::InnerPlusL9, AC::InnerPlusL1L9 }) {
      int fitIndex = trackertrack->GetFitFuzzy(
        AC::Choutko, fitPattern, refitPattern2, particleHypothesis2);
      if (fitIndex >= 0) {
        const auto& trackFit2 = trackertrack->TrackFits().at(fitIndex);
        RigidityValue().at(fitPattern) = trackFit2.Rigidity();
        RigidityInverseError().at(fitPattern) = trackFit2.InverseRigidityError();
        RigidityChi2X().at(fitPattern) = trackFit2.ChiSquareNormalizedX();
        RigidityChi2Y().at(fitPattern) = trackFit2.ChiSquareNormalizedY();
        RigidityResidualYL1().at(fitPattern) = trackFit2.YResidualJLayer1();
        RigidityResidualYL9().at(fitPattern) = trackFit2.YResidualJLayer9();
      }
    }
*/

/*
    assert(RigidityValue().empty());
    assert(RigidityInverseError().empty());
    assert(RigidityChi2X().empty());
    assert(RigidityChi2Y().empty());
    assert(RigidityResidualYL1().empty());
    assert(RigidityResidualYL9().empty());
*/

    for (auto fitAlgorithm : {
        AC::Default, AC::Choutko, AC::Alcaraz, AC::ChikanianF, AC::ChikanianC, AC::Kalman, AC::VertexFit}) {
      int fitIndex = trackertrack->GetFitFuzzy(
       fitAlgorithm , AC::All, refitPattern2, particleHypothesis2);
      if (fitIndex >= 0) {
        const auto& trackFit2 = trackertrack->TrackFits().at(fitIndex);
        RigidityValue().at(fitAlgorithm) = trackFit2.Rigidity();
        RigidityInverseError().at(fitAlgorithm) = trackFit2.InverseRigidityError();
        RigidityChi2X().at(fitAlgorithm) = trackFit2.ChiSquareNormalizedX();
        RigidityChi2Y().at(fitAlgorithm) = trackFit2.ChiSquareNormalizedY();
        RigidityResidualYL1().at(fitAlgorithm) = trackFit2.YResidualJLayer1();
        RigidityResidualYL9().at(fitAlgorithm) = trackFit2.YResidualJLayer9();
      }
    }


        }  ///////////  primaryParticle->TrackerTrack }

     if (const AC::TrackerTrackCoordinates* trackertrackcoordinates=primaryParticle->TrackerTrackCoordinates()){ /////////// { primaryParticle->TrackerTrackCoordinates
          assert(XLayer().empty());
          assert(YLayer().empty());
          XLayer().emplace_back(trackertrackcoordinates->XLayer1());
          XLayer().emplace_back(trackertrackcoordinates->XLayer2());
          XLayer().emplace_back(trackertrackcoordinates->XLayer3());
          XLayer().emplace_back(trackertrackcoordinates->XLayer4());
          XLayer().emplace_back(trackertrackcoordinates->XLayer5());
          XLayer().emplace_back(trackertrackcoordinates->XLayer6());
          XLayer().emplace_back(trackertrackcoordinates->XLayer7());
          XLayer().emplace_back(trackertrackcoordinates->XLayer8());
          XLayer().emplace_back(trackertrackcoordinates->XLayer9());
          YLayer().emplace_back(trackertrackcoordinates->XLayer1());
          YLayer().emplace_back(trackertrackcoordinates->XLayer2());
          YLayer().emplace_back(trackertrackcoordinates->XLayer3());
          YLayer().emplace_back(trackertrackcoordinates->XLayer4());
          YLayer().emplace_back(trackertrackcoordinates->XLayer5());
          YLayer().emplace_back(trackertrackcoordinates->XLayer6());
          YLayer().emplace_back(trackertrackcoordinates->XLayer7());
          YLayer().emplace_back(trackertrackcoordinates->XLayer8());
          YLayer().emplace_back(trackertrackcoordinates->XLayer9());
        }   //////////////////// primaryParticle->TrackerTrackCoordinates }


  } /////////// event.PrimaryParticle } 

 // 16 MVA varabiles
  const AC::Event* rawEvent = event.RawEvent();
  assert(rawEvent);
 
  const Analysis::Particle* particle = event.PrimaryParticle();
  assert(particle);

  const AC::TrackerTrack* trackerTrack = particle->TrackerTrack();
  assert(trackerTrack);  

  const int refitPattern = AC::PGMA + AC::RebuildFromTDV;
  const AC::ParticleHypothesis particleHypothesis = AC::DefaultMass;

  int iFit = trackerTrack->GetFitFuzzy(AC::Choutko, AC::All, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    Rigidity = trackFit.Rigidity();
    InverseRigidityError = trackFit.InverseRigidityError();
    Chi2TrackerX = trackFit.ChiSquareNormalizedX();
    Chi2TrackerY = trackFit.ChiSquareNormalizedY();
    YResidualJLayer1 = trackFit.YResidualJLayer1();
    YResidualJLayer9 = trackFit.YResidualJLayer9();
  }

  iFit = trackerTrack->GetFitFuzzy(AC::Choutko, AC::Inner, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    RigidityInner = trackFit.Rigidity();
    Chi2TrackerXInner = trackFit.ChiSquareNormalizedX();
    Chi2TrackerYInner = trackFit.ChiSquareNormalizedY();
    InverseRigidityErrorInner = trackFit.InverseRigidityError();
  }

  iFit = trackerTrack->GetFitFuzzy(AC::Alcaraz, AC::All, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    RigidityAlcaraz = trackFit.Rigidity();
    InverseRigidityErrorAlcaraz = trackFit.InverseRigidityError();
  }

  iFit = trackerTrack->GetFitFuzzy(AC::ChikanianF, AC::All, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    RigidityChikanian = trackFit.Rigidity();
    Chi2TrackerYChikanian = trackFit.ChiSquareNormalizedY();
    InverseRigidityErrorChikanian = trackFit.InverseRigidityError();
  }

  iFit = trackerTrack->GetFitFuzzy(AC::Choutko, AC::InnerPlusL1, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    RigidityInnerL1 = trackFit.Rigidity();
    Chi2TrackerXInnerL1 = trackFit.ChiSquareNormalizedX();
    Chi2TrackerYInnerL1 = trackFit.ChiSquareNormalizedY();
  }

  iFit = trackerTrack->GetFitFuzzy(AC::Choutko, AC::InnerPlusL9, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    RigidityInnerL9 = trackFit.Rigidity();
    Chi2TrackerXInnerL9 = trackFit.ChiSquareNormalizedX();
    Chi2TrackerYInnerL9 = trackFit.ChiSquareNormalizedY();
  }


    iFit = trackerTrack->GetFitFuzzy(AC::Choutko, AC::UpperHalf, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    RigidityInnerUpperHalf = trackFit.Rigidity();
    Chi2TrackerXInnerUpperHalf = trackFit.ChiSquareNormalizedX();
    Chi2TrackerYInnerUpperHalf = trackFit.ChiSquareNormalizedY();
  }

  iFit = trackerTrack->GetFitFuzzy(AC::Choutko, AC::LowerHalf, refitPattern, particleHypothesis);
  if (iFit >= 0) {
    const AC::TrackerTrackFit& trackFit = trackerTrack->TrackFits().at(iFit);
    RigidityInnerLowerHalf = trackFit.Rigidity();
    Chi2TrackerXInnerLowerHalf = trackFit.ChiSquareNormalizedX();
    Chi2TrackerYInnerLowerHalf = trackFit.ChiSquareNormalizedY();

  }


  if (RigidityInnerL1() + RigidityInnerL9() != 0)
  {  RigidityAsymmetry = (RigidityInnerL1() - RigidityInnerL9()) / (RigidityInnerL1() + RigidityInnerL9());}
  if (RigidityInnerL9() + RigidityInner() != 0)
  {  RigidityAsymmetryL9 = (RigidityInnerL9() - RigidityInner()) / (RigidityInnerL9() + RigidityInner());}
  if (Chi2TrackerYInner() + Chi2TrackerY() != 0)
  {  Chi2TrackerYAsymmetry = (Chi2TrackerYInner() - Chi2TrackerY()) / (Chi2TrackerYInner() + Chi2TrackerY());}

  if (RigidityInner() != 0 && Rigidity() != 0)
  {  InnerMaxSpanRigidityMatching = 100.0 * ((1.0 / RigidityInner()) - (1.0 / Rigidity()));
    InnerMaxSpanRigidityMatchingv2 = 100.0 * RigidityInner() * ((1.0 / RigidityInner()) - (1.0 / Rigidity()));}

  if (RigidityInnerL1() != 0 && RigidityInnerL9() != 0)
  {  L1L9RigidityMatching = 100.0 * ((1.0 / RigidityInnerL1()) - (1.0 / RigidityInnerL9()));
    L1L9RigidityMatchingv2 = 100.0 * RigidityInnerL1() * ((1.0 / RigidityInnerL1()) - (1.0 / RigidityInnerL9()));}

  if (RigidityInnerUpperHalf() != 0 && RigidityInnerLowerHalf() != 0)
  {  L24L58RigidityMatching = 100.0 * ((1.0 / RigidityInnerUpperHalf()) - (1.0 / RigidityInnerLowerHalf()));
    L24L58RigidityMatchingv2 = 100.0 * RigidityInnerUpperHalf() * ((1.0 / RigidityInnerUpperHalf()) - (1.0 / RigidityInnerLowerHalf()));}

  if (Chi2TrackerXInner() > 0)
    {Log10Chi2TrackerXInner = std::log10(Chi2TrackerXInner());}

  if (Chi2TrackerYInner() > 0)
    {Log10Chi2TrackerYInner = std::log10(Chi2TrackerYInner());}

  if (Chi2TrackerX() > 0)
    {Log10Chi2TrackerX = std::log10(Chi2TrackerX());}

  if (Chi2TrackerY() > 0)
    {Log10Chi2TrackerY = std::log10(Chi2TrackerY());}


  {
    const double charge58 = CalculateMeanChargeForLayers(5, 8, TrackerCharges());
    const double charge24 = CalculateMeanChargeForLayers(2, 4, TrackerCharges());
    if (charge58 > 0 && charge24 > 0 && InnerTrackerCharge() > 0)
      TrackerL58L24ChargeAsymmetry = (charge58 - charge24) / InnerTrackerCharge();
  }


  if (TrackerCharges().at(8) > 0)
    TrackerL9Charge = TrackerCharges().at(8);

  {
    const double charge = CalculateMeanChargeForLayers(7, 8, TrackerCharges());
    if (charge > 0) {
      TrackerL78Charge = charge;
    }
  }


  UpperTofCharge = particle->UpperTofCharge();
  LowerTofCharge = particle->LowerTofCharge();
  RichBeta = particle->RichBeta();

  bool pXeOk = false;
  const double pXe = Utilities::SlowControlLookup::Self()->QueryXenonPressure(particle->AnaEvent()->TimeStamp(), pXeOk);
  Analysis::TrdLikelihoodCalculation likelihoodCalculation(particle->TrdHitsFromTrackerTrack(), pXe);

  TrdLogLikelihoodRatioElectronProtonTracker = likelihoodCalculation.ComputeLogLikelihoodRatio(ParticleId::Electron, ParticleId::Proton, particle->Rigidity());
  TrdLogLikelihoodRatioProtonHeliumTracker = likelihoodCalculation.ComputeLogLikelihoodRatio(ParticleId::Proton, ParticleId::Helium, particle->Rigidity());


 // Proton CC MVA
  fProtonChargeConfusionMva->ProcessEvent(event);
  if (fProtonChargeConfusionMva->IsApplicable(std::abs(Rigidity())))
      ProtonCCMVABDT = fProtonChargeConfusionMva->EvaluateClassifier(fProtonChargeConfusionMva->CategoryForEvent(std::abs(Rigidity())));
 
 // Electron CC MVA
  fElectronChargeConfusionMva->ProcessEvent(event);
  if (fElectronChargeConfusionMva->IsApplicable(EcalEnergyElectron()))
      ElectronCCMVABDT = fElectronChargeConfusionMva->EvaluateClassifier(fElectronChargeConfusionMva->CategoryForEvent(EcalEnergyElectron()));
  
 
} ///////////  Analysis::Event }

void RigidityResolutionTree::UpdateInMemoryBranches() {

}
