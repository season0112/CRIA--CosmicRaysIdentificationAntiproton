#include "AnalysisEvent.hh"
#include "AntiprotonLowEnergyTree.hh"

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
#include <TrdKLikelihoodCalculator.hh>
#define INFO_OUT_TAG "AntiprotonLowEnergyTree"
#include "debugging.hh"

AntiprotonLowEnergyTree::AntiprotonLowEnergyTree()
  : IO::TreeInterface("AntiprotonLowEnergyTree", "Antiproton Low Energy Tree") {
  RegisterBranches();
 // Load proton charge-confusion MVA
   const auto* mvaInterfaceProtonCC = Mva::MvaInterface::CreateMvaInterfaceByName("Mva::ProtonChargeConfusionMvaInterface");
      fProtonChargeConfusionMva = mvaInterfaceProtonCC->CreateMvaImplementation();
 // Load electron charge-confusion MVA
   const auto* mvaInterfaceElectronCC = Mva::MvaInterface::CreateMvaInterfaceByName("Mva::ElectronChargeConfusionMvaInterface");
      fElectronChargeConfusionMva = mvaInterfaceElectronCC->CreateMvaImplementation();
}

void AntiprotonLowEnergyTree::Fill(const Analysis::Event& event) {
    TimeStamp = event.TimeStamp().GetSec();
    Run = event.Run();
    EventNumber = event.EventNumber();
    Weight = event.Weight();
    IsEventMC = event.IsMC();
    TriggerFlags = event.TriggerFlags();
    NumberOfEcalShower = event.NumberOfEcalShower();
    TofNumberOfLayers = event.TofNumberOfLayers();

    if (event.IsMC()) {
      MCPrimaryMomentum = event.McMomentum();
      MCParticleID = event.McParticleId();
    }

    if (const Analysis::Particle* primaryParticle = event.PrimaryParticle()) {  // { event.PrimaryParticle

        // ECAL
        EcalEnergyElectron = primaryParticle->EcalEnergyElectron();
        EcalBDT_EnergyD = primaryParticle->EcalEstimator(AC::ECALShower::BDTv7_EnergyD);
        EcalBDT_EnergyD_Smoothed = primaryParticle->EcalEstimator(AC::ECALShower::BDTv7_EnergyD_Smoothed);

        // TRACKER
        if (const AC::TrackerTrack* trackertrack=primaryParticle->TrackerTrack()){  
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
            int iFit = trackertrack->GetFitFuzzy(AC::Choutko, AC::Inner, refitPattern, particleHypothesis);
            if (iFit >= 0) {
              const AC::TrackerTrackFit& trackFit = trackertrack->TrackFits().at(iFit);
              Rigidity = trackFit.Rigidity();
            }
        
            // "good" charge from X-cluster or Y if no "good" X 
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

            ExtrapolatedRichTileIndex = trackertrack->ExtrapolatedRichTileIndex();
            ExtrapolatedRichExpectedPhotoElectronsProton = trackertrack->ExtrapolatedRichExpectedPhotoElectronsProton();
            ExtrapolatedRichExpectedPhotoElectronsPion = trackertrack->ExtrapolatedRichExpectedPhotoElectronsPion();
            ExtrapolatedRichExpectedPhotoElectronsElectron = trackertrack->ExtrapolatedRichExpectedPhotoElectronsElectron();
            CCLikelihood = trackertrack->CCLikelihood();
            CCBDT = trackertrack->CCBDT();
            CCBDTLapp = trackertrack->CCBDTLapp();

            TrdKCharge = trackertrack->TrdKCharge();
            /* TrdKLikelihood
            Analysis::TrdKLikelihoodCalculator lhcalc_trdk;
            std::vector<Analysis::GenericTrdHit> genericHits;
            lhcalc_trdk.PrepareGenericTrdHitK(*trackertrack,genericHits);
            bool pXeOk = false;
            const double pXe = Utilities::SlowControlLookup::Self()->QueryXenonPressure(primaryParticle->AnaEvent()->TimeStamp(), pXeOk);
            TrdKLikelihood = lhcalc_trdk.ComputeElectronProtonLikelihoodK(pXe,std::abs(primaryParticle->Rigidity()),genericHits);
            */
            assert(DepositedEnergyX().empty());
            assert(DepositedEnergyY().empty()); 
            assert(ClusterWidthX().empty());
            assert(ClusterWidthY().empty());
            for(auto& trkhit: trackertrack->ReconstructedHits()){
                ReconstructedHitLayer().emplace_back(trkhit.Layer());
                DepositedEnergyX().emplace_back(trkhit.DepositedEnergyX());
                DepositedEnergyY().emplace_back(trkhit.DepositedEnergyY());
                ClusterWidthX().emplace_back(trkhit.ClusterWidthX());
                ClusterWidthY().emplace_back(trkhit.ClusterWidthY());
                UnbiasedQX().emplace_back(trkhit.UnbiasedQX());
                UnbiasedQY().emplace_back(trkhit.UnbiasedQY());
                ChargeYiJiaXY().emplace_back(trkhit.ChargeYiJiaXY());
            }
        } 
        
        // RICH
        RICHNumberOfHits = event.RawEvent()->RICH().NumberOfHits();
        RichNumberOfRings = event.NumberOfRichRings(); //same as event.RawEvent()->RICH().Rings().size()
        if (primaryParticle->RichRing()) {
            const auto* richRing = primaryParticle->RichRing();
            RichIsNaF = richRing->IsNaF();
            RichBeta = richRing->Beta(); //primaryParticle->RichBeta()
            RichCharge = richRing->ChargeEstimate();
            RICHRINGNumberOfHits = richRing->NumberOfHits();      
            BetaConsistency = richRing->BetaConsistency(); //Beta consistency is the difference between the reconstructed beta of the two reconstruction methods from LIP and CIEMAT. (if available) 
            TileIndex = richRing->TileIndex(); //Aerogel:[0,120], NaF:(128..143), outside of the radiator area:255, default value:254
            NExpectedPhotoElectrons = richRing->NExpectedPhotoElectrons(); /// Number of expected photoelectrons for a Z=1 ring with the reconstruction input parameters of the current event. 
            NPhotoElectrons = richRing->NPhotoElectrons(); //in this ring. 
            NCollectedPhotoElectrons= richRing->NCollectedPhotoElectrons(); //in this event. 
            IsGood = richRing->IsGood();
        }

        // TOF
        UpperTofCharge = primaryParticle->UpperTofCharge();
        LowerTofCharge = primaryParticle->LowerTofCharge();
        TofBeta = primaryParticle->BetaTof();
        TofBetaSize = event.RawEvent()->TOF().Betas().size();
        // BetaH
        BetaConverted = primaryParticle->TofBeta()->BetaConverted();
        TofMassonecharge = std::abs(primaryParticle->Rigidity()) * std::sqrt(pow(1.0 / primaryParticle->TofBeta()->BetaConverted(), 2) - 1.0);
        assert(BetaFromDeDx().empty());
        assert(TOFBETACharges().empty());
        for (unsigned int i = 0; i < primaryParticle->TofBeta()->BetaFromDeDx().size(); ++i){
            BetaFromDeDx().emplace_back(primaryParticle->TofBeta()->BetaFromDeDx().at(i));
        }
        if (BetaFromDeDx().size()==8){
        if ((BetaFromDeDx().at(0)!=1.0) && (BetaFromDeDx().at(2)!=1.0))
            UpperTofBeta = (BetaFromDeDx().at(0)+BetaFromDeDx().at(2))/2;
        if ((BetaFromDeDx().at(0)!=1.0) && (BetaFromDeDx().at(2)==1.0))
            UpperTofBeta = BetaFromDeDx().at(0);
        if ((BetaFromDeDx().at(0)==1.0) && (BetaFromDeDx().at(2)!=1.0))
            UpperTofBeta = BetaFromDeDx().at(2);
        if ((BetaFromDeDx().at(4)!=1.0) && (BetaFromDeDx().at(6)!=1.0))
            LowerTofBeta = (BetaFromDeDx().at(4)+BetaFromDeDx().at(6))/2;
        if ((BetaFromDeDx().at(4)!=1.0) && (BetaFromDeDx().at(6)==1.0))
            LowerTofBeta = BetaFromDeDx().at(4);
        if ((BetaFromDeDx().at(4)==1.0) && (BetaFromDeDx().at(6)!=1.0))
            LowerTofBeta = BetaFromDeDx().at(6);
        }
        for (unsigned int i = 0; i < primaryParticle->TofBeta()->Charges().size(); ++i){
            TOFBETACharges().emplace_back(primaryParticle->TofBeta()->Charges().at(i));  // charge with beta + rigidity correction
            TOFClusterEnergy().emplace_back(primaryParticle->AnaEvent()->RawEvent()->TOF().Clusters().at( primaryParticle->TofBeta()->TOFClusterIndex().at(i) ).EnergyInMeV());
            TOFClusterCharge().emplace_back(primaryParticle->AnaEvent()->RawEvent()->TOF().Clusters().at( primaryParticle->TofBeta()->TOFClusterIndex().at(i) ).Charge()); //Charge estimate from all good PMTs without path length and beta correction
        }

        // TRD
        //const Analysis::Particle::TrdHitsVector& trdHits = primaryParticle->TrdHitsFromTrackerTrack();
        //assert(dEdX_TRD().empty());
        //for (unsigned int i = 0; i < trdHits.size(); ++i) {
        //    const Analysis::TrdHit& trdhit = trdHits.at(i);
        //    dEdX_TRD().emplace_back(trdhit.GetDeDx());   //trdHits.size could large than 20.
        //}
        TrdNumberOfHits = event.NumberOfTrdRawHits();
        TRDVTracksSize = event.RawEvent()->TRD().VTracks().size();
        TRDHTracksSize = event.RawEvent()->TRD().HTracks().size();
        TrdSegmentsXZNumber =  event.TrdSegmentsXZ().size();
        TrdSegmentsYZNumber =  event.TrdSegmentsYZ().size();
        TrdVerticesXZNumber =  event.TrdVerticesXZ().size();
        TrdVerticesYZNumber =  event.TrdVerticesYZ().size();

        const Analysis::TrdTrack* trdTrack = primaryParticle->GetTrdTrack();
        NumberOfHitsOnTrack = trdTrack->NumberOfHitsOnTrack();
        NumberOfLayersWithHit = trdTrack->NumberOfLayersWithHit();
        TrdTrackNumberOfSubLayersXZ = trdTrack->SegmentXZ()->NumberOfSublayersInSegment();
        TrdTrackNumberOfSubLayersYZ = trdTrack->SegmentYZ()->NumberOfSublayersInSegment();

        bool pXeOk = false;
        const double pXe = Utilities::SlowControlLookup::Self()->QueryXenonPressure(primaryParticle->AnaEvent()->TimeStamp(), pXeOk);

        Analysis::TrdLikelihoodCalculation likelihoodCalculation(primaryParticle->TrdHitsFromTrackerTrack(), pXe);
        TrdLogLikelihoodRatioElectronProtonTracker = likelihoodCalculation.ComputeLogLikelihoodRatio(ParticleId::Electron, ParticleId::Proton, primaryParticle->Rigidity());
        TrdLogLikelihoodRatioProtonHeliumTracker = likelihoodCalculation.ComputeLogLikelihoodRatio(ParticleId::Proton, ParticleId::Helium, primaryParticle->Rigidity());
        TrdActiveLayersTracker = likelihoodCalculation.NumberOfActiveLayers();

        Analysis::TrdLikelihoodCalculation likelihoodCalculatorStandalone(primaryParticle->TrdHitsFromTrdTrack(), pXe);
        TrdActiveLayersStandalone = likelihoodCalculatorStandalone.NumberOfActiveLayers();

        Analysis::TrdLikelihoodCalculation likelihoodCalculatorHybrid(primaryParticle->TrdHitsFromTrdAndTrackerTracks(), pXe);
        TrdActiveLayersHybrid = likelihoodCalculatorHybrid.NumberOfActiveLayers();

        // ACC
        ACCHits = event.RawEvent()->ACC().Clusters().size();
        NumberOfClustersMIT = event.RawEvent()->ACC().NumberOfClustersMIT();

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

void AntiprotonLowEnergyTree::UpdateInMemoryBranches() {
}
