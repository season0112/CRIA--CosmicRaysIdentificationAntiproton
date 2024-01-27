#include "ExampleAnalysisTree.hh"

// ACsoft includes
#include "AnalysisEvent.hh"
#include "AcceptanceManager.hh"
#include "AntiprotonBinning.hh"
#include "BinningFunctions.hh"
#include "BinningDefinition.hh"
#include "BinningTools.hh"
#include "ConfigHandler.hh"
#include "CutAttachment.hh"
#include "CutFactory.hh"
#include "EventFactory.hh"
#include "EfficiencyHistograms.hh"
#include "Environment.hh"
#include "FileManager.hh"
#include "McSpectrumScaler.hh"
#include "ObjectManager.hh"
#include "PredefinedBinnings.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "SlowControlLookup.hh"
#include "TreeWriter.hh"
#include "ValueHistograms.hh"

// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include "TreeFormula.hh"
#include <TROOT.h>
#include "TreeWriter.hh"
#include <TApplication.h>
#include <TProof.h>
#include <TAxis.h>

#include "debugging.hh"


  REGISTER_CUT(MyNegativeRigidity,
             "Rigidity < 0",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               if (!p->HasTrackerTrackFit())
                 return true;
               valueForHistograms = p->Rigidity();
               return valueForHistograms < 0;
             },
             &Binning::Predefined::RigidityBinning)


  REGISTER_CUT(MyPositiveRigidity,
             "Rigidity > 0",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               if (!p->HasTrackerTrackFit())
                 return true;
               valueForHistograms = p->Rigidity();
               return valueForHistograms > 0;
             },
             &Binning::Predefined::RigidityBinning)


  REGISTER_CUT(MyPattern0,
             "Tracker Pattern = 0",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               valueForHistograms = (p->TrackerLayerPatternClassification() == 0);
               return valueForHistograms;
             })


  REGISTER_CUT(MyPattern1,
             "Tracker Pattern = 1",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               valueForHistograms = (p->TrackerLayerPatternClassification() == 1);
               return valueForHistograms;
             })


  REGISTER_CUT(MyPattern2,
             "Tracker Pattern = 2",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               valueForHistograms = (p->TrackerLayerPatternClassification() == 2);
               return valueForHistograms;
             })


  REGISTER_CUT(MyPattern3,
             "Tracker Pattern = 3",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               valueForHistograms = (p->TrackerLayerPatternClassification() == 3);
               return valueForHistograms;
             })


  REGISTER_CUT(MyPattern4,
             "Tracker Pattern = 4",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               valueForHistograms = (p->TrackerLayerPatternClassification() == 4);
               return valueForHistograms;
             })


  REGISTER_CUT(MyPattern5,
             "Tracker Pattern = 5",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               valueForHistograms = (p->TrackerLayerPatternClassification() == 5);
               return valueForHistograms;
             })


  REGISTER_CUT(MyPatternMinus1,
             "Tracker Pattern = -1",
             [](const Analysis::Event& ev, double& valueForHistograms) -> bool {
               const Analysis::Particle* p = ev.PrimaryParticle();
               if (!p)
                 return true;
               if (!p->HasTrackerTrack())
                 return true;
               valueForHistograms = (p->TrackerLayerPatternClassification() == -1);
               return valueForHistograms;
             })


