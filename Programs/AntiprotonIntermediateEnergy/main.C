#include "AntiprotonIntermediateEnergyTree.hh"

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

#define INFO_OUT_TAG "AntiprotonIntermediateEnergy"
#include "debugging.hh"

int main(int argc, char** argv) {

  // Workaround to avoid ROOT option parsing
//  static int sNull = 0;
//  TApplication* fApp = new TApplication("Application", &sNull, (char**)0);

  // Command line option handling.
  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);

  config.SetProgramHelpText("AntiprotonIntermediateEnergy",
                            "Illustrates the usage of IO::TreeWriter to write ROOT trees from ACQt files.");

  config.AddHelpExample("Loop over given filelist.", "--inputlist list.txt");

  std::string inputList;
  config.GetValue("OPTIONS", "inputlist", inputList,
                  "List of ACQt input files (full path).");

  std::string resultDirectory;
  config.GetValue("OPTIONS", "resultdir", resultDirectory,
                  "General directory where result files should be stored. Current directory is used if option is not specified.");

  std::string suffix;
  config.GetValue("OPTIONS", "suffix", suffix,
                  "A string identifier to be used in parallel computing, to uniquely identify result files.");
  
  std::string particleID = "";
  config.GetValue("OPTIONS", "particleID", particleID, 
                  "The particleID of the MC dataset is");

  unsigned int maxEntries = 0;
  config.GetValue("OPTIONS","entries", maxEntries,
                  "Number of events to process");

  // Load & parse cut configuration file.
  std::string cutConfigfile = "${MY_ANALYSIS}/Configuration/AntiprotonIntermediateEnergy.cfg";
  Environment::ExpandEnvironmentVariables(cutConfigfile);
  config.Read(cutConfigfile);
  Cuts::SelectionParser selectionParser(config);
 
  // Setup file manager to process ACQt data.
  IO::FileManager fileManager(&config);

  Analysis::EventFactory* eventFactory = Analysis::EventFactory::Create(&config);
  Analysis::Event event;

  // McSpectrumScaler for MC event weights
  Utilities::McSpectrumScaler scaler(&config, resultDirectory, suffix);
  scaler.SetDefaultTargetSpectra();
  eventFactory->RegisterMcSpectrumScaler(&scaler);

  // Load geometry file for the acceptance manager.
  Acceptance::AcceptanceManager acceptanceManager;
  std::string geometryConfigFile = "${MY_ANALYSIS}/Configuration/AntiprotonGeometry.cfg";
  Environment::ExpandEnvironmentVariables(geometryConfigFile);
  eventFactory->RegisterAcceptanceManager(&acceptanceManager);

  // 'AuxiliaryObjectManager' holds all auxiliary histograms / selectors created while processing the ACQt files.
  // NOTE: You should NOT write a TTree together with other histograms etc. in the ROOT file. You most likely
  // want to merge your histograms / selectors from batch jobs, but not the trees. That's why it's a good idea
  // in general to split up in two files: one for holding the tree, one for the rest.
  Utilities::ObjectManager auxiliaryObjectManager("AuxiliaryObjectManager", &config, resultDirectory, suffix);
  auxiliaryObjectManager.SetPrefix("AntiprotonIntermediateEnergy_Auxiliary");

  const Binning::Definition& antiprotonbinning = AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili();
  ParticleId::Species particleId = ParticleId::NoSpecies;
  if (particleID == "") {
    WARN_OUT << "No particleID is given! Please give a particleID." << std::endl;
    return EXIT_FAIL_CONFIG;
  }
  if (particleID != "ISS"){
      if (particleID == "Proton"){
          particleId = ParticleId::Proton;
      }
      if (particleID == "Antiproton"){
          particleId = ParticleId::Antiproton;
      }
      if (particleID == "Electron"){
          particleId = ParticleId::Electron;
      }
      if (particleID == "Deuteron"){
          particleId = ParticleId::Deuteron;
      }
      TH1D* generatedHistogram = new TH1D("GeneratedEvents", "Generated events", antiprotonbinning.NumberOfBins(), &antiprotonbinning.Bins()[0]);
      scaler.FillHistogram(particleId, generatedHistogram, Utilities::KinematicVariable::Momentum);
      auxiliaryObjectManager.Add(generatedHistogram, "McGeneratedSpectrum");
  }

  if (!config.PerformChecksAfterOptionParsing())
    return EXIT_FAIL_CONFIG;

  if (!fileManager.ReadFileList(inputList))
    return EXIT_FAIL_FILEMANAGER;

  // Construct tree manager which will manage the output file to hold the resulting tree.
  IO::TreeWriter treeWriter(new AntiprotonIntermediateEnergyTree, IO::TreeOptions::DontWriteInMemoryBranches);
  std::string treeFileName = Utilities::ObjectManager::MakeStandardRootFileName(resultDirectory, "AntiprotonIntermediateEnergy_Tree", suffix);
  treeWriter.Initialize(treeFileName);

  // Load cut selector(s).
  Cuts::Selector* mcPreselection = auxiliaryObjectManager.Add(selectionParser.GetSelector("McPreselection"));
  Cuts::Selector* badRuns = auxiliaryObjectManager.Add(selectionParser.GetSelector("BadRuns"));
  Cuts::Selector* rTI = auxiliaryObjectManager.Add(selectionParser.GetSelector("RTI"));
  Cuts::Selector* preselectionSelector = auxiliaryObjectManager.Add(selectionParser.GetSelector("Preselection"));
  Cuts::Selector* trackerCuts= auxiliaryObjectManager.Add(selectionParser.GetSelector("TrackerCuts"));
  Cuts::Selector* chargeOne = auxiliaryObjectManager.Add(selectionParser.GetSelector("ChargeOne"));
  Cuts::Selector* qualityCuts = auxiliaryObjectManager.Add(selectionParser.GetSelector("QualityCuts"));

  // In order to fill any cut value histograms, the x axis value must be known.
  // We define this in a generic way using a lamdba function and pass it on
  // to Selector::SetupCommonXAxisInformation() for each selector.

  auto xAxisValue = [](const Analysis::Event& event, double& axisValue) {
//    const auto& ecal = event.RawEvent()->ECAL();
//    axisValue = ecal.TotalEnergy3D();  
    const Analysis::Particle* particle = event.PrimaryParticle();
    if (particle)
      axisValue = std::abs(particle->Rigidity());
  };

  auto xAxisValueMc = [](const Analysis::Event& event, double& axisValueMC) {
    if (event.McMomentum() != 0)
    axisValueMC = std::abs(event.McMomentum());    
  };

  badRuns->SetupCommonXAxisInformation(xAxisValue, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  rTI->SetupCommonXAxisInformation(xAxisValue, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  mcPreselection->SetupCommonXAxisInformation(xAxisValue, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  preselectionSelector->SetupCommonXAxisInformation(xAxisValue, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  trackerCuts->SetupCommonXAxisInformation(xAxisValue, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  chargeOne->SetupCommonXAxisInformation(xAxisValue, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  qualityCuts->SetupCommonXAxisInformation(xAxisValue, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);

  badRuns->SetupCommonXAxisInformationMc(xAxisValueMc, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  rTI->SetupCommonXAxisInformationMc(xAxisValueMc, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  mcPreselection->SetupCommonXAxisInformationMc(xAxisValueMc, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  preselectionSelector->SetupCommonXAxisInformationMc(xAxisValueMc, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  trackerCuts->SetupCommonXAxisInformationMc(xAxisValueMc, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  chargeOne->SetupCommonXAxisInformationMc(xAxisValueMc, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);
  qualityCuts->SetupCommonXAxisInformationMc(xAxisValueMc, "Rigidity / GV", AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili);


  // Begin event loop.
  INFO_OUT_ON_MASTER << "Looping over " << fileManager.GetEntries() << " events..." << std::endl;

  static int sProductionSteps = Analysis::CreateSplineTrack | Analysis::CreateTrdTrack;
  unsigned int nEntries=0; 
  bool firstEvent = true;

  while (fileManager.GetNextEvent()) {
    // Initialize the AcceptanceManager, after the DetectorManager received the run type (needed for MPI).
    if (firstEvent) {
      acceptanceManager.InitSetup(geometryConfigFile);
      firstEvent = false;
    }

    fileManager.DumpEventLoopProgress(20000);

    eventFactory->SetupEmptyEvent(event);
    eventFactory->CreateParticles(event);

    if (!mcPreselection->Passes(event))
      continue;
    if (!badRuns->Passes(event))
      continue;
    if (!rTI->Passes(event))
      continue;
    if (!preselectionSelector->Passes(event))
      continue;
    if (!trackerCuts->Passes(event))
      continue;
    if (!chargeOne->Passes(event))
      continue;

    eventFactory->PerformTrdTracking(event);
    eventFactory->FillParticles(event, sProductionSteps);

    if (!qualityCuts->Passes(event))
      continue;

    if (maxEntries>0 && nEntries==maxEntries)
      break;
    ++nEntries;

    treeWriter.Fill(event);
  }

  // Print Preselection statistics
  mcPreselection->PrintSummary();
  badRuns->PrintSummary();
  rTI->PrintSummary();
  preselectionSelector->PrintSummary();
  trackerCuts->PrintSummary();
  chargeOne->PrintSummary();
  qualityCuts->PrintSummary();
  // Finish writing tree file.
  treeWriter.Finish();

  // Write auxiliary output file.
  auxiliaryObjectManager.WriteToFile();

  return EXIT_SUCCESS;
}
