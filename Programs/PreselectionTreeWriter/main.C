#include "ExampleAnalysisTree.hh"

// ACsoft includes
#include "AnalysisEvent.hh"
#include "ConfigHandler.hh"
#include "EventFactory.hh"
#include "FileManager.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "TreeWriter.hh"
#include "Environment.hh"
#include "ObjectManager.hh"
#include "McSpectrumScaler.hh"

#define INFO_OUT_TAG "PreselectionTreeWriter"
#include "debugging.hh"

int main(int argc, char** argv) {

  // Command line option handling.
  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);

  config.SetProgramHelpText("PreselectionTreeWriter",
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
  
  unsigned int maxEntries =0;
  config.GetValue("OPTIONS","entries",maxEntries,"Number of events to process");

  // Load & parse cut configuration file.
  std::string cutConfigfile = "${MY_ANALYSIS}/Configuration/PreselectionTreeWriterforProtonCC.cfg";
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

  // 'AuxiliaryObjectManager' holds all auxiliary histograms / selectors created while processing the ACQt files.
  // NOTE: You should NOT write a TTree together with other histograms etc. in the ROOT file. You most likely
  // want to merge your histograms / selectors from batch jobs, but not the trees. That's why it's a good idea
  // in general to split up in two files: one for holding the tree, one for the rest.
  Utilities::ObjectManager auxiliaryObjectManager("AuxiliaryObjectManager", &config, resultDirectory, suffix);
  auxiliaryObjectManager.SetPrefix("ExampleAnalysis_Auxiliary");

  if (!config.PerformChecksAfterOptionParsing())
    return EXIT_FAIL_CONFIG;

  if (!fileManager.ReadFileList(inputList))
    return EXIT_FAIL_FILEMANAGER;

  // Construct tree manager which will manage the output file to hold the resulting tree.
  IO::TreeWriter treeWriter(new ExampleAnalysisTree, IO::TreeOptions::DontWriteInMemoryBranches);

  std::string treeFileName = Utilities::ObjectManager::MakeStandardRootFileName(resultDirectory, "ExampleAnalysis_Tree", suffix);
  treeWriter.Initialize(treeFileName);

  // Load cut selector(s).
  Cuts::Selector* mcPreselection = auxiliaryObjectManager.Add(selectionParser.GetSelector("McPreselection"));
  Cuts::Selector* badRuns = auxiliaryObjectManager.Add(selectionParser.GetSelector("BadRuns"));
  Cuts::Selector* rTI = auxiliaryObjectManager.Add(selectionParser.GetSelector("RTI"));

  Cuts::Selector* preselectionSelector = auxiliaryObjectManager.Add(selectionParser.GetSelector("Preselection"));
  Cuts::Selector* trackerCuts= auxiliaryObjectManager.Add(selectionParser.GetSelector("TrackerCuts"));
  Cuts::Selector* chargeOne = auxiliaryObjectManager.Add(selectionParser.GetSelector("ChargeOne"));
  Cuts::Selector* qualityCuts = auxiliaryObjectManager.Add(selectionParser.GetSelector("QualityCuts"));

  // Begin event loop.
  INFO_OUT_ON_MASTER << "Looping over " << fileManager.GetEntries() << " events..." << std::endl;

  static int sProductionSteps = Analysis::CreateSplineTrack | Analysis::CreateTrdTrack;
  unsigned int nEntries=0; 

  while (fileManager.GetNextEvent()) {
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
  preselectionSelector->PrintSummary();
  trackerCuts->PrintSummary();
  chargeOne->PrintSummary();
  qualityCuts->PrintSummary();
  mcPreselection->PrintSummary();
  badRuns->PrintSummary();
  rTI->PrintSummary();
  // Finish writing tree file.
  treeWriter.Finish();

  // Write auxiliary output file.
  auxiliaryObjectManager.WriteToFile();

  return EXIT_SUCCESS;
}
