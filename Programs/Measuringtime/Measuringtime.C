// ACsoft includes
#include "AnalysisEvent.hh"
#include "AcceptanceManager.hh"
#include "BinningDefinition.hh"
#include "BinningFunctions.hh"
#include "BinningTools.hh"
#include "ConfigHandler.hh"
#include "CutFactory.hh"
#include "CutAttachment.hh"
#include "Environment.hh"
#include "EfficiencyHistograms.hh"
#include "EventFactory.hh"
#include "FileManagerController.hh"
#include "FileManager.hh"
#include "GlobalOptions.hh"
#include "MPIEnvironment.hh"
#include "McSpectrumScaler.hh"
#include "MeasuringTime.hh"
#include "ObjectManager.hh"
#include "PredefinedBinnings.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "ValueHistograms.hh"
// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include "TreeFormula.hh"
#include "Utilities.hh"
#include <TROOT.h>
#include "TreeWriter.hh"
#include <TApplication.h>
#include <TProof.h>

#define INFO_OUT_TAG "rti_reader_example"
#include "debugging.hh"
#include "AntiprotonBinning.hh"

int main(int argc, char** argv) {

  //// Command line option handling.
  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);
  config.SetProgramHelpText("Measuringtime",
                            "Illustrates the calculation of total livetime for a given time period based on MeasuringTime class.");
  config.AddHelpExample("Case #1): Loop directly over RTI files.", "--startTime 1305853510 --endTime 1432628592");
  config.AddHelpExample("Case #2): Loop over RTI tree written by rti_writer_example.", "--RTI/TreePattern RTITree*.root --startTime 1305853510 --endTime 1432628592");

  int startTime = -1; 
  int endTime   = -1;
  std::string cutoffmode   = "";
  std::string degree       = "";
  std::string safetyfactor = "";
  config.GetValue("OPTIONS", "startTime"   , startTime,
                  "Start time in seconds as Unix time stamp.");
  config.GetValue("OPTIONS", "endTime"     , endTime,
                  "End time in seconds as Unix time stamp.");
  config.GetValue("OPTIONS", "cutoffmode"  , cutoffmode,
                  "which cutoffmode you use, GEOMETRIC or IGRF.");
  config.GetValue("OPTIONS", "degree"      , degree,
                  "Which degree you use.");
  config.GetValue("OPTIONS", "safetyfactor", safetyfactor,
                  "Which safetyfactor you use.");

  if (cutoffmode == "GEOMETRIC"){
      chdir( (std::string("/hpcwork/jara0052/sichen/Measuringtime/05_2021_GEOMETRIC") + degree + std::string("_") + safetyfactor + std::string("/")).c_str());
  }
  else if (cutoffmode == "IGRF"){
     chdir( (std::string("/hpcwork/jara0052/sichen/Measuringtime/05_2021_IGRF") + degree + std::string("_") + safetyfactor + std::string("/")).c_str());
  }

  Utilities::ObjectManager objectManager(&config, "", "");
  if (!config.PerformChecksAfterOptionParsing())
    return EXIT_FAIL_CONFIG;

  if (startTime == -1 || endTime == -1)
    FATAL_OUT << "You have to specify both --startTime X / --endTime Y." << std::endl;


  //// Preperations to calculate the measuring time:

  // 1. Select a geomagnetic cut-off cut, the field-of-view (here: 40 degress for both positive/negative particle hypothesis) and a safety factor
  Cuts::Cut* cutOffCut;
  if (cutoffmode == "GEOMETRIC"){
      cutOffCut = Cuts::CreateCut("RigidityAboveGeomagneticCutoff", (degree + std::string("PN")).c_str(), std::stoi(safetyfactor));
  }
  else if (cutoffmode == "IGRF"){
      cutOffCut = Cuts::CreateCut("RigidityAboveIGRFCutoff"       , (degree + std::string("PN")).c_str(), std::stoi(safetyfactor));
  }

  // 2. Open the cut config file containing at least "BadRuns" / "RTI" selectors, optionally also "TrdCalibration".
  std::string cutConfigfile = "${MY_ANALYSIS}/Configuration/RTIandBadRuns.cfg";
  Environment::ExpandEnvironmentVariables(cutConfigfile);
  config.Read(cutConfigfile);

  #ifdef ENABLE_MPI
  if (IO::MPIEnvironment::IsMPIEnabled())
    MPI_Init(&argc, &argv);
  #endif

  // NOTE: These two lines are only needed in a standalone program that does NOT loop over ACQt files.
  // If you've instantiated an IO::FileManager before, this is automatically deduced from the file list.
  IO::FileManagerController::Self()->SetRunType(AC::ISSRun);
  IO::FileManagerController::Self()->SetFirstAndLastEventTimes(startTime, endTime);

  const std::vector<std::string>& inputTreeFiles = GlobalOptions::Self()->RTITreeFileName();
  if (inputTreeFiles.empty())
    objectManager.SetPrefix("MeasuringTime_from_RTI_Files");
  else
    objectManager.SetPrefix("MeasuringTime_from_RTI_Tree");

  RTI::MeasuringTime measuringTime(config, objectManager, cutOffCut, AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili());
  measuringTime.ComputeMeasuringTime();

  objectManager.WriteToFile();

  #ifdef ENABLE_MPI
  if (IO::MPIEnvironment::IsMPIEnabled())
    MPI_Finalize();
  #endif

  return EXIT_SUCCESS;
}
