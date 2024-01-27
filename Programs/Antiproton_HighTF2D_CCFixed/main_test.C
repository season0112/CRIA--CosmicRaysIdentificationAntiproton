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
#include <iostream>
#include <string>
#include <cassert>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <vector>
#include <TCanvas.h>
#include <sstream>
#include <string>
#include <unistd.h>

#include "TemplateFitterforProtonCC.hh"
#include "TemplateFitter2DforProtonCC.hh"
#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>
#include <stdlib.h>

#include "TemplateFitter2DforProtonCC_Unbinned.hh"

#include "AntiprotonAnalysisTools.hh"

using namespace std;

#define INFO_OUT_TAG "Antiproton_HighTF2D_Unbinned"
#include "debugging.hh"

int main(int argc, char* argv[]) {

  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);

  config.SetProgramHelpText("Antiproton_HighTF2D_Unbinned",
                            "Template fit for antiproton signal determination in High Rigidity Range.");

  config.AddHelpExample("Antiproton_HighTF2D_Unbinned", "--binningversion 450version --Rigidityofdataset data_negative");

  std::string binningversion = "";
  config.GetValue("OPTIONS", "binningversion", binningversion,
                  "The binningversion is");

  std::string issversion = "";
  config.GetValue("OPTIONS", "issversion", issversion,
                  "The issversion is");

  std::string sign_of_data = "";
  config.GetValue("OPTIONS", "Rigidityofdataset", sign_of_data,
                  "The sign_of_data is:(data_negative or data_positive)");

  if (binningversion == "") {
    WARN_OUT << "No binningversion is given! Please give a binningversion." << std::endl;
    return EXIT_FAIL_CONFIG;
  }

  if (sign_of_data == "") {
    WARN_OUT << "No Rigidityofdataset is given! Please give a Rigidityofdataset." << std::endl;
    return EXIT_FAIL_CONFIG;
  }


/* old
bool Perform2DTemplateFit(SampleIdentifier identifier, TGraph* chargeConfusionGraph, TH2D* dataHistogram,
                            TH2D* electronTemplate, TH2D* positronTemplate,
                            TH2D* ccElectronTemplate, TH2D* ccPositronTemplate,
                            TH2D* protonTemplate, TH2D* ccProtonTemplate,
                            TemplateFit2DResults& fitResults, TemplateFit2DResults* startValues) 
TH2D *dataHistogram
TH2D *electronTemplate
TH2D *positronTemplate
TH2D *ccElectronTemplate
TH2D *ccPositronTemplate
TH2D *protonTemplate
TH2D *ccProtonTemplate

TemplateFit2DFitFunction fitFunction(dataHistogram,
                                     electronTemplate, positronTemplate,
                                     ccElectronTemplate, ccPositronTemplate,
                                     protonTemplate, ccProtonTemplate);

TemplateFit2DParameters fitParameters = GenerateParameters(identifier, chargeConfusionGraph, startValues);
*/


TH2D *template_correct
TH2D *template_confused
TH2D *template_electron_Data
TH2D *data

TemplateFit2DFitFunction fitFunction(data, template_correct, template_confused, template_electron_Data);












  return EXIT_SUCCESS;
}



