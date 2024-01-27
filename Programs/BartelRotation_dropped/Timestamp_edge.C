#include "AntiprotonIntermediateEnergyTree.hh"
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

#include "PredefinedBinnings.hh"
#include "TimeUtilities.hh"
#include "AMSUnixTimes.hh"

#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>

#define INFO_OUT_TAG "Timestamp_edge"
#include "debugging.hh"

#define BeginOfIssTime 1305800000.0
#define EndOfPass7Time 1527491432.0
#define EndOfPass7extTime 1546405413.0


int main(int argc, char* argv[]) {
(void) argc;
(void) argv;


for (int i=2426;i<2521;i++){
    std::cout << Utilities::Time::StartTimeOfBartelsRotation(i) << std::endl;
}

  return EXIT_SUCCESS;
}








