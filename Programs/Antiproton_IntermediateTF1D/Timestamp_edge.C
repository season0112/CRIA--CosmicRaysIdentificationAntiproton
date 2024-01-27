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

#define BeginOfIssTime 1305800000.0 // 05/19/2011 @ 10:13am (UTC)
#define EndOfPass7Time 1527491432.0 // 05/28/2018 @ 7:10am (UTC)
#define EndOfPass7extTime 1546405413.0 // 01/02/2019 @ 5:03am (UTC)
#define sixmonths 15552000 // 6*30*24*60*60
#define oneBartelRotation 2332800 // 

int main(int argc, char* argv[]) {
(void) argc;
(void) argv;

std::cout << "Bartel Rotations:"<< std::endl;
for (int i=2426;i<2580;i++){
    std::cout << Utilities::Time::StartTimeOfBartelsRotation(i) << std::endl;
}

std::cout << "6 Months:"<< std::endl;
for (int j=0;j<25;j++){
    std::cout << 1305800000+j*sixmonths << std::endl;
}
  return EXIT_SUCCESS;
}








