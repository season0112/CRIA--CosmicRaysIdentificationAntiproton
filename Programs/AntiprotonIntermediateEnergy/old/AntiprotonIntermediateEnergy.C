#include "AntiprotonIntermediateEnergy.hh"

// ACsoft includes                                                                                                                     
#include "BinningFunctions.hh"
#include "BinningTools.hh"
#include "ConfigHandler.hh"
#include "Environment.hh"
#include "ObjectManager.hh"
#include "PredefinedBinnings.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "TreeWriter.hh"
#include "ValueHistograms.hh"
#include "CutAttachment.hh"
#include "AcceptanceManager.hh"
#include "EfficiencyHistograms.hh"
#include "McSpectrumScaler.hh"  
#include "Utilities.hh"


// ROOT includes                                                                                                                       
#include <TFile.h>
#include <TH1.h>
#include "TreeFormula.hh"
#include <TROOT.h>

#define INFO_OUT_TAG "AntiprotonIntermediateEnergy"
#include "debugging.hh"

AntiprotonIntermediateEnergy::AntiprotonIntermediateEnergy()
  : IO::TreeReader() {
}

AntiprotonIntermediateEnergy::~AntiprotonIntermediateEnergy() {
}

