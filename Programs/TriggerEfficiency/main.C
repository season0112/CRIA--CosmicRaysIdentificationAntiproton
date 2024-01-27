#include "ExampleAnalysisTree.hh"
#include "BayesUnfoldingWithCutoff.hh"
#include "AnalysisSettings.hh"

// Binning
#include "AntiprotonBinning.hh"

// ACsoft includes
#include "AnalysisEvent.hh"
#include "AcceptanceManager.hh"
#include "BinningDefinition.hh"
#include "BinningFunctions.hh"
#include "BinningTools.hh"
#include "ConfigHandler.hh"
#include "CutFactory.hh"
#include "CutAttachment.hh"
#include "EventFactory.hh"
#include "Environment.hh"
#include "EfficiencyHistograms.hh"
#include "FileManagerController.hh"
#include "FileManager.hh"
#include "GlobalOptions.hh"
#include "MPIEnvironment.hh"
#include "McSpectrumScaler.hh"
#include "ObjectManager.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "ValueHistograms.hh"

#include "TriggerEfficiencyHistograms.hh"

#include "TreeWriter.hh"
#include "TemplateFitter.hh"

#include <iostream>
#include <cassert>
#include <TH2.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCollection.h>
#include <TTree.h>
#include <TText.h>
#include "TemplateFitter2D.hh"
#include <vector>
#include <TCanvas.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include "IterativeUnfolding.hh"
#include "AcceptanceUnfolding.hh"
#include "Utilities.hh"
#include <TEfficiency.h>

#define INFO_OUT_TAG "TriggerEfficiency"
#include "debugging.hh"

int main(int argc, char* argv[]) {

    Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
    config.ReadCommandLine(argc, argv);

    config.SetProgramHelpText("TriggerEfficiency",
                              "Calculate TriggerEfficiency");
    config.AddHelpExample("TriggerEfficiency", "--binningversion 450version");
 
    std::string binningversion = "";
    config.GetValue("OPTIONS", "binningversion", binningversion,
                    "The binningversion is");

    std::string pattern = "";
    config.GetValue("OPTIONS", "pattern", pattern,
                    "The tracker pattern is");

    std::string RigidityRange = "";
    config.GetValue("OPTIONS", "RigidityRange", RigidityRange,
                    "The RigidityRange is");
 
////////////////////////////////////////////////////////////////////////////////////////////

    std::string BinningversionName;
    if (binningversion == "525version") {
      BinningversionName = "525version";
    }
    else if (binningversion == "450version") {
      BinningversionName = "450version";
    }
    std::string PatternName;
    if (pattern == "all"){
      PatternName = "";
    }
    else if (pattern == "-1"){
      PatternName = "_Pattern-1";
    }
    else if (pattern == "0"){
      PatternName = "_Pattern0";
    }
    else if (pattern == "1"){
      PatternName = "_Pattern1";
    }
    else if (pattern == "2"){
      PatternName = "_Pattern2";
    }
    else if (pattern == "3"){
      PatternName = "_Pattern3";
    }
    else if (pattern == "4"){
      PatternName = "_Pattern4";
    }
    else if (pattern == "5"){
      PatternName = "_Pattern5";
    }

    std::string MCUsed = "B1042_pr.pl1.1800_7.6_all_";
    //std::string MCUsed = "B1042_antipr.pl1.1800_7.6_all_";

    TFile *f1 = new TFile();
    if (RigidityRange == "HIGHRANGE"){
      chdir( (std::string(getenv("HPCHIGHENERGYDATADIR")) ).c_str());
      f1 = new TFile ( ( MCUsed + BinningversionName + PatternName + std::string("_Auxiliary.root")).c_str() );
    }

    else if (RigidityRange == "INTERMEDIATERANGE"){
      chdir( (std::string(getenv("HPCINTERMEDIATEDIR")) ).c_str());
      f1 = new TFile ( (MCUsed + std::string("Auxiliary_negative.root")).c_str()  );
    }

    else if (RigidityRange == "LOWRANGE"){
      chdir( (std::string(getenv("HPCLOWENERGYDIR")) ).c_str());
      f1 = new TFile ( (MCUsed + std::string("Auxiliary.root")).c_str()  );

    }


    Cuts::Selector *McPreselection = (Cuts::Selector*)f1->Get("McPreselection");
    Cuts::Selector *Preselection   = (Cuts::Selector*)f1->Get("Preselection");
    Cuts::Selector *QualityCuts    = (Cuts::Selector*)f1->Get("QualityCuts");

    TH1F *PhysicsTriggerHisto           = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->PhysicsTriggerHistogram();
    TH1F *PhysicsTriggerHistoTof        = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->PhysicsTriggerHistogramTof();
    TH1F *PhysicsTriggerHistoEcal       = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->PhysicsTriggerHistogramEcal();
    TH1F *PhysicsTriggerHistoTofAndEcal = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->PhysicsTriggerHistogramTofAndEcal();
    TH1F *UnbiasedTofTriggerHisto        = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->UnbiasedTofTriggerHistogram();
    TH1F *UnbiasedEcalTriggerHisto       = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->UnbiasedEcalTriggerHistogram();
    TH1F *UnbiasedTofAndEcalTriggerHisto = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->UnbiasedTofAndEcalTriggerHistogram();

    TEfficiency *effMC_Preselection   = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->ProduceTriggerEfficiencyMc();
    TEfficiency *effData_Preselection = static_cast<const Cuts::TriggerEfficiencyHistograms*>(Preselection->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->ProduceTriggerEfficiencyData();

    TEfficiency *effMC_QualityCuts   = static_cast<const Cuts::TriggerEfficiencyHistograms*>(QualityCuts->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->ProduceTriggerEfficiencyMc();
    TEfficiency *effData_QualityCuts = static_cast<const Cuts::TriggerEfficiencyHistograms*>(QualityCuts->FindAttachment("Cuts::TriggerEfficiencyHistograms"))->ProduceTriggerEfficiencyData();

    TH1F TriggerEff              = *PhysicsTriggerHisto / (*PhysicsTriggerHisto + 100 * *UnbiasedTofTriggerHisto + 1000 * *UnbiasedEcalTriggerHisto + 90.99 * *UnbiasedTofAndEcalTriggerHisto );
    TH1F TriggerEff_noprescaling = *PhysicsTriggerHisto / (*PhysicsTriggerHisto + *UnbiasedTofTriggerHisto + *UnbiasedEcalTriggerHisto + *UnbiasedTofAndEcalTriggerHisto );


    if (RigidityRange == "HIGHRANGE"){
        TFile *f2 = new TFile( (std::string("TriggerEff_") + MCUsed + std::string(BinningversionName) + PatternName +std::string(".root")).c_str(), "RECREATE");
        PhysicsTriggerHisto ->Write("PhysicsTriggerHisto");
        effMC_Preselection  ->Write("effMC_Preselection");
        effData_Preselection->Write("effData_Preselection");
        effMC_QualityCuts   ->Write("effMC_QualityCuts");
        effData_QualityCuts ->Write("effData_QualityCuts");
        TriggerEff.Write("TriggerEff");
        TriggerEff_noprescaling.Write("TriggerEff_noprescaling");
        f2->Close();
    }
    else if (RigidityRange == "INTERMEDIATERANGE" or RigidityRange == "LOWRANGE"){
        TFile *f2 = new TFile( (std::string("TriggerEff_") + MCUsed + std::string(BinningversionName) + PatternName +std::string(".root")).c_str(), "RECREATE");
        PhysicsTriggerHisto ->Write("PhysicsTriggerHisto");
        effMC_Preselection  ->Write("effMC_Preselection");
        effData_Preselection->Write("effData_Preselection");
        effMC_QualityCuts   ->Write("effMC_QualityCuts");
        effData_QualityCuts ->Write("effData_QualityCuts");
        TriggerEff.Write("TriggerEff");
        TriggerEff_noprescaling.Write("TriggerEff_noprescaling");
        f2->Close();
    }
}








