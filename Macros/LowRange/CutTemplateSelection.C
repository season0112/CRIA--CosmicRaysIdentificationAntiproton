// Useage: root -L /home/bo791269/Software/AntiprotonAnalysis/Libraries/AntiprotonBinning/AntiprotonBinning.C -x CutTemplateSelection.C
#include <sstream>
#include "CutTemplateSelection.hh"

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
#include <TCut.h>
#include <TFile.h>
#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <TCanvas.h>
#include <TF1.h>
#include <sstream>
#include <unistd.h>
#include "Utilities.hh"
#include <TChain.h>
#include <TTree.h>

#include "TBox.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"


void CutTemplateSelection(){

    /*
    //// Cuts and selections Definition
    // Tracker
    std::string trackerpattern = "0";
    TCut PatternCut = (string("Pattern==") + trackerpattern).c_str();
    TCut ExtrapolatedPhotoElectronsFromTracker ="NPhotoElectrons-ExtrapolatedRichExpectedPhotoElectronsProton<20 && NPhotoElectrons>-999 && ExtrapolatedRichTileIndex == TileIndex || NPhotoElectrons==-999 || ExtrapolatedRichTileIndex != TileIndex";
    TCut ProtonCCMVABDCut = "ProtonCCMVABDT > 0.9";
    //Trigger
    TCut TriggerPhysics = "TriggerFlags != 1 && TriggerFlags != 64 && TriggerFlags != 65";
    //TRD
    TCut TrdLikelihoodCut = "TrdLogLikelihoodRatioElectronProtonTracker > 0.7 || TrdLogLikelihoodRatioElectronProtonTracker==-1.5";    // Antiproton:1to1.2, Electron:0.5
    TCut TrdLikelihoodHeProtonCut = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.1";
    //ECAL
    TCut EcalBDT_EnergyDCut = "EcalBDT_EnergyD < -0.9";
    //RICH
    TCut RichBetaCut = "RichBeta==0 || RichIsNaF==1";
    TCut RichNaFCut = "RichIsNaF==1";
    TCut RichAglCut = "RichIsNaF==0";
    TCut RichPhotoElectron = "NPhotoElectrons-NExpectedPhotoElectrons<10  || NPhotoElectrons==-999";
    TCut RichChargeCut = "RichCharge<2 || RichCharge==0";
    //TOF
    //TCut BetaConverted = "BetaConverted<0.9";
    TCut TofMassonecharge = "TofMassonecharge>0.5";
    //TCut TofMassonecharge = "TofMassonecharge>0";
    //TCut TOFBETALikelihood = "TOFBETALikelihood<0.7";
    //PionMVA
    TCut PionMVA = "PionBDT>-999 && PionBDT>-0";
    TCut PionMVA_light = "PionBDT>-999.0 && PionBDT<0.7";

    
    //// Templates Cuts Definition
    TCut AntiprotonCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
    //TCut ElectronCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut ;
    //TCut ElectronCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut ;
    TCut ElectronCut;
    TCut PionCutMC = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && RichBetaCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
    TCut PionCutMCNew = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>0.72";
    TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
    TCut PositiveCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
    TCut PionTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>0.72";
    TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-1)<0.002 && EcalBDT_EnergyD>-999 && EcalBDT_EnergyD>0.9"; 
    */


    //// Load data
    vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());

    //i=0,i<3: 0.8-1.3GV;  i=3,i<6:1.33-1.92GV;  i=4,i<7,1.51-2.15GV;  i=4,i<5,1.51-1.71GV;  i=1,i<16: 1.00-5.37GV;   i=1,i<30: 1.00-18.0GV;
    int indexmin = 1;
    int indexmax = 30;

    for (auto i=indexmin; i<indexmax; i++){ 

        //// Define Rigidity bins 
        std::ostringstream lowedge_stream;
        lowedge_stream << binning.at(i);
        std::string lowedge = lowedge_stream.str();
        std::ostringstream highedge_stream;
        highedge_stream << binning.at(i+1);
        std::string highedge = highedge_stream.str();

        if ( lowedge == "1" || lowedge == "11" || lowedge == "12" || lowedge == "13" || lowedge == "18" ){
            lowedge = to_string_with_precision(binning[i], 1); }
        if ( highedge == "1" || highedge == "11" || highedge == "12" || highedge == "13" || highedge == "18" ){
            highedge = to_string_with_precision(binning[i+1], 1);}

        cout<< "Rigidity bin is " << lowedge << "_" << highedge <<  " GV Now." << endl; 


        //// Load TChain File
        TChain *IssNegativeFile  = new TChain("AntiprotonLowEnergyTree");
        TChain *PionMcFile       = new TChain("AntiprotonLowEnergyTree");
        TChain *ElectronMcFile   = new TChain("AntiprotonLowEnergyTree");
        TChain *AntiprotonMcFile = new TChain("AntiprotonLowEnergyTree");

        //IssNegativeFile->AddFile( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1130_pass7_7.8_all_Tree_negative_") + lowedge + string("_") + highedge + string(".root")) .c_str() );    
        IssNegativeFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1130_pass7_7.8_all_Tree_negative_May2015_") + lowedge + string("_") + highedge + string(".root") ) .c_str() );
        PionMcFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + lowedge + string("_") + highedge + string(".root")) .c_str() );
        ElectronMcFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1091_el.pl1.0_25200_7.6_all_Tree_negative_") + lowedge + string("_") +  highedge + string(".root")) .c_str() );
        AntiprotonMcFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1042_antipr.pl1.1800_7.6_all_Tree_") + lowedge + string("_") + highedge + string(".root")) .c_str() );
    

        //// Build Histogram from TChain
        // Build 1D TrdSegmentsXZNumber Plot 
        IssNegativeFile->Draw("TrdSegmentsXZNumber>>IssNegativeData_TrdSegmentsXZNumber(100, -0.2, 10)","");
        IssNegativeFile->Draw("TrdSegmentsXZNumber>>IssNegativeData_TrdSegmentsXZNumber_withcut(100, -0.2, 10)","TrdSegmentsXZNumber==1");
        TH1F *IssNegativeData_TrdSegmentsXZNumber = (TH1F*)gDirectory->Get("IssNegativeData_TrdSegmentsXZNumber");
        TH1F *IssNegativeData_TrdSegmentsXZNumber_withcut = (TH1F*)gDirectory->Get("IssNegativeData_TrdSegmentsXZNumber_withcut");

        PionMcFile->Draw("TrdSegmentsXZNumber>>PionMC_TrdSegmentsXZNumber(100, -0.2, 10)","");
        PionMcFile->Draw("TrdSegmentsXZNumber>>PionMC_TrdSegmentsXZNumber_withcut(100, -0.2, 10)","TrdSegmentsXZNumber==1");
        TH1F *PionMC_TrdSegmentsXZNumber = (TH1F*)gDirectory->Get("PionMC_TrdSegmentsXZNumber");
        TH1F *PionMC_TrdSegmentsXZNumber_withcut = (TH1F*)gDirectory->Get("PionMC_TrdSegmentsXZNumber_withcut");
        cout<< "Pion TrdSegmentsXZNumber:" << PionMC_TrdSegmentsXZNumber_withcut->GetEntries() / PionMC_TrdSegmentsXZNumber->GetEntries() <<endl;

        ElectronMcFile->Draw("TrdSegmentsXZNumber>>ElectronMC_TrdSegmentsXZNumber(100, -0.2, 10)","");
        ElectronMcFile->Draw("TrdSegmentsXZNumber>>ElectronMC_TrdSegmentsXZNumber_withcut(100, -0.2, 10)","TrdSegmentsXZNumber==1");
        TH1F *ElectronMC_TrdSegmentsXZNumber = (TH1F*)gDirectory->Get("ElectronMC_TrdSegmentsXZNumber");
        TH1F *ElectronMC_TrdSegmentsXZNumber_withcut = (TH1F*)gDirectory->Get("ElectronMC_TrdSegmentsXZNumber_withcut");
        cout<< "Electron TrdSegmentsXZNumber:" << ElectronMC_TrdSegmentsXZNumber_withcut->GetEntries() / ElectronMC_TrdSegmentsXZNumber->GetEntries() <<endl;

        AntiprotonMcFile->Draw("TrdSegmentsXZNumber>>AntiprotonMC_TrdSegmentsXZNumber(100, -0.2, 10)","");
        AntiprotonMcFile->Draw("TrdSegmentsXZNumber>>AntiprotonMC_TrdSegmentsXZNumber_withcut(100, -0.2, 10)","TrdSegmentsXZNumber==1");
        TH1F *AntiprotonMC_TrdSegmentsXZNumber = (TH1F*)gDirectory->Get("AntiprotonMC_TrdSegmentsXZNumber");
        TH1F *AntiprotonMC_TrdSegmentsXZNumber_withcut = (TH1F*)gDirectory->Get("AntiprotonMC_TrdSegmentsXZNumber_withcut");    
        cout<< "Antiproton TrdSegmentsXZNumber:" << AntiprotonMC_TrdSegmentsXZNumber_withcut->GetEntries() / AntiprotonMC_TrdSegmentsXZNumber->GetEntries() <<endl;

        // Build 1D TrdSegmentsYZNumber Plot
        IssNegativeFile->Draw("TrdSegmentsYZNumber>>IssNegativeData_TrdSegmentsYZNumber(100, -0.2, 10)","");
        IssNegativeFile->Draw("TrdSegmentsYZNumber>>IssNegativeData_TrdSegmentsYZNumber_withcut(100, -0.2, 10)","TrdSegmentsYZNumber==1");
        TH1F *IssNegativeData_TrdSegmentsYZNumber = (TH1F*)gDirectory->Get("IssNegativeData_TrdSegmentsYZNumber");
        TH1F *IssNegativeData_TrdSegmentsYZNumber_withcut = (TH1F*)gDirectory->Get("IssNegativeData_TrdSegmentsYZNumber_withcut");

        PionMcFile->Draw("TrdSegmentsYZNumber>>PionMC_TrdSegmentsYZNumber(100, -0.2, 10)","");
        PionMcFile->Draw("TrdSegmentsYZNumber>>PionMC_TrdSegmentsYZNumber_withcut(100, -0.2, 10)","TrdSegmentsYZNumber==1");
        TH1F *PionMC_TrdSegmentsYZNumber = (TH1F*)gDirectory->Get("PionMC_TrdSegmentsYZNumber");
        TH1F *PionMC_TrdSegmentsYZNumber_withcut = (TH1F*)gDirectory->Get("PionMC_TrdSegmentsYZNumber_withcut");
        cout<< "Pion TrdSegmentsYZNumber:"<< PionMC_TrdSegmentsYZNumber_withcut->GetEntries() / PionMC_TrdSegmentsYZNumber->GetEntries() <<endl;

        ElectronMcFile->Draw("TrdSegmentsYZNumber>>ElectronMC_TrdSegmentsYZNumber(100, -0.2, 10)","");
        ElectronMcFile->Draw("TrdSegmentsYZNumber>>ElectronMC_TrdSegmentsYZNumber_withcut(100, -0.2, 10)","TrdSegmentsYZNumber==1");
        TH1F *ElectronMC_TrdSegmentsYZNumber = (TH1F*)gDirectory->Get("ElectronMC_TrdSegmentsYZNumber");
        TH1F *ElectronMC_TrdSegmentsYZNumber_withcut = (TH1F*)gDirectory->Get("ElectronMC_TrdSegmentsYZNumber_withcut");
        cout<< "Electron TrdSegmentsYZNumber:"<< ElectronMC_TrdSegmentsYZNumber_withcut->GetEntries() / ElectronMC_TrdSegmentsYZNumber->GetEntries() <<endl;

        AntiprotonMcFile->Draw("TrdSegmentsYZNumber>>AntiprotonMC_TrdSegmentsYZNumber(100, -0.2, 10)","");
        AntiprotonMcFile->Draw("TrdSegmentsYZNumber>>AntiprotonMC_TrdSegmentsYZNumber_withcut(100, -0.2, 10)","TrdSegmentsYZNumber==1");
        TH1F *AntiprotonMC_TrdSegmentsYZNumber = (TH1F*)gDirectory->Get("AntiprotonMC_TrdSegmentsYZNumber");
        TH1F *AntiprotonMC_TrdSegmentsYZNumber_withcut = (TH1F*)gDirectory->Get("AntiprotonMC_TrdSegmentsYZNumber_withcut");
        cout<< "Antiproton TrdSegmentsYZNumber:"<< AntiprotonMC_TrdSegmentsYZNumber_withcut->GetEntries() / AntiprotonMC_TrdSegmentsYZNumber->GetEntries()<<endl;

        // Build 1D TrdNumberOfHits Plot
        IssNegativeFile->Draw("TrdNumberOfHits>>IssNegativeData_TrdNumberOfHits(100, 0, 200)","");
        IssNegativeFile->Draw("TrdNumberOfHits>>IssNegativeData_TrdNumberOfHits_withcut(100, 0, 200)","TrdNumberOfHits<35");
        TH1F *IssNegativeData_TrdNumberOfHits = (TH1F*)gDirectory->Get("IssNegativeData_TrdNumberOfHits");
        TH1F *IssNegativeData_TrdNumberOfHits_withcut = (TH1F*)gDirectory->Get("IssNegativeData_TrdNumberOfHits_withcut");

        PionMcFile->Draw("TrdNumberOfHits>>PionMC_TrdNumberOfHits(100, 0, 200)","");
        PionMcFile->Draw("TrdNumberOfHits>>PionMC_TrdNumberOfHits_withcut(100, 0, 200)","TrdNumberOfHits<35");
        TH1F *PionMC_TrdNumberOfHits = (TH1F*)gDirectory->Get("PionMC_TrdNumberOfHits");
        TH1F *PionMC_TrdNumberOfHits_withcut = (TH1F*)gDirectory->Get("PionMC_TrdNumberOfHits_withcut");
        cout<< "Pion TrdNumberOfHits:"<< PionMC_TrdNumberOfHits_withcut->GetEntries() / PionMC_TrdNumberOfHits->GetEntries() <<endl;

        ElectronMcFile->Draw("TrdNumberOfHits>>ElectronMC_TrdNumberOfHits(100, 0, 200)","");
        ElectronMcFile->Draw("TrdNumberOfHits>>ElectronMC_TrdNumberOfHits_withcut(100, 0, 200)","TrdNumberOfHits<35");
        TH1F *ElectronMC_TrdNumberOfHits = (TH1F*)gDirectory->Get("ElectronMC_TrdNumberOfHits");
        TH1F *ElectronMC_TrdNumberOfHits_withcut = (TH1F*)gDirectory->Get("ElectronMC_TrdNumberOfHits_withcut");
        cout<< "Electron TrdNumberOfHits:"<< ElectronMC_TrdNumberOfHits_withcut->GetEntries() / ElectronMC_TrdNumberOfHits->GetEntries() << endl;

        AntiprotonMcFile->Draw("TrdNumberOfHits>>AntiprotonMC_TrdNumberOfHits(100, 0, 200)","");
        AntiprotonMcFile->Draw("TrdNumberOfHits>>AntiprotonMC_TrdNumberOfHits_withcut(100, 0, 200)","TrdNumberOfHits<35");
        TH1F *AntiprotonMC_TrdNumberOfHits = (TH1F*)gDirectory->Get("AntiprotonMC_TrdNumberOfHits");
        TH1F *AntiprotonMC_TrdNumberOfHits_withcut = (TH1F*)gDirectory->Get("AntiprotonMC_TrdNumberOfHits_withcut");
        cout<< "Antiproton TrdNumberOfHits:"<< AntiprotonMC_TrdNumberOfHits_withcut->GetEntries() / AntiprotonMC_TrdNumberOfHits->GetEntries()<<endl;

       // Build 1D TRDVTracksSize Plot
        IssNegativeFile->Draw("TRDVTracksSize>>IssNegativeData_TRDVTracksSize(100, -1, 10)","");
        IssNegativeFile->Draw("TRDVTracksSize>>IssNegativeData_TRDVTracksSize_withcut(100, -1, 10)","TRDVTracksSize==1");
        TH1F *IssNegativeData_TRDVTracksSize = (TH1F*)gDirectory->Get("IssNegativeData_TRDVTracksSize");
        TH1F *IssNegativeData_TRDVTracksSize_withcut = (TH1F*)gDirectory->Get("IssNegativeData_TRDVTracksSize_withcut");

        PionMcFile->Draw("TRDVTracksSize>>PionMC_TRDVTracksSize(100, -1, 10)","");
        PionMcFile->Draw("TRDVTracksSize>>PionMC_TRDVTracksSize_withcut(100, -1, 10)","TRDVTracksSize==1");
        TH1F *PionMC_TRDVTracksSize = (TH1F*)gDirectory->Get("PionMC_TRDVTracksSize");
        TH1F *PionMC_TRDVTracksSize_withcut = (TH1F*)gDirectory->Get("PionMC_TRDVTracksSize_withcut");
        cout<< "Pion TRDVTracksSize:"<< PionMC_TRDVTracksSize_withcut->GetEntries() / PionMC_TRDVTracksSize->GetEntries() <<endl;

        ElectronMcFile->Draw("TRDVTracksSize>>ElectronMC_TRDVTracksSize(100, -1, 10)","");
        ElectronMcFile->Draw("TRDVTracksSize>>ElectronMC_TRDVTracksSize_withcut(100, -1, 10)","TRDVTracksSize==1");
        TH1F *ElectronMC_TRDVTracksSize = (TH1F*)gDirectory->Get("ElectronMC_TRDVTracksSize");
        TH1F *ElectronMC_TRDVTracksSize_withcut = (TH1F*)gDirectory->Get("ElectronMC_TRDVTracksSize_withcut");
        cout<< "Electron TRDVTracksSize:"<< ElectronMC_TRDVTracksSize_withcut->GetEntries() / ElectronMC_TRDVTracksSize->GetEntries() << endl;

        AntiprotonMcFile->Draw("TRDVTracksSize>>AntiprotonMC_TRDVTracksSize(100, -1, 10)","");
        AntiprotonMcFile->Draw("TRDVTracksSize>>AntiprotonMC_TRDVTracksSize_withcut(100, -1, 10)","TRDVTracksSize==1");
        TH1F *AntiprotonMC_TRDVTracksSize = (TH1F*)gDirectory->Get("AntiprotonMC_TRDVTracksSize");
        TH1F *AntiprotonMC_TRDVTracksSize_withcut = (TH1F*)gDirectory->Get("AntiprotonMC_TRDVTracksSize_withcut");
        cout<< "Antiproton TRDVTracksSize:"<< AntiprotonMC_TRDVTracksSize_withcut->GetEntries() / AntiprotonMC_TRDVTracksSize->GetEntries()<<endl;

        // Build 1D ACCHits Plot
        IssNegativeFile->Draw("ACCHits>>IssNegativeData_ACCHits(100, -1, 10)","");
        IssNegativeFile->Draw("ACCHits>>IssNegativeData_ACCHits_withcut(100, -1, 10)","ACCHits==0");
        TH1F *IssNegativeData_ACCHits = (TH1F*)gDirectory->Get("IssNegativeData_ACCHits");
        TH1F *IssNegativeData_ACCHits_withcut = (TH1F*)gDirectory->Get("IssNegativeData_ACCHits_withcut");

        PionMcFile->Draw("ACCHits>>PionMC_ACCHits(100, -1, 10)","");
        PionMcFile->Draw("ACCHits>>PionMC_ACCHits_withcut(100, -1, 10)","ACCHits==0");
        TH1F *PionMC_ACCHits = (TH1F*)gDirectory->Get("PionMC_ACCHits");
        TH1F *PionMC_ACCHits_withcut = (TH1F*)gDirectory->Get("PionMC_ACCHits_withcut");
        cout<< "Pion ACCHits:"<< PionMC_ACCHits_withcut->GetEntries() / PionMC_ACCHits->GetEntries() <<endl;

        ElectronMcFile->Draw("ACCHits>>ElectronMC_ACCHits(100, -1, 10)","");
        ElectronMcFile->Draw("ACCHits>>ElectronMC_ACCHits_withcut(100, -1, 10)","ACCHits==0");
        TH1F *ElectronMC_ACCHits = (TH1F*)gDirectory->Get("ElectronMC_ACCHits");
        TH1F *ElectronMC_ACCHits_withcut = (TH1F*)gDirectory->Get("ElectronMC_ACCHits_withcut");
        cout<< "Electron ACCHits:"<< ElectronMC_ACCHits_withcut->GetEntries() / ElectronMC_ACCHits->GetEntries() << endl;

        AntiprotonMcFile->Draw("ACCHits>>AntiprotonMC_ACCHits(100, -1, 10)","");
        AntiprotonMcFile->Draw("ACCHits>>AntiprotonMC_ACCHits_withcut(100, -1, 10)","ACCHits==0");
        TH1F *AntiprotonMC_ACCHits = (TH1F*)gDirectory->Get("AntiprotonMC_ACCHits");
        TH1F *AntiprotonMC_ACCHits_withcut = (TH1F*)gDirectory->Get("AntiprotonMC_ACCHits_withcut");
        cout<< "Antiproton ACCHits:"<< AntiprotonMC_ACCHits_withcut->GetEntries() / AntiprotonMC_ACCHits->GetEntries()<<endl;
        
        // Build 1D EcalBDT Plot
        IssNegativeFile->Draw("EcalBDT_EnergyD>>IssNegativeData_EcalBDT_EnergyD(100, -1, 1)","");
        IssNegativeFile->Draw("EcalBDT_EnergyD>>IssNegativeData_EcalBDT_EnergyD_withcut(100, -1, 1)","EcalBDT_EnergyD>-1");
        TH1F *IssNegativeData_EcalBDT_EnergyD = (TH1F*)gDirectory->Get("IssNegativeData_EcalBDT_EnergyD");
        TH1F *IssNegativeData_EcalBDT_EnergyD_withcut = (TH1F*)gDirectory->Get("IssNegativeData_EcalBDT_EnergyD_withcut");

        PionMcFile->Draw("EcalBDT_EnergyD>>PionMC_EcalBDT_EnergyD(100, -1, 1)","");
        PionMcFile->Draw("EcalBDT_EnergyD>>PionMC_EcalBDT_EnergyD_withcut(100, -1, 1)","EcalBDT_EnergyD>-1");
        TH1F *PionMC_EcalBDT_EnergyD = (TH1F*)gDirectory->Get("PionMC_EcalBDT_EnergyD");
        TH1F *PionMC_EcalBDT_EnergyD_withcut = (TH1F*)gDirectory->Get("PionMC_EcalBDT_EnergyD_withcut");
        cout<< "Pion EcalBDT_EnergyD:"<< PionMC_EcalBDT_EnergyD_withcut->GetEntries() / PionMC_EcalBDT_EnergyD->GetEntries() <<endl;

        ElectronMcFile->Draw("EcalBDT_EnergyD>>ElectronMC_EcalBDT_EnergyD(100, -1, 1)","");
        ElectronMcFile->Draw("EcalBDT_EnergyD>>ElectronMC_EcalBDT_EnergyD_withcut(100, -1, 1)","EcalBDT_EnergyD>-1");
        TH1F *ElectronMC_EcalBDT_EnergyD = (TH1F*)gDirectory->Get("ElectronMC_EcalBDT_EnergyD");
        TH1F *ElectronMC_EcalBDT_EnergyD_withcut = (TH1F*)gDirectory->Get("ElectronMC_EcalBDT_EnergyD_withcut");
        cout<< "Electron EcalBDT_EnergyD:"<< ElectronMC_EcalBDT_EnergyD_withcut->GetEntries() / ElectronMC_EcalBDT_EnergyD->GetEntries() << endl;

        AntiprotonMcFile->Draw("EcalBDT_EnergyD>>AntiprotonMC_EcalBDT_EnergyD(100, -1, 1)","");
        AntiprotonMcFile->Draw("EcalBDT_EnergyD>>AntiprotonMC_EcalBDT_EnergyD_withcut(100, -1, 1)","EcalBDT_EnergyD>-1");
        TH1F *AntiprotonMC_EcalBDT_EnergyD = (TH1F*)gDirectory->Get("AntiprotonMC_EcalBDT_EnergyD");
        TH1F *AntiprotonMC_EcalBDT_EnergyD_withcut = (TH1F*)gDirectory->Get("AntiprotonMC_EcalBDT_EnergyD_withcut");
        cout<< "Antiproton EcalBDT_EnergyD:"<< AntiprotonMC_EcalBDT_EnergyD_withcut->GetEntries() / AntiprotonMC_EcalBDT_EnergyD->GetEntries()<<endl;

        // Build 1D TrdLogLikelihoodRatioProtonHeliumTracker Plot
        IssNegativeFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker(100, 0, 1.5)","");
        IssNegativeFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker_withcut(100, 0, 1.5)","TrdLogLikelihoodRatioProtonHeliumTracker>-1.5");
        TH1F *IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker = (TH1F*)gDirectory->Get("IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker");
        TH1F *IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker_withcut = (TH1F*)gDirectory->Get("IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker_withcut");

        PionMcFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>PionMC_TrdLogLikelihoodRatioProtonHeliumTracker(100, 0, 1.5)","");
        PionMcFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>PionMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut(100, 0, 1.5)","TrdLogLikelihoodRatioProtonHeliumTracker>-1.5");
        TH1F *PionMC_TrdLogLikelihoodRatioProtonHeliumTracker = (TH1F*)gDirectory->Get("PionMC_TrdLogLikelihoodRatioProtonHeliumTracker");
        TH1F *PionMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut = (TH1F*)gDirectory->Get("PionMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut");
        cout<< "Pion TrdLogLikelihoodRatioProtonHeliumTracker:"<< PionMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut->GetEntries() / PionMC_TrdLogLikelihoodRatioProtonHeliumTracker->GetEntries() <<endl;

        ElectronMcFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker(100, 0, 1.5)","");
        ElectronMcFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut(100, 0, 1.5)","TrdLogLikelihoodRatioProtonHeliumTracker>-1.5");
        TH1F *ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker = (TH1F*)gDirectory->Get("ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker");
        TH1F *ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut = (TH1F*)gDirectory->Get("ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut");
        cout<< "Electron TrdLogLikelihoodRatioProtonHeliumTracker:"<< ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut->GetEntries() / ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker->GetEntries() << endl;

        AntiprotonMcFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker(100, 0, 1.5)","");
        AntiprotonMcFile->Draw("TrdLogLikelihoodRatioProtonHeliumTracker>>AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut(100, 0, 1.5)","TrdLogLikelihoodRatioProtonHeliumTracker>-1.5");
        TH1F *AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker = (TH1F*)gDirectory->Get("AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker");
        TH1F *AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut = (TH1F*)gDirectory->Get("AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut");
        cout<< "Antiproton TrdLogLikelihoodRatioProtonHeliumTracker:"<< AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker_withcut->GetEntries() / AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker->GetEntries()<<endl;

     
        //// Scaling
        Double_t factor = 10000;

        IssNegativeData_TrdSegmentsXZNumber->Scale(factor/IssNegativeData_TrdSegmentsXZNumber->GetEntries());
        PionMC_TrdSegmentsXZNumber->Scale(factor/PionMC_TrdSegmentsXZNumber->GetEntries());
        ElectronMC_TrdSegmentsXZNumber->Scale(factor/ElectronMC_TrdSegmentsXZNumber->GetEntries());
        AntiprotonMC_TrdSegmentsXZNumber->Scale(factor/AntiprotonMC_TrdSegmentsXZNumber->GetEntries());

        IssNegativeData_TrdSegmentsYZNumber->Scale(factor/IssNegativeData_TrdSegmentsYZNumber->GetEntries());
        PionMC_TrdSegmentsYZNumber->Scale(factor/PionMC_TrdSegmentsYZNumber->GetEntries());
        ElectronMC_TrdSegmentsYZNumber->Scale(factor/ElectronMC_TrdSegmentsYZNumber->GetEntries());
        AntiprotonMC_TrdSegmentsYZNumber->Scale(factor/AntiprotonMC_TrdSegmentsYZNumber->GetEntries());

        IssNegativeData_TrdNumberOfHits->Scale(factor/IssNegativeData_TrdNumberOfHits->GetEntries());
        PionMC_TrdNumberOfHits->Scale(factor/PionMC_TrdNumberOfHits->GetEntries());
        ElectronMC_TrdNumberOfHits->Scale(factor/ElectronMC_TrdNumberOfHits->GetEntries());
        AntiprotonMC_TrdNumberOfHits->Scale(factor/AntiprotonMC_TrdNumberOfHits->GetEntries());

        IssNegativeData_TRDVTracksSize->Scale(factor/IssNegativeData_TRDVTracksSize->GetEntries());
        PionMC_TRDVTracksSize->Scale(factor/PionMC_TRDVTracksSize->GetEntries());
        ElectronMC_TRDVTracksSize->Scale(factor/ElectronMC_TRDVTracksSize->GetEntries());
        AntiprotonMC_TRDVTracksSize->Scale(factor/AntiprotonMC_TRDVTracksSize->GetEntries());

        IssNegativeData_ACCHits->Scale(factor/IssNegativeData_ACCHits->GetEntries());
        PionMC_ACCHits->Scale(factor/PionMC_ACCHits->GetEntries());
        ElectronMC_ACCHits->Scale(factor/ElectronMC_ACCHits->GetEntries());
        AntiprotonMC_ACCHits->Scale(factor/AntiprotonMC_ACCHits->GetEntries());

        IssNegativeData_EcalBDT_EnergyD->Scale(factor/IssNegativeData_EcalBDT_EnergyD->GetEntries());
        PionMC_EcalBDT_EnergyD->Scale(factor/PionMC_EcalBDT_EnergyD->GetEntries());
        ElectronMC_EcalBDT_EnergyD->Scale(factor/ElectronMC_EcalBDT_EnergyD->GetEntries());
        AntiprotonMC_EcalBDT_EnergyD->Scale(factor/AntiprotonMC_EcalBDT_EnergyD->GetEntries());

        IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker->Scale(factor/IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker->GetEntries());
        PionMC_TrdLogLikelihoodRatioProtonHeliumTracker->Scale(factor/PionMC_TrdLogLikelihoodRatioProtonHeliumTracker->GetEntries());
        ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker->Scale(factor/ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker->GetEntries());
        AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker->Scale(factor/AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker->GetEntries());


        //// Plot
        // Plot 1D TrdSegmentsXZNumber 
        Plot_1D_TrdSegmentsXZNumber(binning, i, AntiprotonMC_TrdSegmentsXZNumber, IssNegativeData_TrdSegmentsXZNumber, PionMC_TrdSegmentsXZNumber, ElectronMC_TrdSegmentsXZNumber);
        // Plot 1D TrdSegmentsYZNumber
        Plot_1D_TrdSegmentsYZNumber(binning, i, AntiprotonMC_TrdSegmentsYZNumber, IssNegativeData_TrdSegmentsYZNumber, PionMC_TrdSegmentsYZNumber, ElectronMC_TrdSegmentsYZNumber); 
        // Plot 1D TrdNumberOfHits 
        Plot_1D_TrdNumberOfHits(binning, i, AntiprotonMC_TrdNumberOfHits, IssNegativeData_TrdNumberOfHits, PionMC_TrdNumberOfHits, ElectronMC_TrdNumberOfHits);
        // Plot 1D TRDVTracksSize 
        Plot_1D_TRDVTracksSize(binning, i, AntiprotonMC_TRDVTracksSize, IssNegativeData_TRDVTracksSize, PionMC_TRDVTracksSize, ElectronMC_TRDVTracksSize);
        // Plot 1D ACCHits
        Plot_1D_ACCHits(binning, i, AntiprotonMC_ACCHits, IssNegativeData_ACCHits, PionMC_ACCHits, ElectronMC_ACCHits);
        // Plot 1D EcalBDT_EnergyD
        Plot_1D_EcalBDT_EnergyD(binning, i, AntiprotonMC_EcalBDT_EnergyD, IssNegativeData_EcalBDT_EnergyD, PionMC_EcalBDT_EnergyD, ElectronMC_EcalBDT_EnergyD);
        // Plot 1D TrdLogLikelihoodRatioProtonHeliumTracker
        Plot_1D_TrdLogLikelihoodRatioProtonHeliumTracker(binning, i, AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker, IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker, PionMC_TrdLogLikelihoodRatioProtonHeliumTracker, ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker);

        delete IssNegativeFile;
        delete PionMcFile;
        delete ElectronMcFile;
        delete AntiprotonMcFile;

    }

}
