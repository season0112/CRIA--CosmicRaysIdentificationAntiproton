// Useage: root -L /home/bo791269/Software/AntiprotonAnalysis/Libraries/AntiprotonBinning/AntiprotonBinning.C -x PionTemplateSelection.C
#include <sstream>

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

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


void PionTemplateSelection(){
    //// Cuts and selections
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
    
    ///////////////////////////////////////////////////////////// 
    /////////////////// Cuts for templates //////////////////////
    ////////////////////////////////////////////////////////////
    TCut AntiprotonCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
    //TCut ElectronCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut ;
    //TCut ElectronCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut ;
    TCut ElectronCut;
    TCut PionCutMC = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && RichBetaCut && TriggerPhysics && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
    TCut PionCutMCNew = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>0.72";
    TCut NegativeCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;
    TCut PositiveCut = TrdLikelihoodCut && TrdLikelihoodHeProtonCut && ProtonCCMVABDCut && EcalBDT_EnergyDCut && RichBetaCut && TriggerPhysics && TofMassonecharge && ExtrapolatedPhotoElectronsFromTracker && RichPhotoElectron && RichChargeCut;

    TCut PionTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>0.72 && TrdNumberOfHits>42";
    TCut ElectronTemplateDataCut = "RichIsNaF==0 && abs(RichBeta-1)<0.002 && EcalBDT_EnergyD>-999 && EcalBDT_EnergyD>0.9"; 

    //////////////////////////////////////////////////////////
    //// Load data ///////////////////////////////////////////
    //////////////////////////////////////////////////////////
    vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());

    //i=0,i<3: 0.8-1.3GV;  i=3,i<6:1.33-1.92;  i=4,i<7,1.51-2.15;  i=4,i<5,1.51-1.71;
    int indexmin = 1;
    int indexmax = 16;

    for (auto i=indexmin; i<indexmax; i++){ 
        TChain *IssNegativeFile = new TChain("AntiprotonLowEnergyTree");
        TChain *PionMcFile = new TChain("AntiprotonLowEnergyTree");
        TChain *ElectronMcFile = new TChain("AntiprotonLowEnergyTree");
        TChain *AntiprotonMcFile = new TChain("AntiprotonLowEnergyTree");

        std::ostringstream lowedge_stream;
        lowedge_stream << binning.at(i);
        std::string lowedge = lowedge_stream.str();

        std::ostringstream highedge_stream;
        highedge_stream << binning.at(i+1);
        std::string highedge = highedge_stream.str();

        cout<< lowedge <<endl;
        cout<< highedge <<endl;

        double R_min = binning.at(i);
        double R_max = binning.at(i+1);

        //double R_min = binning.at(indexmin);
        //double R_max = binning.at(indexmax); 
        //cout << R_min << endl;
        //cout << R_max << endl;

        if ( lowedge == "1" || lowedge == "11" || lowedge == "12" || lowedge == "13" || lowedge == "18" ){
            lowedge = to_string_with_precision(binning[i], 1); }
        if ( highedge == "1" || highedge == "11" || highedge == "12" || highedge == "13" || highedge == "18" ){
            highedge = to_string_with_precision(binning[i+1], 1);}

        //IssNegativeFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1130_pass7_7.8_all_Tree_negative_") + lowedge + string("_") + highedge + string(".root")) .c_str() );    
        IssNegativeFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1130_pass7_7.8_all_Tree_negative_May2015_") + lowedge + string("_") + highedge + string(".root")) .c_str() );
        PionMcFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + lowedge + string("_") + highedge + string(".root")) .c_str() );
        ElectronMcFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1091_el.pl1.0_25200_7.6_all_Tree_negative_") + lowedge + string("_") +  highedge + string(".root")) .c_str() );
        AntiprotonMcFile->AddFile( ( string(getenv("HPCLOWENERGYDIR")) + string("/totalall/B1042_antipr.pl1.1800_7.6_all_Tree_") + lowedge + string("_") + highedge + string(".root")) .c_str() );
    

    /////////////////////////////////////////////////////////////////////////////////////////
    
    //// Build 1D TrdLikelihood Plot  (Before cuts, goal:show iss negative data components)
    IssNegativeFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker>>IssNegativeDataTrd(100, -0.2, 1.7)","");
    TH1F *IssNegativeDataTrd = (TH1F*)gDirectory->Get("IssNegativeDataTrd");

    PionMcFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker>>PionTrdMc(100, -0.2, 1.7)","");
    TH1F *PionTrdMc = (TH1F*)gDirectory->Get("PionTrdMc");

    ElectronMcFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker>>ElectronTrdMc(100, -0.2, 1.7)","");
    TH1F *ElectronTrdMc = (TH1F*)gDirectory->Get("ElectronTrdMc");

    AntiprotonMcFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker>>AntiprotonTrdMc(100, -0.2, 1.7)","");
    TH1F *AntiprotonTrdMc = (TH1F*)gDirectory->Get("AntiprotonTrdMc");
    

    //// Build 1D EcalBDT_EnergyD Plot  (Before cuts, goal:show iss negative data components)
    IssNegativeFile->Draw("EcalBDT_EnergyD>>IssNegativeDataEcalBDT_EnergyD(100, -1.0, 1.0)","");
    TH1F *IssNegativeDataEcalBDT_EnergyD = (TH1F*)gDirectory->Get("IssNegativeDataEcalBDT_EnergyD");

    PionMcFile->Draw("EcalBDT_EnergyD>>PionEcalBDT_EnergyDMc(100, -1.0, 1.0)","");
    TH1F *PionEcalBDT_EnergyDMc = (TH1F*)gDirectory->Get("PionEcalBDT_EnergyDMc");

    ElectronMcFile->Draw("EcalBDT_EnergyD>>ElectronEcalBDT_EnergyDMc(100, -1.0, 1.0)","");
    TH1F *ElectronEcalBDT_EnergyDMc = (TH1F*)gDirectory->Get("ElectronEcalBDT_EnergyDMc");     

    AntiprotonMcFile->Draw("EcalBDT_EnergyD>>AntiprotonEcalBDT_EnergyDMc(100, -1.0, 1.0)","");
    TH1F *AntiprotonEcalBDT_EnergyDMc = (TH1F*)gDirectory->Get("AntiprotonEcalBDT_EnergyDMc");   
    

    //// Build 1D template(1/TOFBeta)    (Beta = Z*R/m/gamma, Z=1, m_pion=139.57 MeV/c^2)
    IssNegativeFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) >>PionTemplateData_tofbeta(100, -0.6, 0.3)", PionTemplateDataCut);
    TH1F *PionTemplateData_tofbeta = (TH1F*)gDirectory->Get("PionTemplateData_tofbeta");

    IssNegativeFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) >>ElectronTemplateData_tofbeta(100, -0.6, 0.3)", ElectronTemplateDataCut);
    TH1F *ElectronTemplateData_tofbeta = (TH1F*)gDirectory->Get("ElectronTemplateData_tofbeta");

    PionMcFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) >>PionTemplateMcNew_tofbeta(100, -0.6, 0.3)", PionCutMCNew);
    TH1F *PionTemplateMcNew_tofbeta = (TH1F*)gDirectory->Get("PionTemplateMcNew_tofbeta");

    PionMcFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) >>PionTemplateMcOld_tofbeta(100, -0.6, 0.3)", PionCutMC);
    TH1F *PionTemplateMcOld_tofbeta = (TH1F*)gDirectory->Get("PionTemplateMcOld_tofbeta");

    ElectronMcFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) >>ElectronTemplateMc_tofbeta(100, -0.6, 0.3)", ElectronCut);
    TH1F *ElectronTemplateMc_tofbeta = (TH1F*)gDirectory->Get("ElectronTemplateMc_tofbeta");


    //// Build 1D template(Trdlikelihood)   
    IssNegativeFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker >>PionTemplateData_trdlikelihood(100, -0.2, 1.7)", PionTemplateDataCut);
    TH1F *PionTemplateData_trdlikelihood = (TH1F*)gDirectory->Get("PionTemplateData_trdlikelihood");

    IssNegativeFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker >>ElectronTemplateData_trdlikelihood(100, -0.2, 1.7)", ElectronTemplateDataCut);
    TH1F *ElectronTemplateData_trdlikelihood = (TH1F*)gDirectory->Get("ElectronTemplateData_trdlikelihood");

    PionMcFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker >>PionTemplateMcNew_trdlikelihood(100, -0.2, 1.7)", PionCutMCNew);
    TH1F *PionTemplateMcNew_trdlikelihood = (TH1F*)gDirectory->Get("PionTemplateMcNew_trdlikelihood");

    PionMcFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker >>PionTemplateMcOld_trdlikelihood(100, -0.2, 1.7)", PionCutMC);
    TH1F *PionTemplateMcOld_trdlikelihood = (TH1F*)gDirectory->Get("PionTemplateMcOld_trdlikelihood");

    ElectronMcFile->Draw("TrdLogLikelihoodRatioElectronProtonTracker >>ElectronTemplateMc_trdlikelihood(100, -0.2, 1.7)", ElectronCut);
    TH1F *ElectronTemplateMc_trdlikelihood = (TH1F*)gDirectory->Get("ElectronTemplateMc_trdlikelihood");


    //// Build 1D RichBeta Plot
    IssNegativeFile->Draw("RichBeta>>IssNegativeDataRichBeta(100, 0.94, 1.02)","RichIsNaF==0");
    TH1F *IssNegativeDataRichBeta = (TH1F*)gDirectory->Get("IssNegativeDataRichBeta");

    PionMcFile->Draw("RichBeta>>PionRichBetaMc(100, 0.94, 1.02)","RichIsNaF==0");
    TH1F *PionRichBetaMc = (TH1F*)gDirectory->Get("PionRichBetaMc");

    ElectronMcFile->Draw("RichBeta>>ElectronRichBetaMc(100, 0.94, 1.02)","RichIsNaF==0");
    TH1F *ElectronRichBetaMc = (TH1F*)gDirectory->Get("ElectronRichBetaMc");

    AntiprotonMcFile->Draw("RichBeta>>AntiprotonRichBetaMc(100, 0.94, 1.02)","RichIsNaF==0");
    TH1F *AntiprotonRichBetaMc = (TH1F*)gDirectory->Get("AntiprotonRichBetaMc");        


    //// Build 1D RichBeta-RichBeta_Pion Plot  (to show data composition)
    IssNegativeFile->Draw(" RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2)  >>IssNegativeDataRichBeta_PionAssumption(100, -0.02, 0.025)","RichIsNaF==0 && TrdLogLikelihoodRatioElectronProtonTracker>0.72");
    TH1F *IssNegativeDataRichBeta_PionAssumption = (TH1F*)gDirectory->Get("IssNegativeDataRichBeta_PionAssumption");
    
    PionMcFile->Draw(" RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2) >>PionRichBetaMc_PionAssumption(100, -0.02, 0.025)","RichIsNaF==0 && TrdLogLikelihoodRatioElectronProtonTracker>0.72");
    TH1F *PionRichBetaMc_PionAssumption = (TH1F*)gDirectory->Get("PionRichBetaMc_PionAssumption");

    ElectronMcFile->Draw(" RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2) >>ElectronRichBetaMc_PionAssumption(100, -0.02, 0.025)","RichIsNaF==0 && TrdLogLikelihoodRatioElectronProtonTracker>0.72");
    TH1F *ElectronRichBetaMc_PionAssumption = (TH1F*)gDirectory->Get("ElectronRichBetaMc_PionAssumption");


    //// Build 1D RichBeta-RichBeta_Electron Plot (to show data composition)
    IssNegativeFile->Draw(" RichBeta-1.0  >>IssNegativeDataRichBeta_ElectronAssumption(100, -0.02, 0.025)","RichIsNaF==0 && TrdLogLikelihoodRatioElectronProtonTracker<0.55");
    TH1F *IssNegativeDataRichBeta_ElectronAssumption = (TH1F*)gDirectory->Get("IssNegativeDataRichBeta_ElectronAssumption");

    PionMcFile->Draw(" RichBeta-1.0 >>PionRichBetaMc_ElectronAssumption(100, -0.02, 0.025)","RichIsNaF==0 && TrdLogLikelihoodRatioElectronProtonTracker<0.55");
    TH1F *PionRichBetaMc_ElectronAssumption = (TH1F*)gDirectory->Get("PionRichBetaMc_ElectronAssumption");

    ElectronMcFile->Draw(" RichBeta-1.0 >>ElectronRichBetaMc_ElectronAssumption(100, -0.02, 0.025)","RichIsNaF==0 && TrdLogLikelihoodRatioElectronProtonTracker<0.55");
    TH1F *ElectronRichBetaMc_ElectronAssumption = (TH1F*)gDirectory->Get("ElectronRichBetaMc_ElectronAssumption");


    //// Build 2D template (1/TOFBETA vs Rigidity)
    IssNegativeFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) : abs(Rigidity) >> PionTemplateData2D(100, R_min, R_max, 100, -0.6, 0.5)",PionTemplateDataCut);
    TH2F *PionTemplateData2D = (TH2F*)gDirectory->Get("PionTemplateData2D");

    IssNegativeFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) : abs(Rigidity) >>ElectronTemplateData2D(100, R_min, R_max, 100, -0.6, 0.5)",ElectronTemplateDataCut);
    TH2F *ElectronTemplateData2D = (TH2F*)gDirectory->Get("ElectronTemplateData2D");

    PionMcFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) : abs(Rigidity) >>PionTemplateMcNew2D(100, R_min, R_max, 100, -0.6, 0.5)",PionCutMCNew);
    TH2F *PionTemplateMcNew2D = (TH2F*)gDirectory->Get("PionTemplateMcNew2D");

    PionMcFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) : abs(Rigidity) >>PionTemplateMcOld2D(100, R_min, R_max, 100, -0.6, 0.5)", PionCutMC);
    TH2F *PionTemplateMcOld2D = (TH2F*)gDirectory->Get("PionTemplateMcOld2D");

    ElectronMcFile->Draw("1/TofBeta-1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2)))) : abs(Rigidity) >>ElectronTemplateMc2D(100, R_min, R_max, 100, -0.6, 0.5)", ElectronCut);
    TH2F *ElectronTemplateMc2D = (TH2F*)gDirectory->Get("ElectronTemplateMc2D");

 
    //// Scaling
    Double_t factor = 1000;
    PionTemplateData_tofbeta->Scale(factor/PionTemplateData_tofbeta->GetEntries());
    PionTemplateMcNew_tofbeta->Scale(factor/PionTemplateMcNew_tofbeta->GetEntries());
    PionTemplateMcOld_tofbeta->Scale(factor/PionTemplateMcOld_tofbeta->GetEntries());
    ElectronTemplateMc_tofbeta->Scale(factor/ElectronTemplateMc_tofbeta->GetEntries());
    ElectronTemplateData_tofbeta->Scale(factor/ElectronTemplateData_tofbeta->GetEntries());

    PionTemplateData_trdlikelihood->Scale(factor/PionTemplateData_trdlikelihood->GetEntries());
    PionTemplateMcNew_trdlikelihood->Scale(factor/PionTemplateMcNew_trdlikelihood->GetEntries());
    PionTemplateMcOld_trdlikelihood->Scale(factor/PionTemplateMcOld_trdlikelihood->GetEntries());
    ElectronTemplateMc_trdlikelihood->Scale(factor/ElectronTemplateMc_trdlikelihood->GetEntries());
    ElectronTemplateData_trdlikelihood->Scale(factor/ElectronTemplateData_trdlikelihood->GetEntries());


    PionEcalBDT_EnergyDMc->Scale(IssNegativeDataEcalBDT_EnergyD->GetEntries()*0.1/PionEcalBDT_EnergyDMc->GetEntries());
    ElectronEcalBDT_EnergyDMc->Scale(IssNegativeDataEcalBDT_EnergyD->GetEntries()*0.9/ElectronEcalBDT_EnergyDMc->GetEntries());
    //AntiprotonEcalBDT_EnergyDMc->Scale(IssNegativeDataEcalBDT_EnergyD->GetEntries()/AntiprotonEcalBDT_EnergyDMc->GetEntries());

    PionTrdMc->Scale(40*factor/PionTrdMc->GetEntries());
    ElectronTrdMc->Scale(400*factor/ElectronTrdMc->GetEntries());
    AntiprotonTrdMc->Scale(factor/AntiprotonTrdMc->GetEntries());

    IssNegativeDataRichBeta->Scale(factor/IssNegativeDataRichBeta->GetEntries());
    PionRichBetaMc->Scale(factor/PionRichBetaMc->GetEntries());
    ElectronRichBetaMc->Scale(factor/ElectronRichBetaMc->GetEntries());

    //IssNegativeDataRichBeta_PionAssumption->Scale(factor/IssNegativeDataRichBeta_PionAssumption->GetEntries());
    PionRichBetaMc_PionAssumption->Scale(factor*20/PionRichBetaMc_PionAssumption->GetEntries());
    ElectronRichBetaMc_PionAssumption->Scale(factor*12/ElectronRichBetaMc_PionAssumption->GetEntries());

    //IssNegativeDataRichBeta_ElectronAssumption->Scale(factor/IssNegativeDataRichBeta_ElectronAssumption->GetEntries());
    PionRichBetaMc_ElectronAssumption->Scale(factor*16/PionRichBetaMc_ElectronAssumption->GetEntries());
    ElectronRichBetaMc_ElectronAssumption->Scale(factor*150/ElectronRichBetaMc_ElectronAssumption->GetEntries());

    /*
    PionTemplateData2D->Scale(factor/PionTemplateData2D->GetEntries());
    PionTemplateMcNew2D->Scale(factor/PionTemplateMcNew2D->GetEntries());
    PionTemplateMcOld2D->Scale(factor/PionTemplateMcOld2D->GetEntries());
    ElectronTemplateMc2D->Scale(factor/ElectronTemplateMc2D->GetEntries());
    ElectronTemplateData2D->Scale(factor/ElectronTemplateData2D->GetEntries());
    */


    //// Plot 1D ISSMC compare (1/TOFBeta) 
    TCanvas cplot1D_TOFBETA_ISSMCCompare("cplot1D_TOFBETA_ISSMCCompare","cplot1D_TOFBETA_ISSMCCompare",1000,500);

    //ElectronTemplateData_tofbeta->SetMarkerStyle(7);
    ElectronTemplateData_tofbeta->SetMarkerColor(4);
    ElectronTemplateData_tofbeta->SetLineColor(4);
    ElectronTemplateData_tofbeta->SetTitle("");
    gPad->SetBottomMargin(0.18);
    ElectronTemplateData_tofbeta->GetXaxis()->SetTitle("#frac{1}{#bf{#beta}_{TOF}} - #frac{1}{#bf{#beta}(R,m_{p})}");
    ElectronTemplateData_tofbeta->GetXaxis()->SetTitleSize(0.05);
    ElectronTemplateData_tofbeta->GetXaxis()->SetTitleOffset(1.2);
    ElectronTemplateData_tofbeta->Draw("HIST");

    //PionTemplateData_tofbeta->SetMarkerStyle(2);
    PionTemplateData_tofbeta->SetMarkerColor(1);
    PionTemplateData_tofbeta->SetLineColor(1);
    PionTemplateData_tofbeta->Draw("HIST same");

    //PionTemplateMcNew_tofbeta->SetMarkerStyle(8);
    PionTemplateMcNew_tofbeta->SetMarkerColor(2);
    PionTemplateMcNew_tofbeta->SetLineColor(2);
    PionTemplateMcNew_tofbeta->Draw("HIST same");

    //ElectronTemplateMc_tofbeta->SetMarkerStyle(3);
    ElectronTemplateMc_tofbeta->SetMarkerColor(3);
    ElectronTemplateMc_tofbeta->SetLineColor(3);
    ElectronTemplateMc_tofbeta->Draw("HIST same");

    
    PionTemplateMcOld_tofbeta->SetMarkerStyle(7);
    PionTemplateMcOld_tofbeta->SetMarkerColor(10);
    PionTemplateMcOld_tofbeta->Draw("HIST same");
    

    TLegend * leg_tofbeta = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
    leg_tofbeta->AddEntry(PionTemplateData_tofbeta,"PionTemplateData","lp");
    leg_tofbeta->AddEntry(PionTemplateMcNew_tofbeta,"PionTemplateMcNew","lp");  // New way to select Pion MC. 
    leg_tofbeta->AddEntry(ElectronTemplateData_tofbeta,"ElectronTemplateData","lp");
    leg_tofbeta->AddEntry(ElectronTemplateMc_tofbeta,"ElectronTemplateMc","lp");
    leg_tofbeta->AddEntry(PionTemplateMcOld_tofbeta,"PionTemplateMcOld","lp");
    gStyle->SetLegendTextSize(0.03);
    leg_tofbeta->Draw();

    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gStyle->SetErrorX(0);
    cplot1D_TOFBETA_ISSMCCompare.SaveAs( (std::string("cplot1D_TOFBETA_ISSMCCompare") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );
    cplot1D_TOFBETA_ISSMCCompare.Close();


    //// Plot 1D ISSMC compare (TRDlikelihood)
    TCanvas cplot1D_TRDlikelihood_ISSMCCompare("cplot1D_TRDlikelihood_ISSMCCompare","cplot1D_TRDlikelihood_ISSMCCompare",1000,500);

    //ElectronTemplateData_trdlikelihood->SetMarkerStyle(7);
    ElectronTemplateData_trdlikelihood->SetMarkerColor(4);
    ElectronTemplateData_trdlikelihood->SetLineColor(4);
    ElectronTemplateData_trdlikelihood->SetTitle("");
    gPad->SetBottomMargin(0.18);
    ElectronTemplateData_trdlikelihood->GetXaxis()->SetTitle("#Lambda_{TRD}");
    ElectronTemplateData_trdlikelihood->GetXaxis()->SetTitleSize(0.05);
    ElectronTemplateData_trdlikelihood->GetXaxis()->SetTitleOffset(1.2);
    ElectronTemplateData_trdlikelihood->Draw("HIST");

    //PionTemplateData_trdlikelihood->SetMarkerStyle(2);
    PionTemplateData_trdlikelihood->SetMarkerColor(1);
    PionTemplateData_trdlikelihood->SetLineColor(1);
    PionTemplateData_trdlikelihood->Draw("HIST same");

    //PionTemplateMcNew_trdlikelihood->SetMarkerStyle(8);
    PionTemplateMcNew_trdlikelihood->SetMarkerColor(2);
    PionTemplateMcNew_trdlikelihood->SetLineColor(2);
    PionTemplateMcNew_trdlikelihood->Draw("HIST same");

    //ElectronTemplateMc_trdlikelihood->SetMarkerStyle(3);
    ElectronTemplateMc_trdlikelihood->SetMarkerColor(3);
    ElectronTemplateMc_trdlikelihood->SetLineColor(3);
    ElectronTemplateMc_trdlikelihood->Draw("HIST same");

    TLegend * leg_trd = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
    leg_trd->AddEntry(PionTemplateData_trdlikelihood,"PionTemplateData","lp");
    leg_trd->AddEntry(PionTemplateMcNew_trdlikelihood,"PionTemplateMc","lp"); // New Way to select.
    leg_trd->AddEntry(ElectronTemplateData_trdlikelihood,"ElectronTemplateData","lp");
    leg_trd->AddEntry(ElectronTemplateMc_trdlikelihood,"ElectronTemplateMc","lp");
    gStyle->SetLegendTextSize(0.03);
    leg_trd->Draw();

    gStyle->SetOptStat(0);
    gPad->SetLogy();
    gStyle->SetErrorX(0);
    cplot1D_TRDlikelihood_ISSMCCompare.SaveAs(  (std::string("cplot1D_TRDlikelihood_ISSMCCompare") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );
    cplot1D_TRDlikelihood_ISSMCCompare.Close();


    //// Plot 1D (Trdlikelihood)
    TCanvas cplot_issnegativedata("cplot_issnegativedata","cplot_issnegativedata",1000,500);
    IssNegativeDataTrd->SetTitle("TrdLogLikelihood (ISS Negative Rigidity Data)");
    IssNegativeDataTrd->GetXaxis()->SetTitle("#Lambda_{TRD}");
    IssNegativeDataTrd->GetXaxis()->SetTitleSize(0.05);
    IssNegativeDataTrd->GetXaxis()->SetTitleOffset(0.8);
    IssNegativeDataTrd->SetLineColor(kBlack);
    IssNegativeDataTrd->Draw("");
    PionTrdMc->SetLineColor(kRed);
    PionTrdMc->Draw("same HIST");
    ElectronTrdMc->SetLineColor(kBlue);
    ElectronTrdMc->Draw("same HIST");
    AntiprotonTrdMc->SetLineColor(kGreen);
    AntiprotonTrdMc->Draw("same HIST");
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    TLegend* legend = new TLegend(0.62, 0.68, 0.88, 0.87, NULL, "brNDC");
    legend->SetFillColor(kWhite);
    legend->SetLineColor(kBlack);
    legend->AddEntry(PionTrdMc, "Pion (MC)", "lp");
    legend->AddEntry(ElectronTrdMc, "Electron (MC)", "lp");
    legend->AddEntry(AntiprotonTrdMc, "Antiproton (MC)", "lp");
    legend->AddEntry(IssNegativeDataTrd, "ISS negative rigidity data", "lp");
    legend->Draw();
    cplot_issnegativedata.SaveAs( (std::string("IssNegativeData_Trdlikelihood") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );
    cplot_issnegativedata.Close();


    //// Plot 1D (EcalBDT_EnergyD)
    TCanvas cplot_issnegativedata_EcalBDT_EnergyD("cplot_issnegativedata_EcalBDT_EnergyD","cplot_issnegativedata_EcalBDT_EnergyD",1000,500);
    IssNegativeDataEcalBDT_EnergyD->SetTitle("TrdLogLikelihood (ISS Negative Rigidity Data)");
    IssNegativeDataEcalBDT_EnergyD->GetXaxis()->SetTitle("EcalBDT");
    IssNegativeDataEcalBDT_EnergyD->GetXaxis()->SetTitleSize(0.05);
    IssNegativeDataEcalBDT_EnergyD->GetXaxis()->SetTitleOffset(0.8);
    IssNegativeDataEcalBDT_EnergyD->SetLineColor(kBlack);
    IssNegativeDataEcalBDT_EnergyD->Draw("");
    PionEcalBDT_EnergyDMc->SetLineColor(kRed);
    PionEcalBDT_EnergyDMc->Draw("same HIST");
    ElectronEcalBDT_EnergyDMc->SetLineColor(kBlue);
    ElectronEcalBDT_EnergyDMc->Draw("same HIST");
    //AntiprotonEcalBDT_EnergyDMc->SetLineColor(kGreen);
    //AntiprotonEcalBDT_EnergyDMc->Draw("same HIST");
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);
    TLegend* legend_EcalBDT_EnergyD = new TLegend(0.62, 0.68, 0.88, 0.87, NULL, "brNDC");
    legend_EcalBDT_EnergyD->SetFillColor(kWhite);
    legend_EcalBDT_EnergyD->SetLineColor(kBlack);
    legend_EcalBDT_EnergyD->AddEntry(PionEcalBDT_EnergyDMc, "Pion (MC)", "lp");
    legend_EcalBDT_EnergyD->AddEntry(ElectronEcalBDT_EnergyDMc, "Electron (MC)", "lp");
    //legend_EcalBDT_EnergyD->AddEntry(AntiprotonEcalBDT_EnergyDMc, "Antiproton (MC)", "lp");
    legend_EcalBDT_EnergyD->AddEntry(IssNegativeDataEcalBDT_EnergyD, "ISS negative rigidity data", "lp");
    legend_EcalBDT_EnergyD->Draw();
    cplot_issnegativedata_EcalBDT_EnergyD.SaveAs( ( std::string("IssNegativeData_EcalBDT_EnergyD") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );
    cplot_issnegativedata_EcalBDT_EnergyD.Close();


    //// Plot 1D (RichBeta-PionAssumption)
    TCanvas cplot_richbeta_PionAssumption("cplot_richbeta_PionAssumption","cplot_richbeta_PionAssumption",1000,500);
    IssNegativeDataRichBeta_PionAssumption->SetTitle("RichBeta-_Beta_{Pion} (ISS Negative Rigidity Data)");
    IssNegativeDataRichBeta_PionAssumption->GetXaxis()->SetTitle("#beta_{RICH}-#beta_{Pion}");
    IssNegativeDataRichBeta_PionAssumption->GetXaxis()->SetTitleSize(0.05);
    IssNegativeDataRichBeta_PionAssumption->GetXaxis()->SetTitleOffset(0.8);
    IssNegativeDataRichBeta_PionAssumption->SetMarkerColor(1);
    IssNegativeDataRichBeta_PionAssumption->SetLineColor(1); //1 is black.
    IssNegativeDataRichBeta_PionAssumption->SetMarkerSize(0.04);
    IssNegativeDataRichBeta_PionAssumption->Draw("HIST");

    PionRichBetaMc_PionAssumption->SetMarkerStyle(50);
    PionRichBetaMc_PionAssumption->SetMarkerColor(kRed);
    PionRichBetaMc_PionAssumption->SetMarkerSize(0.04);
    PionRichBetaMc_PionAssumption->SetLineColor(kRed);
    PionRichBetaMc_PionAssumption->Draw("HIST same");

    ElectronRichBetaMc_PionAssumption->SetMarkerStyle(70);
    ElectronRichBetaMc_PionAssumption->SetMarkerColor(kBlue);
    ElectronRichBetaMc_PionAssumption->SetMarkerSize(0.04);
    ElectronRichBetaMc_PionAssumption->SetLineColor(kBlue);
    ElectronRichBetaMc_PionAssumption->Draw("HIST same");

    TLegend* legend_richbeta_PionAssumption = new TLegend(0.62, 0.68, 0.88, 0.87, NULL, "brNDC");
    legend_richbeta_PionAssumption->SetFillColor(kWhite);
    legend_richbeta_PionAssumption->SetLineColor(kBlack);
    legend_richbeta_PionAssumption->AddEntry(IssNegativeDataRichBeta_PionAssumption, "ISS negative rigidity data", "lp");
    legend_richbeta_PionAssumption->AddEntry(PionRichBetaMc_PionAssumption, "Pion (MC)", "lp");
    legend_richbeta_PionAssumption->AddEntry(ElectronRichBetaMc_PionAssumption, "Electron (MC)", "lp");
    legend_richbeta_PionAssumption->Draw();

    gPad->SetLogy();
    gStyle->SetErrorX(0);
    cplot_richbeta_PionAssumption.SaveAs( (std::string("IssNegativeData_RichBeta-PionAssumption")  + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );
    cplot_richbeta_PionAssumption.Close();


    //// Plot 1D (RichBeta-ElectronAssumption)
    TCanvas cplot_richbeta_ElectronAssumption("cplot_richbeta_ElectronAssumption","cplot_richbeta_ElectronAssumption",1000,500);
    IssNegativeDataRichBeta_ElectronAssumption->SetTitle("RichBeta-1.0 (ISS Negative Rigidity Data)");
    IssNegativeDataRichBeta_ElectronAssumption->GetXaxis()->SetTitle("#beta_{RICH}-1.0");
    IssNegativeDataRichBeta_ElectronAssumption->GetXaxis()->SetTitleSize(0.05);
    IssNegativeDataRichBeta_ElectronAssumption->GetXaxis()->SetTitleOffset(0.8);
    IssNegativeDataRichBeta_ElectronAssumption->SetMarkerColor(1);
    IssNegativeDataRichBeta_ElectronAssumption->SetLineColor(1); //1 is black.
    IssNegativeDataRichBeta_ElectronAssumption->SetMarkerSize(0.04);
    IssNegativeDataRichBeta_ElectronAssumption->Draw("HIST");

    PionRichBetaMc_ElectronAssumption->SetMarkerStyle(50);
    PionRichBetaMc_ElectronAssumption->SetMarkerColor(kRed);
    PionRichBetaMc_ElectronAssumption->SetMarkerSize(0.04);
    PionRichBetaMc_ElectronAssumption->SetLineColor(kRed);
    PionRichBetaMc_ElectronAssumption->Draw("HIST same");

    ElectronRichBetaMc_ElectronAssumption->SetMarkerStyle(70);
    ElectronRichBetaMc_ElectronAssumption->SetMarkerColor(kBlue);
    ElectronRichBetaMc_ElectronAssumption->SetMarkerSize(0.04);
    ElectronRichBetaMc_ElectronAssumption->SetLineColor(kBlue);
    ElectronRichBetaMc_ElectronAssumption->Draw("HIST same");

    TLegend* legend_richbeta_ElectronAssumption = new TLegend(0.62, 0.68, 0.88, 0.87, NULL, "brNDC");
    legend_richbeta_ElectronAssumption->SetFillColor(kWhite);
    legend_richbeta_ElectronAssumption->SetLineColor(kBlack);
    legend_richbeta_ElectronAssumption->AddEntry(IssNegativeDataRichBeta_ElectronAssumption, "ISS negative rigidity data", "lp");
    legend_richbeta_ElectronAssumption->AddEntry(PionRichBetaMc_ElectronAssumption, "Pion (MC)", "lp");
    legend_richbeta_ElectronAssumption->AddEntry(ElectronRichBetaMc_ElectronAssumption, "Electron (MC)", "lp");
    legend_richbeta_ElectronAssumption->Draw();

    gPad->SetLogy();
    gStyle->SetErrorX(0);
    cplot_richbeta_ElectronAssumption.SaveAs( (std::string("IssNegativeData_RichBeta-ElectronAssumption") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );
    cplot_richbeta_ElectronAssumption.Close();


    //// Plot 1D (RichBeta)
    TCanvas cplot_richbeta("cplot_richbeta","cplot_richbeta",1000,500);
    IssNegativeDataRichBeta->SetTitle("RichBeta (ISS Negative Rigidity Data)");
    IssNegativeDataRichBeta->GetXaxis()->SetTitle("#beta_{RICH}");
    IssNegativeDataRichBeta->GetXaxis()->SetTitleSize(0.05);
    IssNegativeDataRichBeta->GetXaxis()->SetTitleOffset(0.8);
    IssNegativeDataRichBeta->SetMarkerColor(1);
    IssNegativeDataRichBeta->SetLineColor(1); //1 is black.
    IssNegativeDataRichBeta->SetMarkerSize(0.04);
    IssNegativeDataRichBeta->Draw("HIST");

    PionRichBetaMc->SetMarkerStyle(50);
    PionRichBetaMc->SetMarkerColor(kRed);
    PionRichBetaMc->SetMarkerSize(0.04);
    PionRichBetaMc->SetLineColor(kRed);
    PionRichBetaMc->Draw("HIST same");

    ElectronRichBetaMc->SetMarkerStyle(70);
    ElectronRichBetaMc->SetMarkerColor(kBlue);
    ElectronRichBetaMc->SetMarkerSize(0.04);
    ElectronRichBetaMc->SetLineColor(kBlue);
    ElectronRichBetaMc->Draw("HIST same");

    TLegend* legend_richbeta = new TLegend(0.22, 0.68, 0.48, 0.87, NULL, "brNDC");
    legend_richbeta->SetFillColor(kWhite);
    legend_richbeta->SetLineColor(kBlack);
    legend_richbeta->AddEntry(IssNegativeDataRichBeta, "ISS negative rigidity data", "lp");
    legend_richbeta->AddEntry(PionRichBetaMc, "Pion (MC)", "lp");
    legend_richbeta->AddEntry(ElectronRichBetaMc, "Electron (MC)", "lp");
    legend_richbeta->Draw();

    gPad->SetLogy();
    gStyle->SetErrorX(0);
    cplot_richbeta.SaveAs( (std::string("IssNegativeData_RichBeta") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );
    cplot_richbeta.Close();



    //// Plot 2D (1/TOFBeta : Rigidity)
    TCanvas cplot_PionTemplateData2D("cplot_PionTemplateData2D","cplot_PionTemplateData2D",1000,500);
    PionTemplateData2D->Draw("COLZ");
    cplot_PionTemplateData2D.SaveAs( (std::string("PionTemplateData2D") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) +  std::string(".pdf")).c_str() );
    cplot_PionTemplateData2D.Close();

    TCanvas cplot_PionTemplateMcNew2D("cplot_PionTemplateMcNew2D","cplot_PionTemplateMcNew2D",1000,500);
    PionTemplateMcNew2D->Draw("COLZ");
    cplot_PionTemplateMcNew2D.SaveAs( (std::string("PionTemplateMcNew2D") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) +  std::string(".pdf")).c_str() ) ;
    cplot_PionTemplateMcNew2D.Close();

    TCanvas cplot_ElectronTemplateMc2D("cplot_ElectronTemplateMc2D","cplot_ElectronTemplateMc2D",1000,500);
    ElectronTemplateMc2D->Draw("COLZ");
    cplot_ElectronTemplateMc2D.SaveAs( (std::string("ElectronTemplateMc2D") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) +  std::string(".pdf")).c_str() );
    cplot_ElectronTemplateMc2D.Close();

    TCanvas cplot_ElectronTemplateData2D("cplot_ElectronTemplateData2D","cplot_ElectronTemplateData2D",1000,500);
    ElectronTemplateData2D->Draw("COLZ");
    cplot_ElectronTemplateData2D.SaveAs( (std::string("ElectronTemplateData2D") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) +  std::string(".pdf")).c_str() );
    cplot_ElectronTemplateData2D.Close();

    delete IssNegativeFile;
    delete PionMcFile;
    delete ElectronMcFile;
    delete AntiprotonMcFile;

    //PionTemplateData2D->Draw("COLZ same");
    //PionTemplateMcNew2D->Draw("COLZ same");
    //ElectronTemplateMc2D->Draw("COLZ same");

    /*
    TFile ratioroot("TemplatesPlots.root", "RECREATE");
    PionTemplateData2D->Write("PionTemplateData2D");
    PionTemplateMcNew2D->Write("PionTemplateMcNew2D");
    ElectronTemplateMc2D->Write("ElectronTemplateMc2D");
    ElectronTemplateData2D->Write("ElectronTemplateData2D");
    ratioroot.Close();
    */

    /*
    //// Fit in Effective_Acceptance_Ratio
    TF1  *function1 = new TF1("function1","[0]*x+[1]",1.51,1.71);
    eff_A->Fit(function1);
    TF1 *fittedfuction1 = eff_A->GetFunction("function1");
    TGraph *gFitFunction = new TGraph(fittedfuction1);
    */
    }

}
