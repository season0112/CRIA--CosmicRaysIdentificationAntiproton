// Useage: root -L /home/bo791269/Software/AntiprotonAnalysis/Libraries/AntiprotonBinning/AntiprotonBinning.C -x TempalteCheck_Intermediate.C
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

void TempalteCheck_Intermediate(){

    //// Load RichbetaCut
    vector<double> v_richbetacut;
    ifstream richbeta;
    richbeta.open( string("/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v5.0/total/RichBetaCutValue_eff_0.8_1.txt"));
    string cutvalue;
    assert(v_richbetacut.empty());
    while (getline(richbeta, cutvalue)) {
    v_richbetacut.push_back(stod(cutvalue));
    }
    richbeta.close();

    //// Load Data 
    vector<double> binning (AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili().Bins());

    TChain *felectron = new TChain("AntiprotonLowEnergyTree");
    TChain *fpass7_negative = new TChain("AntiprotonLowEnergyTree");
    TChain *fpion = new TChain("AntiprotonLowEnergyTree");

    int indexmin = 16;
    int indexmax = 17;

    for (auto i=indexmin; i<indexmax; i++){

        // Cut Defination
        TCut RichBetaCut = (string("RichBeta<") + to_string(v_richbetacut.at( (i-indexmin)/1 )) ).c_str();
        TCut PatternCut = ("Pattern==0 || Pattern==1 || Pattern==2 || Pattern==4");
        TCut TrdSegmentsXZNumberCut = "TrdSegmentsXZNumber==1";
        TCut TrdSegmentsYZNumberCut = "TrdSegmentsYZNumber==1";
        TCut TRDVTracksSizeCut = "TRDVTracksSize==1";
        TCut TrdNumberOfHitsCut = "TrdNumberOfHits<35";
        TCut ACCHitsCut = "ACCHits==0";


        // Templates selectron
        TCut ElectronTemplateDataCut = PatternCut && "EcalBDT_EnergyD>0.5" && "abs(RichBeta-(-Rigidity)/sqrt(0.000511^2+Rigidity^2))<0.002" && "TrdSegmentsXZNumber==1" && "TrdSegmentsYZNumber==1" && "TrdNumberOfHits<35" && "RichBeta>0"; //"RichBeta>0"?
        TCut PionTemplateDataCut = PatternCut && "TrdLogLikelihoodRatioElectronProtonTracker>-1.5 && TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>70 && TrdLogLikelihoodRatioProtonHeliumTracker<0.02 && TrdLogLikelihoodRatioProtonHeliumTracker>-1";


        // Print current Rigidity bin
        std::ostringstream lowedge_stream;
        lowedge_stream << binning.at(i);
        std::string lowedge = lowedge_stream.str();
        cout<< lowedge <<endl;

        std::ostringstream highedge_stream;
        highedge_stream << binning.at(i+1);
        std::string highedge = highedge_stream.str();
        cout<< highedge <<endl;


        // Load root files
        felectron->AddFile((string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v10.0/totalall/B1091_el.pl1.0_25200_7.6_all_Tree_negative_") + lowedge + string("_") + highedge + string(".root")).c_str());
        fpass7_negative->AddFile((string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v10.0/totalall/B1130_pass7_7.8_all_Tree_negative_May2015_") + lowedge + string("_") + highedge + string(".root")).c_str());
        fpion->AddFile((string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v10.0/totalall/B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + lowedge + string("_") + highedge + string(".root")).c_str());    


        // apply cut and draw histograms //

        // Electron histogram (MC)
        felectron->Draw("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_electron(50,-2.0, 0)", ElectronTemplateDataCut);
        TH1F *template_electron = (TH1F*)gDirectory->Get("th1f_electron");
        template_electron->SetTitle("Electron");

        // ISS negaitve data histogram
        fpass7_negative->Draw("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_ElectronDataTemplate(50,-2.0, 0)", ElectronTemplateDataCut);
        TH1F *template_electron_Data =  (TH1F*)gDirectory->Get("th1f_ElectronDataTemplate");
        template_electron_Data->SetTitle("Electron");    

        fpass7_negative->Draw("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_PionDataTemplate(50,-2.0, 0)", PionTemplateDataCut);
        TH1F *template_pion_Data =  (TH1F*)gDirectory->Get("th1f_PionDataTemplate");
        template_pion_Data->SetTitle("Pion");

        // Pion histogram (MC)
        fpion->Draw("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_pion(50, -2.0, 0)", PionTemplateDataCut);
        TH1F *template_pion =  (TH1F*)gDirectory->Get("th1f_pion");
        template_pion->SetTitle("Pion");


        //// Scaling
        Double_t factor = 1000;
        template_electron->Scale(factor/template_electron->GetEntries());
        template_electron_Data->Scale(factor/template_electron_Data->GetEntries());
        template_pion->Scale(factor/template_pion->GetEntries());
        template_pion_Data->Scale(factor/template_pion_Data->GetEntries());

        //// Plot 1D (Trdlikelihood) (Electron)
        TCanvas cplot_TRD_Electron("cplot_TRD_Electron","cplot_TRD_Electron",1000,500);
        /*
        IssNegativeDataTrd->SetTitle("TrdLogLikelihood (ISS Negative Rigidity Data)");
        IssNegativeDataTrd->GetXaxis()->SetTitle("#Lambda_{TRD}");
        IssNegativeDataTrd->GetXaxis()->SetTitleSize(0.05);
        IssNegativeDataTrd->GetXaxis()->SetTitleOffset(0.8);
        IssNegativeDataTrd->SetLineColor(kBlack);
        */
        template_electron->Draw("HIST");
        template_electron->SetLineColor(kRed);
        template_electron_Data->Draw("same");
        template_electron_Data->SetLineColor(kBlue);
        gPad->SetLogy();
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);
        TLegend* legend = new TLegend(0.12, 0.68, 0.38, 0.87, NULL, "brNDC");
        legend->SetFillColor(kWhite);
        legend->SetLineColor(kBlack);
        legend->AddEntry(template_electron, "Electron (MC)", "lp");
        legend->AddEntry(template_electron_Data, "Electron (Data)", "lp");
        legend->Draw();
        cplot_TRD_Electron.SaveAs("ElectronTemplate_Trdlikelihood.pdf");
        cplot_TRD_Electron.Close();

        //// Plot 1D (Trdlikelihood) (Pion)
        TCanvas cplot_TRD_Pion("cplot_TRD_Pion","cplot_TRD_Pion",1000,500);
        template_pion->Draw("HIST");
        template_pion->SetLineColor(kRed);
        template_pion_Data->Draw("same");
        template_pion_Data->SetLineColor(kBlue);
        gPad->SetLogy();
        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);
        TLegend* legend_pion = new TLegend(0.12, 0.68, 0.38, 0.87, NULL, "brNDC");
        legend_pion->SetFillColor(kWhite);
        legend_pion->SetLineColor(kBlack);
        legend_pion->AddEntry(template_pion, "Pion (MC)", "lp");
        legend_pion->AddEntry(template_pion_Data, "Pion (Data)", "lp");
        legend_pion->Draw();
        cplot_TRD_Pion.SaveAs("PionTemplate_Trdlikelihood.pdf");
        cplot_TRD_Pion.Close();

    }

}






