// ACsoft includes
#include "AnalysisEvent.hh"
#include "AcceptanceManager.hh"
#include "AntiprotonBinning.hh"
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
#include "Utilities.hh"
#include "ValueHistograms.hh"
#include <string>

// ROOT includes
#include <math.h>
#include <TFile.h>
#include <TH1.h>
#include "TreeFormula.hh"
#include "Utilities.hh"
#include <TROOT.h>
#include "TreeWriter.hh"
#include <TApplication.h>
#include <TProof.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TAttMarker.h>
#include <TAttLine.h>
#include "TStyle.h"
#include "TColor.h"
#include "TF2.h"
#include "TExec.h"
#include "TCanvas.h"


#define INFO_OUT_TAG "FluxCalculation"
#include "debugging.hh"

//// Produce flux, right now only on 16.6-259GV, last two points binnings are different, need to be correctd. 
//// Proton Published flux result: because bin different, value can use "eval", but error  above 80 needs to be correct.

// low range 1-3.29GV: 1, 1.16, 1.33, 1.51, 1.71, 1.92, 2.15, 2.4, 2.67, 2.97, 3.29. 10 bins in total.
// // tomorary solution for 2.97-16.6GV flux, 2.97, 3.29, 3.64, 4.02, 4.43, 4.88, 5.37, 5.9, 6.47, 7.09, 7.76, 8.48, 9.26, 10.1, 11, 12, 13, 14.1, 15.3, 16.6. 19 bins in total.
// // 2.97-4.43 GV v4.0 cut0.985 (4 points), 4.43- 5.37 GV v4.0 cut 0.990 (2 points),  5.37-7.76 GV V1.0 cut 0.997(4 points),

// Error has some problems?

int main(int argc, char** argv) {

  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);

  config.SetProgramHelpText("FluxCalculation",
                            "Calculate Flux");

  config.AddHelpExample("FluxCalculation", "--binningversion 450version --issversion pass7.8");

  std::string binningversion = "";
  config.GetValue("OPTIONS", "binningversion", binningversion,
                  "The binningversion is");

  std::string issversion = "";
  config.GetValue("OPTIONS", "issversion", issversion,
                  "The issversion is:");

  if (binningversion == "") {
    WARN_OUT << "No binningversion is given! Please give a binningversion." << std::endl;
    return EXIT_FAIL_CONFIG;
  }

  if (issversion == "") {
    WARN_OUT << "No issversion is given! Please give a issversion." << std::endl;
    return EXIT_FAIL_CONFIG;
  }

chdir( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/ISS_anylsis/unfolding/") + std::string(issversion) + std::string("/") + std::string(binningversion) ).c_str());

/////////////////////////////
//// Load Files  ////////////
////////////////////////////
TFile *f0 = new TFile("ratio_TH2D.root");
TH1D *hproton_number = (TH1D*)f0->Get("proton_number");
TH1D *hantiproton_number = (TH1D*)f0->Get("antiproton_number");

TFile *f1 = new TFile("unfolded_results.root");
TH1D *hAntiproton_unfolded = (TH1D*)f1->Get("hAntiproton_unfolded");
TH1D *hProton_unfolded = (TH1D*)f1->Get("hProton_unfolded");
TH1D *hProton_raw = (TH1D*)f1->Get("hProton_raw");
TH1D *hAntiproton_raw = (TH1D*)f1->Get("hAntiproton_raw");
TH1D *Statistic_error = (TH1D*)f1->Get("Statistic_error");
TH1D *Delta_antiproton = (TH1D*)f1->Get("Delta_antiproton");

TFile *f2 = new TFile( (std::string("/hpcwork/jara0052/sichen/Measuringtime/MeasuringTime_") + std::string(issversion) + std::string(".root")).c_str() );
TH1D *hIntegratedMeasuringTimeOverCutOff = (TH1D*)f2->Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff");

TFile *f3 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_") + std::string(binningversion) + std::string(".root")).c_str() );   //  old:EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_v7.0.root
TGraphAsymmErrors *geffectiveAcceptanceAfterAllCuts = (TGraphAsymmErrors*)f3->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
TFile *f7 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_") + std::string(binningversion) + std::string(".root")).c_str() );   //  old:EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_v7.0.root
TGraphAsymmErrors *geffectiveAcceptanceAfterAllCuts_proton = (TGraphAsymmErrors*)f7->Get("QualityCuts/effectiveAcceptanceAfterAllCuts"); 

TFile *f6 = new TFile( (std::string(getenv("MY_ANALYSIS")) + std::string("/ReferenceFiles/AntiprotonFlux/AntiprotonFluxMIT2018.root")).c_str() );
TGraphAsymmErrors *AntiprotonFluxMIT2018 = (TGraphAsymmErrors*)f6->Get("g_Flux_pbar");

//TFile *f4 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/TriggerEff_B1042_antipr.pl1.1800_7.6_all.root")).c_str() );
  TFile *f4 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/TriggerEff_B1042_antipr.pl1.1800_7.6_all_") + std::string(binningversion)  + std::string(".root")).c_str() );
//TFile *f4 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/TriggerEff_old.root")).c_str() );

TH1F *TriggerEff = (TH1F*)f4->Get("TriggerEff");
TH1F *TriggerEff_noprescaling = (TH1F*)f4->Get("TriggerEff_noprescaling");
TGraphErrors gTriggerEff_noprescaling = TGraphErrors(TriggerEff_noprescaling);
for (int s=1; s <= TriggerEff_noprescaling->GetNbinsX(); s++){
    gTriggerEff_noprescaling.SetPointError(s-1, TriggerEff_noprescaling->GetBinWidth(s)/2, 0.0);
}

/////////////////////////////
////// Published result//////
/////////////////////////////
std::vector<double> subrange_450( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+27, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+59); //from 14.1

///////////////////////////
//////  Some Set Up ///////
///////////////////////////
TH1D *hAntiproton_unfolded_withspectrumindex = new TH1D(*hAntiproton_unfolded); 
TH1D *hProton_unfolded_withspectrumindex = new TH1D(*hProton_unfolded); 
TH1D *hProton_raw_withspectrumindex = new TH1D(*hProton_unfolded); 
TH1D *hIntegratedMeasuringTimeOverCutOff_shortcut = new TH1D(*hAntiproton_unfolded); 
TH1D *heffectiveAcceptanceAfterAllCuts_shortcut = new TH1D(*hAntiproton_unfolded);
TH1D *heffectiveAcceptanceAfterAllCuts_shortcut_proton = new TH1D(*hAntiproton_unfolded); 
TH1D *heffectiveAcceptanceAfterAllCutsError_shortcut = new TH1D(*hAntiproton_unfolded); 
TH1D *heffectiveAcceptanceAfterAllCutsError_shortcut_proton = new TH1D(*hAntiproton_unfolded); 
TH1D *hTriggerEff_shortcut = new TH1D(*hAntiproton_unfolded); 
TH1D *hhrigiditybin = new TH1D(*hAntiproton_unfolded); 
auto binning = AntiprotonNewBinning::NewBinning::AntiprotonBinning525_zhili();
TH1D *hStatistical = new TH1D(*hAntiproton_unfolded); 
TH1D *hStatistical_proton = new TH1D(*hAntiproton_unfolded); 

for (int k = 1; k <= 32; ++k){
  hAntiproton_unfolded_withspectrumindex->SetBinContent(k, hAntiproton_unfolded->GetBinContent(k) * pow(hAntiproton_unfolded->GetXaxis()->GetBinCenter(k), 2.7));
  hProton_unfolded_withspectrumindex->SetBinContent(k,hProton_unfolded->GetBinContent(k) * pow(hProton_unfolded->GetXaxis()->GetBinCenter(k), 2.7));
  hProton_raw_withspectrumindex->SetBinContent(k,hProton_raw->GetBinContent(k) * pow(hProton_unfolded->GetXaxis()->GetBinCenter(k), 2.7));
  hIntegratedMeasuringTimeOverCutOff_shortcut->SetBinContent(k, hIntegratedMeasuringTimeOverCutOff->GetBinContent(k+27)); //i+29:from 16.6, i+27:from 14.1
  hTriggerEff_shortcut->SetBinContent(k, TriggerEff_noprescaling->GetBinContent(k+27)); // i+29:from 16.6, i+27:from 14.1
  heffectiveAcceptanceAfterAllCuts_shortcut->SetBinContent(k, geffectiveAcceptanceAfterAllCuts->Eval(hAntiproton_unfolded->GetXaxis()->GetBinCenter(k)));
  heffectiveAcceptanceAfterAllCuts_shortcut_proton->SetBinContent(k, geffectiveAcceptanceAfterAllCuts_proton->Eval(hProton_unfolded->GetXaxis()->GetBinCenter(k)));
  heffectiveAcceptanceAfterAllCutsError_shortcut->SetBinContent(k, geffectiveAcceptanceAfterAllCuts->GetY()[k+26]); //index29 is 17.3(16.6-18.0)
  heffectiveAcceptanceAfterAllCutsError_shortcut_proton->SetBinContent(k, geffectiveAcceptanceAfterAllCuts_proton->GetY()[k+26]); //index29 is 17.3(16.6-18.0)
  hhrigiditybin->SetBinContent(k, binning.BinWidth(k+27)); // to check
  hStatistical->SetBinContent(k,  sqrt(pow((heffectiveAcceptanceAfterAllCutsError_shortcut->GetBinContent(k)*10000 / pow(heffectiveAcceptanceAfterAllCuts_shortcut->GetBinContent(k)*10000,2) * hAntiproton_unfolded->GetBinContent(k) / hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(k) / hhrigiditybin->GetBinWidth(k) * pow(hAntiproton_unfolded->GetXaxis()->GetBinCenter(k), 2.7)),2) + pow((Delta_antiproton->GetBinContent(k) / heffectiveAcceptanceAfterAllCuts_shortcut->GetBinContent(k)*10000 / hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(k) / hhrigiditybin->GetBinWidth(k)),2)));
  hStatistical_proton->SetBinContent(k, hProton_unfolded->GetBinContent(k) /hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(k) /hhrigiditybin->GetBinWidth(k) / pow(heffectiveAcceptanceAfterAllCuts_shortcut_proton->GetBinContent(k),2) * heffectiveAcceptanceAfterAllCutsError_shortcut_proton->GetBinContent(k) * pow(hAntiproton_unfolded->GetXaxis()->GetBinCenter(k), 2.7 /1000) );
}

////////////////////////////////////////////////////////
/////////////////       Flux Results       /////////////
////////////////////////////////////////////////////////
// 1. My Antiproton Flux: N/T/E/Rbin*CCapplicable  (deltaA is sqrt(A) for now)

TGraphErrors *gflux = new TGraphErrors();
TH1D *hStatistical_relative = new TH1D(*hAntiproton_unfolded); //32 points 14.1-525.
for (int l=1; l <= 32; ++l){
//    gflux->SetPoint(l-1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_525.at(29+l-1), hAntiproton_unfolded_withspectrumindex->GetBinContent(l) / hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(l) / heffectiveAcceptanceAfterAllCuts_shortcut->GetBinContent(l) *10000 / hhrigiditybin->GetBinWidth(l) / hTriggerEff_shortcut->GetBinContent(l) );
    gflux->SetPoint(l-1, (hAntiproton_unfolded_withspectrumindex->GetBinLowEdge(l) + hAntiproton_unfolded_withspectrumindex->GetBinLowEdge(l+1))/2, hAntiproton_unfolded_withspectrumindex->GetBinContent(l) / hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(l) / heffectiveAcceptanceAfterAllCuts_shortcut->GetBinContent(l) *10000 / hhrigiditybin->GetBinWidth(l) / hTriggerEff_shortcut->GetBinContent(l) );
    std::cout<< "l"  <<std::endl;
    /*
    std::cout<<hAntiproton_unfolded_withspectrumindex->GetBinContent(l) <<std::endl;
    std::cout<< hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(l) <<std::endl;
    std::cout<< heffectiveAcceptanceAfterAllCuts_shortcut->GetBinContent(l)<<std::endl;
    std::cout<< hhrigiditybin->GetBinWidth(l)<<std::endl;
    std::cout<< hTriggerEff_shortcut->GetBinContent(l)<<std::endl;
    */
    std::cout<< hAntiproton_unfolded_withspectrumindex->GetBinContent(l) / hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(l) / heffectiveAcceptanceAfterAllCuts_shortcut->GetBinContent(l) *10000 / hhrigiditybin->GetBinWidth(l) / hTriggerEff_shortcut->GetBinContent(l) <<std::endl;
    gflux->SetPointError(l-1, 0.0, hStatistical->GetBinContent(l)*pow(hAntiproton_unfolded->GetXaxis()->GetBinCenter(l), 2.7) );
    hStatistical_relative->SetBinContent(l, hStatistical->GetBinContent(l) *pow(hAntiproton_unfolded->GetXaxis()->GetBinCenter(l), 2.7)  /gflux->GetY()[l-1] *100);
}


// 2. My Proton Flux: N/T/E/Rbin*CCapplicable (deltaA is sqrt(A) for now)
//ProtonFlux is devided by 1000 to show in same plot. eventcount is without R**2.7.
TGraphErrors *gprotonflux = new TGraphErrors(); // TGraphErrors constructor  from THistogram : th1d Xaxis and bincontent index 1 becomes TGraphErrors Xaxis index 0.  see iterator.
for (int c=1; c <= 32; c++){
    gprotonflux->SetPoint(c-1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_525.at(27+c-1), hProton_unfolded_withspectrumindex->GetBinContent(c) / hIntegratedMeasuringTimeOverCutOff_shortcut->GetBinContent(c) / heffectiveAcceptanceAfterAllCuts_shortcut_proton->GetBinContent(c) *10000 / hhrigiditybin->GetBinWidth(c) / hTriggerEff_shortcut->GetBinContent(c) / 1000);
    gprotonflux->SetPointError(c-1, 0.0, hStatistical_proton->GetBinContent(c) );
}


//// 3. Antiproton Flux 2016 Paper
TGraphErrors *gantiprotonfluxpapaer_withspectrumindex = new TGraphErrors(); // TGraphErrors constructor  from THistogram : th1d Xaxis and bincontent index 1 becomes TGraphErrors Xaxis index 0.  see iterator.
TH1D *PublishedFAntiprotonluxError = new TH1D("","", 31, subrange_450.data());
TH1D *PublishedFAntiprotonluxStatisticError = new TH1D("","", 31, subrange_450.data());
TH1D *PublishedFAntiprotonluxStatisticError_relative = new TH1D("","", 31, subrange_450.data());

for (int q = 1; q <= 31; ++q){
    gantiprotonfluxpapaer_withspectrumindex->SetPoint(q-1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.at(27+q-1), AntiprotonNewBinning::AntiprotonResults::PublishedFluxPRL.at(q-1)  * pow( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.at(27+q-1) , 2.7) ); 
    gantiprotonfluxpapaer_withspectrumindex->SetPointError(q-1, 0.0, AntiprotonNewBinning::AntiprotonResults::PublishedFluxErrorPRL.at(q-1)* pow( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinningCenter_450.at(27+q-1) , 2.7) );
    PublishedFAntiprotonluxError->SetBinContent(q,  AntiprotonNewBinning::AntiprotonResults::PublishedFluxErrorPRL.at(q-1));
    PublishedFAntiprotonluxStatisticError->SetBinContent(q, AntiprotonNewBinning::AntiprotonResults::PublishedFluxStatisticErrorPRL.at(q-1));
    PublishedFAntiprotonluxStatisticError_relative->SetBinContent(q, AntiprotonNewBinning::AntiprotonResults::PublishedFluxStatisticRelativeErrorPRL.at(q-1));
}


// 4.  Proton Flux 2015 Paper
TGraphErrors *gprotonfluxpapaer_withspectrumindex = new TGraphErrors(); // TGraphErrors constructor  from THistogram : th1d Xaxis and bincontent index 1 becomes TGraphErrors Xaxis index 0.  see iterator.  from 16.6
for (int e = 1; e <= 45; ++e){
    gprotonfluxpapaer_withspectrumindex->SetPoint(e-1, Binning::Predefined::AmsProtonPaper2015Binning().Value(28+e-1), AntiprotonNewBinning::AntiprotonResults::PublishedProtonFluxPRL.at(26+e-1)* pow( Binning::Predefined::AmsProtonPaper2015Binning().Value(28+e-1),2.7)/1000 );
    gprotonfluxpapaer_withspectrumindex->SetPointError(e-1, 0.0, AntiprotonNewBinning::AntiprotonResults::PublishedProtonFluxErrorPRL.at(26+e-1) * pow( Binning::Predefined::AmsProtonPaper2015Binning().Value(28+e-1) , 2.7) /pow(1000,2.7) );
}


// 5. Antiprton Flux compare
TGraph *flux_diff =new TGraph();
for (int e = 1; e <= 30; ++e){   
    flux_diff->SetPoint(e-1, gflux->GetX()[e-1], (gflux->GetY()[e-1]-gantiprotonfluxpapaer_withspectrumindex->GetY()[e-1])/gflux->GetY()[e-1]);
}


//////////////////////////////
//////        plot        ////
//////////////////////////////
TCanvas *c1 = new TCanvas();
gStyle->SetOptStat("00000000");
hIntegratedMeasuringTimeOverCutOff_shortcut->Draw("HIST");
gPad->SetLogx();
hIntegratedMeasuringTimeOverCutOff_shortcut->GetXaxis()->SetMoreLogLabels();
c1->SaveAs("plot_flux/MeasuringTime.pdf");

TCanvas *c2 = new TCanvas();
gStyle->SetOptStat("00000000");
heffectiveAcceptanceAfterAllCuts_shortcut->Draw("HIST");
gPad->SetLogx();
heffectiveAcceptanceAfterAllCuts_shortcut->GetXaxis()->SetMoreLogLabels();
c2->SaveAs("plot_flux/EffectiveAcceptance.pdf");

TCanvas *c22 = new TCanvas();
gStyle->SetOptStat("00000000");
gTriggerEff_noprescaling.Draw("AP *");
gTriggerEff_noprescaling.GetXaxis()->SetLimits(4,700);
gTriggerEff_noprescaling.GetYaxis()->SetLimits(0.5,1.0);
gTriggerEff_noprescaling.GetYaxis()->SetRangeUser(0.5,1.0);
gPad->SetLogx();
gTriggerEff_noprescaling.GetXaxis()->SetMoreLogLabels();
gTriggerEff_noprescaling.GetXaxis()->SetTitle("|Rigidity| (GV)");
gTriggerEff_noprescaling.GetYaxis()->SetTitle("Trigger Efficiency");
c22->SaveAs("plot_flux/TriggerEff_noprescaling.pdf");

TCanvas *c3 = new TCanvas();
TPad *p32 = new TPad("p32","p32",0.,0.,1.,0.3);
p32->Draw();
p32->SetTopMargin(0.001);
p32->SetBottomMargin(0.3);
p32->SetLogx();
TPad *p1 = new TPad("p1","p1",0.,0.3,1.,1.0);
p1->Draw();
p1->SetBottomMargin(0.001);
p1->cd();
p1->SetLogx();
gflux->Draw("AP");
gantiprotonfluxpapaer_withspectrumindex->Draw("same *");
gantiprotonfluxpapaer_withspectrumindex->SetMarkerStyle(kStar);
gflux->SetMarkerStyle(kCircle);
gflux->SetMarkerColor(kRed);
gantiprotonfluxpapaer_withspectrumindex->SetMarkerColor(kBlue);
gflux->SetLineColor(kRed);
gantiprotonfluxpapaer_withspectrumindex->SetLineColor(kBlue);
gflux->GetXaxis()->SetMoreLogLabels();
gflux->GetXaxis()->SetTitle("|Rigidity| (GV)");
gflux->GetYaxis()->SetTitle("#Phi  R^{2.7}");
gflux->GetXaxis()->SetLimits(31.1,525);
gflux->GetYaxis()->SetRangeUser(0.1,4.7);
TLegend *legend = new TLegend(0.50,0.70,0.78,0.84);
legend->AddEntry(gflux,"This analysis (Statistical Error only)","lpf");
legend->AddEntry(gantiprotonfluxpapaer_withspectrumindex,"2016PRL paper","lpf");
legend->Draw();
p32->cd();
flux_diff->GetXaxis()->SetMoreLogLabels();
flux_diff->GetXaxis()->SetLimits(31.1,525);
flux_diff->SetTitle("");
flux_diff->GetXaxis()->SetTitleSize(0.4);
flux_diff->GetXaxis()->SetLabelSize(0.09);
flux_diff->GetYaxis()->SetLabelSize(0.09);
flux_diff->Draw("AP *");
c3->SaveAs("plot_flux/Flux.pdf");

TCanvas *c31 = new TCanvas();
hStatistical->Draw("hist");
PublishedFAntiprotonluxStatisticError->Draw("same");
c31->SaveAs("plot_flux/StatisticalError.pdf");

TCanvas *c32 = new TCanvas();
hStatistical_relative->Draw("hist");
PublishedFAntiprotonluxStatisticError_relative->Draw("same");
TLegend *legend_sta_rela = new TLegend(0.20,0.70,0.40,0.84);
legend_sta_rela->AddEntry(hStatistical_relative, "This analysis","lpf");
legend_sta_rela->AddEntry(PublishedFAntiprotonluxStatisticError_relative,"2016PRL paper","lpf");
legend_sta_rela->Draw();
hStatistical_relative->SetLineColor(1);
PublishedFAntiprotonluxStatisticError_relative->SetLineColor(2);
hStatistical_relative->GetXaxis()->SetTitle("|Rigidity| (GV)");
hStatistical_relative->GetYaxis()->SetTitle("Relative Statistical Error");
c32->SaveAs("plot_flux/StatisticalError_relative.pdf");


TCanvas *c4 = new TCanvas();
TPad *p2 = new TPad();
gPad->SetLogx();
gprotonflux->Draw("AP");
gprotonfluxpapaer_withspectrumindex->Draw("same *");
gprotonflux->SetMarkerStyle(kCircle);
gprotonfluxpapaer_withspectrumindex->SetMarkerStyle(kStar);
gprotonflux->SetMarkerColor(kRed);
gprotonfluxpapaer_withspectrumindex->SetMarkerColor(kBlue);
gprotonflux->SetLineColor(kRed);
gprotonfluxpapaer_withspectrumindex->SetLineColor(kBlue);
gprotonflux->GetXaxis()->SetMoreLogLabels();
gprotonflux->GetXaxis()->SetTitle("|Rigidity| (GV)");
gprotonflux->GetYaxis()->SetTitle("#Phi#upointR^{2.7}#upoint10^{3}");
gprotonflux->GetYaxis()->SetTitleSize(0.04);
gprotonflux->GetYaxis()->SetRangeUser(0.0,26.0);
TLegend *legend2 = new TLegend(0.45,0.63,0.78,0.84);
legend2->AddEntry(gprotonflux,"This analysis (Statistical Error only)","lpf");
legend2->AddEntry(gprotonfluxpapaer_withspectrumindex,"2015 PRL Paper","lpf");
legend2->Draw();
p2->Update();
c4->SaveAs("plot_flux/ProtonFlux.pdf");


/////////////////////////////
////// Result File ///////////
/////////////////////////////
TFile *results = new TFile("flux.root","RECREATE");
gantiprotonfluxpapaer_withspectrumindex->Write("gantiprotonfluxpapaer_withspectrumindex");
gprotonfluxpapaer_withspectrumindex->Write("gprotonfluxpapaer_withspectrumindex");
hStatistical_relative->Write("hStatistical_relative");
heffectiveAcceptanceAfterAllCuts_shortcut->Write("heffectiveAcceptanceAfterAllCuts_shortcut");
results->Close();





  return EXIT_SUCCESS;
}




