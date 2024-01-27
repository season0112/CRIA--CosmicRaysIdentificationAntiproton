void DataMcComparision_TofBeta()
{
/*
TFile *Data = new TFile("B1130_pass7_7.8_all_Tree_positive_May2015_36.1_38.9.root");
TTree *treedata= (TTree*) Data->Get("ExampleAnalysisTree");
treedata->Draw("RichBeta>>test(500,0.96,1.05)","RichBeta>0.96");
TH1F *test = (TH1F*)gDirectory->Get("test");
Double_t factor = 1000;
test->Scale(factor/test->GetEntries());
test->SetMarkerStyle(8);
test->Draw("");

TFile *Mc = new TFile("B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_36.1_38.9.root");
TTree *treemc= (TTree*) Mc->Get("ExampleAnalysisTree");
treemc->Draw("RichBeta>>test2(500,0.96,1.05) same","RichBeta>0.96");
TH1F *test2 = (TH1F*)gDirectory->Get("test2");
test2->Scale(factor/test2->GetEntries());
test2->SetMarkerStyle(7);
test2->Draw("");

test->Draw("same");
*/


TCanvas TofBetaplot("TofBetaplot","TofBetaplot",1000,500);

TFile *Data = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1130_pass7_7.8_all_Tree_positive_May2015_1_1.16.root");
TTree *treedata= (TTree*) Data->Get("AntiprotonLowEnergyTree");
treedata->Draw("TofBeta>>test(500,0.50,1.05)","TofBeta>0.50");
TH1F *test = (TH1F*)gDirectory->Get("test");
Double_t factor = 1000;
test->Scale(factor/test->GetEntries());
test->SetMarkerStyle(8);
test->SetMarkerColor(kBlue);
test->Draw("");

TFile *Mc = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v9.0/totalall/B1042_antipr.pl1.1800_7.6_all_Tree_1_1.16.root");
TTree *treemc= (TTree*) Mc->Get("AntiprotonLowEnergyTree");
treemc->Draw("TofBeta>>test2(500,0.50,1.05) same","TofBeta>0.50");
TH1F *test2 = (TH1F*)gDirectory->Get("test2");
test2->Scale(factor/test2->GetEntries());
test2->SetMarkerStyle(7);
test2->SetMarkerColor(kRed);
test2->Draw("");

test->Draw("same");

TofBetaplot.SaveAs("TofBetaplot.pdf");





return 0;
}
