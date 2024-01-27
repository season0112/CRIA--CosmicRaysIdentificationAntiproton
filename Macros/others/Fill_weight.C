void Fill_weight(){

TFile *ISSdata = new TFile("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/B1130_pass7_7.8_all_Tree_positive_175_211.root");
TFile *B1220proton = new TFile("/hpcwork/jara0052/sichen/AntiprotonHighEnergy_v2.0/B1220_pr.pl1phpsa.l1o9.flux.2016000.4_00_7.8_all_Tree_positive_175_211.root");

TTree *ISSdataTree = (TTree*)ISSdata->Get("ExampleAnalysisTree");
TTree *B1220protonTree = (TTree*)B1220proton->Get("ExampleAnalysisTree");
TH1F* h1 = new TH1F("h1", "Without Weight; Riigidity; N", 100, 170, 220);
TH1F* h2 = new TH1F("h2", "With Weight;    Riigidity; N", 100, 170, 220);
TH1F* h3 = new TH1F("h3", "With Weight3;   Riigidity; N", 100, 170, 220);
TH1F* h4 = new TH1F("h4", "With Weight4;   Riigidity; N", 100, 170, 220);
TH1F* h_ISS = new TH1F("h_ISS",       ";   Riigidity; N", 100, 170, 220);

//// Fill MC
Float_t v_Rigidity;
Double_t v_Weight;

B1220protonTree->SetBranchAddress("Rigidity", &v_Rigidity);
B1220protonTree->SetBranchAddress("Weight", &v_Weight);

int entries = B1220protonTree->GetEntries();
double totalWeight=0;
for( int i=0;i<entries;i++){
    B1220protonTree->GetEntry(i);
    h1->Fill(v_Rigidity);
    h2->Fill(v_Rigidity, v_Weight);
    h3->Fill(v_Rigidity*v_Weight);
    h4->Fill(v_Rigidity, 1.0);
    totalWeight = totalWeight + v_Weight;
}

/*
double totalcontent1=0, totalcontent2=0, totalcontent3=0;
for( int i=0;i<entries;i++){
    totalcontent1 = totalcontent1 + h1->GetBinContent(i);
    totalcontent2 = totalcontent2 + h2->GetBinContent(i);
    totalcontent3 = totalcontent3 + h3->GetBinContent(i);
}

cout<< totalcontent1 <<endl;
cout<< totalcontent2 <<endl;
cout<< totalcontent3 <<endl;
cout<< totalWeight <<endl;
*/

/*
for( int i=0;i<entries;i++){
    h1->SetBinContent(i, h1->GetBinContent(i));
    h2->SetBinContent(i, h2->GetBinContent(i)/totalWeight*entries);
    h3->SetBinContent(i, h3->GetBinContent(i)/totalWeight*entries);
}
*/

//// Fill ISS data
Float_t v_Rigidity_ISS;

ISSdataTree->SetBranchAddress("Rigidity", &v_Rigidity_ISS);
int entries_ISS = ISSdataTree->GetEntries();
for( int i=0;i<entries_ISS;i++){
    ISSdataTree->GetEntry(i);
    h_ISS->Fill(v_Rigidity_ISS);
}



////
double scale_h1 = 1/h1->Integral();
h1->Scale(scale_h1);
double scale_h2 = 1/h2->Integral();
h2->Scale(scale_h2);
double scale_h3 = 1/h3->Integral();
h3->Scale(scale_h3);
double scale_h4 = 1/h4->Integral();
h4->Scale(scale_h4);
double scale_h_ISS = 1/h_ISS->Integral();
h_ISS->Scale(scale_h_ISS);



/*
TFile *Fill_test = new TFile("Fill_test.root","RECREATE");
h1->Write("h1");
h2->Write("h2");
h3->Write("h3");
Fill_test->Close();
*/

TCanvas* c = new TCanvas();
h1->SetMarkerColor(2);
h1->SetLineColor(2);
h2->SetMarkerColor(3);
h2->SetLineColor(3);
h3->SetMarkerColor(4);
h3->SetLineColor(4);
h4->SetMarkerColor(6);
h4->SetLineColor(6);
h_ISS->SetMarkerColor(1);
h_ISS->SetLineColor(1);

h3->SetMarkerColor(4);
h3->SetMarkerStyle(16);
h3->SetMarkerSize(0.8);
h4->SetMarkerColor(6);
h4->SetMarkerStyle(16);
h4->SetMarkerSize(0.8);

h1->Draw("HIST");
h2->Draw("HIST same");
h3->Draw("HIST same");
h4->Draw("P same");
h_ISS->Draw("HIST same");

TLegend *leg = new TLegend(0.25,0.28,0.45,0.48);
leg->AddEntry(h1,"No Reweighting","lpf");
leg->AddEntry(h2,"Reweighting","lpf");
leg->AddEntry(h3,"Wrong Reweighting","lpf");
leg->AddEntry(h4,"Reweighting with 1.0","lpf");
leg->AddEntry(h_ISS,"ISS data","lpf");
leg->SetTextFont(62);
leg->SetTextSize(0.027);
leg->Draw();

c->Print("Fill_Reweighting.pdf");

}
