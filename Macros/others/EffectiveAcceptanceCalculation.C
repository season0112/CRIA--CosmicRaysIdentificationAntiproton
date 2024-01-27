
void EffectiveAcceptanceCalculation(){
TFile file("Auxiliary.root");

Cuts::Selector *QualityCuts = (Cuts::Selector*)file.Get("QualityCuts");
Cuts::Selector *McPreselection = (Cuts::Selector*)file.Get("McPreselection");

TH1D *generated = (TH1D*)file.Get("McGeneratedSpectrum/GeneratedEvents");


TH2F *pass = static_cast<const Cuts::EfficiencyHistograms*>(QualityCuts->GetCut(2)->FindAttachment("Cuts::EfficiencyHistograms"))->PassedHistogramMc();
TH2F *total = static_cast<const Cuts::EfficiencyHistograms*>(McPreselection->GetCut(0)->FindAttachment("Cuts::EfficiencyHistograms"))->TotalHistogramMc();

TH1D *passx = pass->ProjectionX();
TH1D *totalx = total->ProjectionX();

TH1D ratio2 = *passx / *generated;
TEfficiency *pEff = new TEfficiency(*passx,*generated);

TH1D acceptance = 3.9 * 3.9 * 3.1415926 * 10000 * ratio2;
TGraphAsymmErrors *graph = pEff->CreateGraph();
Double_t a,b;
for (int i=0; i<100;i++) {
graph->GetPoint(i,a,b);
graph->SetPoint(i, a, b * 3.9 * 3.9 * 3.1415926 * 10000);
graph->SetPointError(i, graph->GetErrorXlow(i), graph->GetErrorXhigh(i), 3.9 * 3.9 * 3.1415926 * 10000 *graph->GetErrorYlow(i), 3.9 * 3.9 * 3.1415926 * 10000 *graph->GetErrorYhigh(i));
}



TCanvas *c = new TCanvas;
acceptance.Draw("hist");
acceptance.GetXaxis()->SetRangeUser(16.6,525);
acceptance.GetXaxis()->SetMoreLogLabels();
gPad->SetLogx(1);
acceptance.SetTitle("");
gStyle->SetOptStat("00000000");
gPad->Print("passtotalratio.png");

TFile *f1 = new TFile("effective_acceptance.root","RECREATE");
acceptance.Write();
pEff->Write();
graph->Write("graph");
f1->Close();


}
