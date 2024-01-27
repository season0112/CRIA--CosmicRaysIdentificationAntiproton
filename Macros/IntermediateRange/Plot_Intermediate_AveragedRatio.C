
int Plot_Intermediate_AveragedRatio(){

//// Load data
// Load Low range
TFile *f1 = new TFile( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_results2016paper.root")).c_str() );
TGraphErrors *unfolded_2016 = (TGraphErrors*)f1->Get("gRatio_unfolded");
TGraphErrors *raw_2016 = (TGraphErrors*)f1->Get("gRatio_Raw");
TGraphErrors *gPublishedRatio = (TGraphErrors*)f1->Get("gPublishedRatio");
TGraphErrors *g_chi2_tof_2016 = (TGraphErrors*)f1->Get("g_chi2_tof");

TFile *f2 = new TFile( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultspass7.8.root")).c_str() );
TGraphErrors *unfolded_pass78 = (TGraphErrors*)f2->Get("gRatio_unfolded");
TGraphErrors *raw_pass78 = (TGraphErrors*)f2->Get("gRatio_Raw");
TGraphErrors *g_chi2_tof_pass78 = (TGraphErrors*)f1->Get("g_chi2_tof");

// Load Intermediate range
TFile *IntermediateFile2016 = new TFile( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_2016paper.root")).c_str() ) ;
TGraphErrors *g_IntermediateResult_2016 = (TGraphErrors*)IntermediateFile2016->Get("gRatio_unfolded");
TGraph *g_chi2dof_intermediate_2016 = (TGraphErrors*)IntermediateFile2016->Get("g_chi2dof_intermediate");
IntermediateFile2016->Close();

TFile *IntermediateFilepass78 = new TFile( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_pass7.8.root")).c_str() );
TGraphErrors *g_IntermediateResult_pass78 = (TGraphErrors*)IntermediateFilepass78->Get("gRatio_unfolded");
TGraph *g_chi2dof_intermediate_pass78 = (TGraphErrors*)IntermediateFilepass78->Get("g_chi2dof_intermediate");
IntermediateFilepass78->Close();

//// Plot
// Plot pass78 ratio 
TCanvas c1("c1","c1",1000,500);
unfolded_pass78->SetMarkerStyle(15);
unfolded_pass78->SetMarkerColor(6);
gPublishedRatio->SetMarkerStyle(15);
gPublishedRatio->SetMarkerColor(1);
g_IntermediateResult_pass78->SetMarkerStyle(15);
g_IntermediateResult_pass78->SetMarkerColor(4);

gPublishedRatio->Draw("AP");
TAxis * xaxis = gPublishedRatio->GetXaxis();
TAxis * yaxis = gPublishedRatio->GetYaxis();
gPublishedRatio->SetTitle("");
xaxis->SetTitle("Rigidity (GV)");
yaxis->SetTitle("Antiproton to Proton ratio");
xaxis->SetLimits(1.0,600);
yaxis->SetRangeUser(0.000000001,0.00028);
xaxis->SetTitleFont(62);
yaxis->SetTitleFont(62);
xaxis->SetTitleSize(0.035);
yaxis->SetTitleSize(0.04);
xaxis->SetLabelFont(62);
xaxis->SetLabelSize(0.05);
yaxis->SetLabelFont(62);
yaxis->SetLabelSize(0.05);

unfolded_pass78->Draw("same P");
g_IntermediateResult_pass78->Draw("same P");

gPad->SetLogx();
xaxis->SetMoreLogLabels();

TLegend *legend1 = new TLegend(0.45,0.15,0.78,0.35);
legend1->AddEntry(unfolded_pass78,"LowRange","lpf");
legend1->AddEntry(g_IntermediateResult_pass78,"IntermediateRange","lpf");
legend1->AddEntry(gPublishedRatio,"PRL paper 2016","lpf");
legend1->SetTextSize(0.05);
legend1->SetTextFont(62);
legend1->Draw();

c1.Update();
c1.SaveAs( "Intermediate_ratio.pdf" );




// Plot chi2
TCanvas c_chi2("c_chi2","c_chi2",1000,500);
TPad *pad_chi2pass78 = new TPad("pad_chi2pass78", "", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
pad_chi2pass78->Draw();
pad_chi2pass78->cd();
pad_chi2pass78->SetBottomMargin(0.2);
pad_chi2pass78->SetLeftMargin(0.2);

g_chi2dof_intermediate_2016->Draw("AP");

g_chi2dof_intermediate_2016->SetMarkerStyle(15);
g_chi2dof_intermediate_2016->SetMarkerColor(1);
TAxis * xaxis_chi2_tof_pass78 = g_chi2dof_intermediate_2016->GetXaxis();
TAxis * yaxis_chi2_tof_pass78 = g_chi2dof_intermediate_2016->GetYaxis();
xaxis_chi2_tof_pass78->SetTitle("Rigidity (GV)");
yaxis_chi2_tof_pass78->SetTitle("Chi2/dof");
xaxis_chi2_tof_pass78->SetTitleFont(62);
yaxis_chi2_tof_pass78->SetTitleFont(62);
xaxis_chi2_tof_pass78->SetTitleSize(0.06);
yaxis_chi2_tof_pass78->SetTitleSize(0.06);
xaxis_chi2_tof_pass78->SetLabelFont(62);
yaxis_chi2_tof_pass78->SetLabelFont(62);
xaxis_chi2_tof_pass78->SetLabelSize(0.05);
yaxis_chi2_tof_pass78->SetLabelSize(0.05);

c_chi2.SaveAs("chi2_tof_pass78.pdf");

return EXIT_SUCCESS;
}

