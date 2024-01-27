
void Plot_below10_ratio(){

//std::string compareversion = "2016paper";
std::string compareversion = "PhysicsReport";

// Load PhysicsReport result
TFile *f_phyreport = new TFile( (string(getenv("MY_ANALYSIS")) + string("/ReferenceFiles/AntiprotonToProtonRatio/AntiprotonToProtonRatio_PhysicsReport/ssdc_canvas.root")).c_str() );
TGraphErrors *gPhyrReortRatio = (TGraphErrors*)f_phyreport->Get("graph1");

//// Load data
// Load Low range
TFile *f1 = new TFile( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_results2016paper.root")).c_str() );
TGraphErrors *g_LowResult_2016          = (TGraphErrors*)f1->Get("gRatio_unfolded");   // Get("gRatio_Raw")
TGraphErrors *raw_g_LowResult_2016Range = (TGraphErrors*)f1->Get("gRatio_Raw");
TGraphErrors *gPublishedRatio           = (TGraphErrors*)f1->Get("gPublishedRatio");
TGraphErrors *g_chi2_tof_2016Range      = (TGraphErrors*)f1->Get("g_chi2_tof");

TFile *f2 = new TFile( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultspass7.8.root")).c_str() );
TGraphErrors *g_LowResult_pass78     = (TGraphErrors*)f2->Get("gRatio_unfolded");
TGraphErrors *raw_g_LowResult_pass78 = (TGraphErrors*)f2->Get("gRatio_Raw");
TGraphErrors *g_chi2_tof_pass78      = (TGraphErrors*)f2->Get("g_chi2_tof");

TFile *f3 = new TFile( (string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultsPhysicsReport.root")).c_str() );
TGraphErrors *g_LowResult_PhysicsReportRange     = (TGraphErrors*)f3->Get("gRatio_unfolded");
TGraphErrors *raw_g_LowResult_PhysicsReportRange = (TGraphErrors*)f3->Get("gRatio_Raw");
TGraphErrors *g_chi2_tof_PhysicsReportRange      = (TGraphErrors*)f3->Get("g_chi2_tof");

// Load Intermediate range
TFile *IntermediateFile2016 = new TFile( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_2016paper.root")).c_str() ) ;
TGraphErrors *g_IntermediateResult_2016 = (TGraphErrors*)IntermediateFile2016->Get("gRatio_unfolded");  // Get("gRatio_Raw")
IntermediateFile2016->Close();

TFile *IntermediateFilepass78 = new TFile( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_pass7.8.root")).c_str() );
TGraphErrors *g_IntermediateResult_pass78 = (TGraphErrors*)IntermediateFilepass78->Get("gRatio_unfolded");
IntermediateFilepass78->Close();

TFile *IntermediateFilePhysicsReport = new TFile( (string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_PhysicsReport.root")).c_str() );
TGraphErrors *g_IntermediateResult_PhysicsReport = (TGraphErrors*)IntermediateFilePhysicsReport->Get("gRatio_unfolded");
IntermediateFilePhysicsReport->Close();


//General define the overlap range.
int LowRangeRemoveAtEnd            = 3;
int IntermediateRangeRemoveAtBegin = 4;
int IntermediateRangeRemoveAtEnd   = 2;
int HighRangeRemovedAtBegin        = 1;

// plot ratio (overlaped range NOT included)
for (int i=0; i<LowRangeRemoveAtEnd; i++){
    g_LowResult_2016->RemovePoint(g_LowResult_2016->GetN()-1);
    g_LowResult_pass78->RemovePoint(g_LowResult_pass78->GetN()-1);
    g_LowResult_PhysicsReportRange->RemovePoint(g_LowResult_PhysicsReportRange->GetN()-1);
}
for (int i=0; i<IntermediateRangeRemoveAtBegin; i++){
    g_IntermediateResult_pass78->RemovePoint(0);
    g_IntermediateResult_2016->RemovePoint(0);
    g_IntermediateResult_PhysicsReport->RemovePoint(0);
}


//////////////
//// Plot ////
////////////// 

TGraphErrors *gReferenceRatio;
TGraphErrors *gLowResult;
TGraphErrors *gIntermediateResult;
if (compareversion == "2016paper"){
    gReferenceRatio     = gPublishedRatio;
    gLowResult          = g_LowResult_2016;
    gIntermediateResult = g_IntermediateResult_2016;
}
else if (compareversion == "PhysicsReport"){
    gReferenceRatio     = gPhyrReortRatio;
    gLowResult          = g_LowResult_PhysicsReportRange;
    gIntermediateResult = g_IntermediateResult_PhysicsReport;
}

TCanvas c1("c1","c1",1000,500);
TPad *pad_ratio = new TPad("pad_ratio", "The pad 80% of the height", 0.0, 0.245, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
TPad *pad_residual = new TPad("pad_residual", "The pad 20% of the height", 0.0, 0.0, 1.0, 0.3, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
pad_ratio->Draw();
pad_residual->Draw();

// upper plot: ratio plot
pad_ratio->cd();
gStyle->SetPadBorderMode(0);
gStyle->SetFrameBorderMode(0);

g_LowResult_pass78->SetMarkerStyle(15);
g_LowResult_pass78->SetMarkerColor(2);
g_LowResult_pass78->SetLineColor(2);
g_IntermediateResult_pass78->SetMarkerStyle(15);
g_IntermediateResult_pass78->SetMarkerColor(2);
gLowResult->SetMarkerStyle(15);
gLowResult->SetMarkerColor(4);
gIntermediateResult->SetMarkerStyle(15);
gIntermediateResult->SetMarkerColor(4);
gReferenceRatio->SetMarkerStyle(21);
gReferenceRatio->SetMarkerColor(1);
gReferenceRatio->SetLineColor(1);

g_LowResult_pass78->Draw("AP");
gReferenceRatio->Draw("same P");
g_IntermediateResult_pass78->Draw("same P");
gLowResult->Draw("same P");
gIntermediateResult->Draw("same P");
g_LowResult_pass78->Draw("same P");

TAxis * xaxis_ratio = g_LowResult_pass78->GetXaxis();
TAxis * yaxis_ratio = g_LowResult_pass78->GetYaxis();
xaxis_ratio->SetTitle("Rigidity (GV)");
yaxis_ratio->SetTitle("Antiproton to Proton ratio");
xaxis_ratio->SetLimits(0.5,10);
yaxis_ratio->SetRangeUser(0,0.0002);
xaxis_ratio->SetTitleFont(62);
yaxis_ratio->SetTitleFont(62);
xaxis_ratio->SetTitleSize(0.05);
yaxis_ratio->SetTitleSize(0.05);
xaxis_ratio->SetLabelFont(62);
yaxis_ratio->SetLabelFont(62);
xaxis_ratio->SetLabelSize(0.07);
yaxis_ratio->SetLabelSize(0.07);
xaxis_ratio->SetTitleOffset(1.0);
yaxis_ratio->SetTitleOffset(0.8);

xaxis_ratio->SetLabelOffset(999); // remove x axis
xaxis_ratio->SetLabelSize(0); // remove x label

TLegend *legend1 = new TLegend(0.48,0.2,0.9,0.4);
legend1->AddEntry(g_LowResult_pass78,"This analysis (All time range)","lpf");
if (compareversion == "2016paper"){
    legend1->AddEntry(gReferenceRatio,"PRL paper 2016","lpf");
    legend1->AddEntry(gLowResult,"This analysis (PRL paper 2016 range)","lpf");
}
if (compareversion == "PhysicsReport"){
    legend1->AddEntry(gReferenceRatio,"PhysicsReport 2021","lpf");
    legend1->AddEntry(gLowResult,"This analysis (PhysicsReport 2021 range)","lpf");
}

legend1->SetTextSize(0.05);
legend1->SetTextFont(62);
legend1->Draw();

// lower plot
pad_residual->cd();
pad_residual->SetBottomMargin(0.3);
pad_residual->SetLeftMargin(0.1);

TGraphErrors *Residual_Low_pass78          = new TGraphErrors(g_LowResult_pass78->GetN());
TGraphErrors *Residual_Intermediate_pass78 = new TGraphErrors(g_IntermediateResult_pass78->GetN());
TGraphErrors *Residual_Low_2016            = new TGraphErrors(gLowResult->GetN());
TGraphErrors *Residual_Intermediate_2016   = new TGraphErrors(gIntermediateResult->GetN());


for (int e = 0; e < g_LowResult_pass78->GetN(); e++){
    Residual_Low_pass78->SetPoint(e, g_LowResult_pass78->GetX()[e], (g_LowResult_pass78->GetY()[e] - gReferenceRatio->GetY()[e])/gReferenceRatio->GetY()[e]);
    Residual_Low_2016  ->SetPoint(e, gLowResult->GetX()[e], (gLowResult->GetY()[e] - gReferenceRatio->GetY()[e])/gReferenceRatio->GetY()[e]); // gLowResult and pass78 have same N.
}

int IntermediateStartIndex = g_LowResult_pass78->GetN();
for (int e = 0; e < g_IntermediateResult_pass78->GetN(); e++){
    Residual_Intermediate_pass78->SetPoint(e, g_IntermediateResult_pass78->GetX()[e], (g_IntermediateResult_pass78->GetY()[e] - gReferenceRatio->GetY()[e+IntermediateStartIndex]) / gReferenceRatio->GetY()[e+IntermediateStartIndex]);
    Residual_Intermediate_2016  ->SetPoint(e, gIntermediateResult->GetX()[e], (gIntermediateResult->GetY()[e] - gReferenceRatio->GetY()[e+IntermediateStartIndex]) / gReferenceRatio->GetY()[e+IntermediateStartIndex]);
}

Residual_Low_pass78->Draw("AP");

Residual_Low_pass78->SetTitle("");
Residual_Low_pass78->SetMarkerStyle(15);
Residual_Low_pass78->SetMarkerColor(2);
Residual_Intermediate_pass78->SetMarkerStyle(15);
Residual_Intermediate_pass78->SetMarkerColor(2);
Residual_Low_2016->SetMarkerStyle(15);
Residual_Low_2016->SetMarkerColor(4);
Residual_Intermediate_2016->SetMarkerStyle(15);
Residual_Intermediate_2016->SetMarkerColor(4);

TAxis * xaxis_residuallow = Residual_Low_pass78->GetXaxis();
TAxis * yaxis_residuallow = Residual_Low_pass78->GetYaxis();
xaxis_residuallow->SetLimits(0.5,10);
yaxis_residuallow->SetRangeUser(-0.25,0.35);
xaxis_residuallow->SetTitle("Rigidity (GV)");
if (compareversion == "2016paper"){
    yaxis_residuallow->SetTitle("#frac{This - PRL}{PRL}");
}
else if (compareversion == "PhysicsReport"){
    yaxis_residuallow->SetTitle("#frac{This - PhyRep}{PhyRep}");
}
xaxis_residuallow->SetMoreLogLabels();
xaxis_residuallow->SetTitleFont(62);
yaxis_residuallow->SetTitleFont(62);
xaxis_residuallow->SetTitleSize(0.15);
yaxis_residuallow->SetTitleSize(0.1);
xaxis_residuallow->SetLabelFont(62);
yaxis_residuallow->SetLabelFont(62);
xaxis_residuallow->SetLabelSize(0.15);
yaxis_residuallow->SetLabelSize(0.1);
yaxis_residuallow->SetTitleOffset(0.4);
xaxis_residuallow->SetTickSize(0.1);
xaxis_residuallow->SetTitleOffset(0.7);

Residual_Intermediate_pass78->Draw("same P");
Residual_Low_2016->Draw("same P");
Residual_Intermediate_2016->Draw("same P");

TLine line(0.5, 0, 10, 0);
line.Draw("same");
TLine line2(0.5, -0.1, 10, -0.1);
line2.SetLineStyle(2);
line2.Draw("same");
TLine line3(0.5, 0.1, 10, 0.1);
line3.SetLineStyle(2);
line3.Draw("same");

if (compareversion == "2016paper"){
    c1.SaveAs("below10ratio_NOT_overlapped_2016paper.pdf");
}
else if (compareversion == "PhysicsReport"){
    c1.SaveAs("below10ratio_NOT_overlapped_PhysicsReport.pdf");
}

}
