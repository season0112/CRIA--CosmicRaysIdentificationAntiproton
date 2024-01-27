// Example Usage: root -L /home/bo791269/Software/AntiprotonAnalysis/Libraries/AntiprotonBinning/AntiprotonBinning.C -x Plot_CheckTwoVersionResult.C
#include "Utilities.hh"
#include "Plot_FullRangeRatio.hh"

void Plot_CheckTwoVersionResult(){

TFile *CheckResult_new         = new TFile("CheckResult_new.root");
TFile *CheckResult_lastversion = new TFile("CheckResult_lastversion.root");
TFile *CheckResult_old         = new TFile("CheckResult_old.root");

TGraph *g_Antiproton_number_unfolded_Low_new          = (TGraph*) CheckResult_new->Get("g_Antiproton_number_unfolded_Low");
TGraph *g_Antiproton_number_unfolded_Intermediate_new = (TGraph*) CheckResult_new->Get("g_Antiproton_number_unfolded_Intermediate");
TGraph *g_Antiproton_number_unfolded_High_new         = (TGraph*) CheckResult_new->Get("g_Antiproton_number_unfolded_High");
TGraphErrors *g_LowResult_new                         = (TGraphErrors*) CheckResult_new->Get("g_LowResul");
TGraphErrors *g_Intermediate_new                      = (TGraphErrors*) CheckResult_new->Get("g_Intermediate");
TGraphErrors *g_HighResult_new                        = (TGraphErrors*) CheckResult_new->Get("g_HighResult");

TGraph *g_Antiproton_number_unfolded_Low_lastversion          = (TGraph*) CheckResult_lastversion->Get("g_Antiproton_number_unfolded_Low");
TGraph *g_Antiproton_number_unfolded_Intermediate_lastversion = (TGraph*) CheckResult_lastversion->Get("g_Antiproton_number_unfolded_Intermediate");
TGraph *g_Antiproton_number_unfolded_High_lastversion         = (TGraph*) CheckResult_lastversion->Get("g_Antiproton_number_unfolded_High");

TGraphErrors *g_LowResult_old                         = (TGraphErrors*) CheckResult_old->Get("g_LowResul");
TGraphErrors *g_Intermediate_old                      = (TGraphErrors*) CheckResult_old->Get("g_Intermediate");
TGraphErrors *g_HighResult_old                        = (TGraphErrors*) CheckResult_old->Get("g_HighResult");




// Plot numbers
TCanvas c_number("c_number","c_number",1000,500);

g_Antiproton_number_unfolded_Low_new->Draw("AP");
g_Antiproton_number_unfolded_Intermediate_new->Draw("same P");
//g_Antiproton_number_unfolded_High_new->Draw("same P");
g_Antiproton_number_unfolded_Low_lastversion->Draw("same P");
g_Antiproton_number_unfolded_Intermediate_lastversion->Draw("same P");
//g_Antiproton_number_unfolded_High_lastversion->Draw("same P");

g_Antiproton_number_unfolded_Low_new->SetMarkerStyle(15);
g_Antiproton_number_unfolded_Low_new->SetMarkerColor(2);
g_Antiproton_number_unfolded_Intermediate_new->SetMarkerColor(2);
g_Antiproton_number_unfolded_High_new->SetMarkerColor(2);

g_Antiproton_number_unfolded_Low_lastversion->SetMarkerStyle(15);
g_Antiproton_number_unfolded_Low_lastversion->SetMarkerColor(4);
g_Antiproton_number_unfolded_Intermediate_lastversion->SetMarkerColor(4);
g_Antiproton_number_unfolded_High_lastversion->SetMarkerColor(4);


TAxis * xaxis_number = g_Antiproton_number_unfolded_Low_new->GetXaxis();
TAxis * yaxis_number = g_Antiproton_number_unfolded_Low_new->GetYaxis();
xaxis_number->SetLimits(1, 20);
yaxis_number->SetRangeUser(1, 90000);
xaxis_number->SetTitle("Rigidity (GV)");
yaxis_number->SetTitle("Antiproton Numbers");
xaxis_number->SetTitleFont(62);
yaxis_number->SetTitleFont(62);
xaxis_number->SetTitleSize(0.045);
yaxis_number->SetTitleSize(0.045);
xaxis_number->SetLabelFont(62);
xaxis_number->SetLabelSize(0.05);
yaxis_number->SetLabelFont(62);
yaxis_number->SetLabelSize(0.05);

TLegend *legend_number = new TLegend(0.60, 0.15, 0.80, 0.42);
legend_number->AddEntry(g_Antiproton_number_unfolded_Low_lastversion, "before", "p");
legend_number->AddEntry(g_Antiproton_number_unfolded_Low_new, "after" , "p");
legend_number->SetTextSize(0.04);
legend_number->SetTextFont(62);
legend_number->Draw();

g_Antiproton_number_unfolded_Low_new->SetTitle("");
gPad->SetLogy();
gPad->SetLogx();

c_number.SaveAs("number_compare.pdf");


// Plot ratio
int lownumber = g_LowResult_old->GetN();
int internumber = g_Intermediate_old->GetN();
int highnumber = g_HighResult_old->GetN();
cout<< "lownumber: " << lownumber <<endl;
cout<< "internumber: " << internumber <<endl;
cout<< "highnumber: " << highnumber <<endl;

TGraphErrors *g_oldraio = new TGraphErrors( lownumber+internumber+highnumber );

for (int i=0; i<lownumber; i++){
    g_oldraio->SetPoint(i, g_LowResult_old->GetX()[i], g_LowResult_old->GetY()[i] ); 
    cout<< "low:" << g_LowResult_old->GetX()[i] << endl;
}

for (int i=0; i<internumber; i++){
    g_oldraio->SetPoint(lownumber+i, g_Intermediate_old->GetX()[i], g_Intermediate_old->GetY()[i] );
    cout<< "Intermediate:" << g_Intermediate_old->GetX()[i] << endl;
}

for (int i=0; i<highnumber; i++){
    g_oldraio->SetPoint(lownumber+internumber+i, g_HighResult_old->GetX()[i], g_HighResult_old->GetY()[i] );
    cout<< "High:" << g_HighResult_old->GetX()[i] << endl;
}

std::string issversion = "PhyRep2021";
Plot_Ratio_CompareWithReference_NotOverlapped(g_oldraio, g_oldraio, g_HighResult_new, g_LowResult_new, g_Intermediate_new, issversion); // (Overlapped range NOT Included)


TCanvas c_plot("c_plot","c_plot",1000,500);

g_LowResult_old->Draw("");
g_Intermediate_old->Draw("same");
g_HighResult_old->Draw("same");
g_oldraio->Draw("same");

TAxis * xaxis = g_LowResult_old->GetXaxis();
TAxis * yaxis = g_LowResult_old->GetYaxis();
xaxis->SetLimits(1, 600);
yaxis->SetRangeUser(0.000000001,0.00028);
xaxis->SetTitle("Rigidity (GV)");
yaxis->SetTitle("Antiproton Numbers Ratio (This analysis / PhyRep)");
xaxis->SetTitleFont(62);
yaxis->SetTitleFont(62);
xaxis->SetTitleSize(0.045);
yaxis->SetTitleSize(0.045);
xaxis->SetLabelFont(62);
xaxis->SetLabelSize(0.05);
yaxis->SetLabelFont(62);
yaxis->SetLabelSize(0.05);

gPad->SetLogx();

c_plot.SaveAs("test.pdf");




}






