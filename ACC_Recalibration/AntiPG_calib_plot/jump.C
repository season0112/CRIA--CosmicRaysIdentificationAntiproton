void jump()
{
   TFile *f_new = new TFile("AntiPG_calib_2.root");

   TH1D *T1_new = (TH1D*)f_new->Get("Timecal1");

//   jump plot begin:1384905600 11/20/2013
//   jump plot end:1386201600, 12/05/2013
//   new histogram begin:1310000000 07/07/2011,
//   old histogram begin:1306000000 05/21/2011
//   old actual data ending:1383264000 11/01/2013, 
//   old histogram endling:1400608000 05/20/2014, 
//   new histogram endling:1562288000 07/05/2019,

   TCanvas * c_t1_new = new TCanvas;
   gStyle->SetOptStat(0);
   T1_new->SetTitle("Time Recalibration");
   T1_new->GetYaxis()->SetRangeUser(-1, 2);
   T1_new->GetYaxis()->SetTitle("[Time(Acc) - Expected(dist,beta=-1,Time(Tof))] (ns)");
   T1_new->GetXaxis()->SetRangeUser(1384905600, 1386201600);
   T1_new->GetXaxis()->SetTimeFormat("%d\/%m\/%Y");
   T1_new->GetXaxis()->SetLabelSize(0.02);

   T1_new->Draw("colz");
   c_t1_new->SaveAs("jump.pdf");

}









