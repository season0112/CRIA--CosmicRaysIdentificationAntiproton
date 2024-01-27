#include "slices.C"


void makeTimeCalibrations(int sector) {
 gROOT->Reset(); 


gStyle->SetCanvasColor(kWhite);
gStyle->SetPadGridY(kTRUE);
gStyle->SetPadGridX(kTRUE);

//gStyle->SetOptStat(0);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(1);
  gStyle->SetStatX(0.87);
  gStyle->SetStatH(0.08);
   gStyle->SetStatY(0.95);
 // gStyle->SetStatW(0.19);

gStyle->SetLabelSize(0.05);

gStyle->SetNumberContours(511);
gStyle->SetPalette(1);
 gStyle->SetErrorX(0.2);


 int ips = sector;
 if (ips<1 || ips>8) {
   cout << " wrong sector number" << endl;
   return;
 }

 TCut cleanti = "acchi2<30 && npairs>0 && adc0>0 && adc1>0";
 TCut trdacc = "abs(phimiss-phibar)<15 && abs(AntiCrossZ)<35";
 TCut zacc = "abs(AntiCrossZ-cal_zeta)<15";


 TCut evclean = trdacc && cleanti && zacc; 
 TCanvas* ctime = new TCanvas("ctime"," ",1200,800);


 double t0 = 1310.e6;
 int ny = 8;


   TChain* tree = (TChain*) new TChain("tree");
   tree->Add(Form("../REDUCEDNTUPLE/results/*.root"));  
   tree->SetAlias("toftime","0.5*(Tusedtof0+Tusedtof1)");
   TH2D* gg1 = new TH2D("gg1","gg1",365*ny,t0,24*3600*365*ny+t0,250,-5,5);
   SetXaxisTimeHist(gg1);
   tree->Draw("uncal_time-toftime-dista/(29.98):run>>gg1",Form("sector==%d",ips) && evclean,"colz");   
   TH1D* gg1x = (TH1D*) gg1->ProjectionX("gg1x");
   FindMaxX(gg1,gg1x,0,0,6,10.);

   TFile mfout(Form("OUTCAL/Timecal%d.root",ips),"recreate");
   gg1x->Write(Form("Timecal%d",ips));
   mfout.Close();
  
   SliceNormalizeX(gg1,-1,(TGraph*)1,(TGraph*)1,0,0);
   gg1->Draw("colz");
   gg1->GetYaxis()->SetTitle("[Time(Acc) - Expected(dist,beta=-1,Time(Tof))] (ns)");
   gg1x->Draw("same");

   ctime->SaveAs(Form("png_cal/Timecal%d.png",ips));
   ctime->SaveAs(Form("png_cal/Timecal%d.root",ips));
    
 return;
}
