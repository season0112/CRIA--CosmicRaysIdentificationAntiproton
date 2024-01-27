#include "slices.C"


void makeChargeCalibrations(int sector) {
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

 int ips=sector;
 if(ips<1 || ips>8) {
   std::cout << " error sector" << std::endl;
   return;
 }


 TCut cleanti = "acchi2<30 && npairs>0 && adc0>0 && adc1>0";
 TCut trdacc = "abs(phimiss-phibar)<15 && abs(AntiCrossZ)<35";
 TCut tacc = "abs(uncal_time-toftime-dista/(29.98))<2";
 TCut zacc = "abs(uncal_zeta-AntiCrossZ)<15";



 TCut evclean = trdacc && cleanti && tacc && zacc; 


 TCanvas* czeta = new TCanvas("czeta"," ",1200,800);


 double t0 = 1310.e6;
 int ny = 8;


   TChain* tree = (TChain*) new TChain("tree");
   tree->Add(Form("../REDUCEDNTUPLE/OUTROOT_nominal/out*.root"));
   tree->Add(Form("../REDUCEDNTUPLE/OUTROOT_ttcs_off/out*.root"));
   std::cout<< "tree" << tree->GetEntries() <<std::endl;
//   TTree *tree = ch->GetTree();
//   std::cout<< "tree" << tree <<std::endl;
   tree->SetAlias("toftime","0.5*(Tusedtof0+Tusedtof1)");
   tree->SetAlias("qtot","0.49*(0.5*(Qusedtof0+Qusedtof1)+1.15*QTRD)");

   TH2D* gg2_1 = new TH2D("gg2_1","gg2_1",365*ny,t0,24*3600*365*ny+t0,200,0,5.9);
   SetXaxisTimeHist(gg2_1);
   tree->Draw("uncal_qside1:run>>gg2_1",Form("(qtot<2.2 && qtot>1.8) && sector==%d",ips) && evclean,"colz");   
   TH1D* gg2_1x = (TH1D*) gg2_1->ProjectionX("gg2_1x");
   FindMaxX(gg2_1,gg2_1x,0,0,6,10.);
   SliceNormalizeX(gg2_1,-1,(TGraph*)1,(TGraph*)1,0,0);
   gg2_1->Draw("colz");
   gg2_1->GetYaxis()->SetTitle("ACC side 1 charge for He");
   gg2_1x->Draw("same");
   czeta->SaveAs(Form("test/QTimecal2_%d_1.png",ips));
   czeta->SaveAs(Form("test/QTimecal2_%d_1.root",ips));

   TH2D* gg2_0 = new TH2D("gg2_0","gg2_0",365*ny,t0,24*3600*365*ny+t0,200,0,5.9);
   SetXaxisTimeHist(gg2_0);
   tree->Draw("uncal_qside0:run>>gg2_0",Form("(qtot<2.2 && qtot>1.8) && sector==%d",ips) && evclean,"colz");   
   TH1D* gg2_0x = (TH1D*) gg2_0->ProjectionX("gg2_0x");
   FindMaxX(gg2_0,gg2_0x,0,0,6,10.);
   SliceNormalizeX(gg2_0,-1,(TGraph*)1,(TGraph*)1,0,0);
   gg2_0->Draw("colz");
   gg2_0->GetYaxis()->SetTitle("ACC side 0 charge for He");
   gg2_0x->Draw("same");
   czeta->SaveAs(Form("test/QTimecal2_%d_0.png",ips));
   czeta->SaveAs(Form("test/QTimecal2_%d_0.root",ips));

  
 return;
}
