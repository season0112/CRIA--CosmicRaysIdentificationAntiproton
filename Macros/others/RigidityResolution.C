#include <iostream>
#include <fstream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"

void RigidityResolution(){

  Double_t x[34] = {162. ,  192. ,  222. ,  252. ,  282. ,  312. ,  342. ,  372. ,
        402. ,  432. ,  462. ,  492. ,  522. ,  552. ,  582. ,  612. ,
        642. ,  672. ,  702. ,  732. ,  762. ,  792. ,  822. ,  852. ,
        882. ,  912. ,  942. ,  972. ,  993.5, 1250. , 1750. , 2250. ,
       2750. , 3250.};
   Double_t yresolution[34];
   Double_t yresolution_ratio[34];
/*
   Double_t x[185] = {152.,  162.,  172.,  182.,  192.,  202.,  212.,  222.,  232.,
        242.,  252.,  262.,  272.,  282.,  292.,  302.,  312.,  322.,
        332.,  342.,  352.,  362.,  372.,  382.,  392.,  402.,  412.,
        422.,  432.,  442.,  452.,  462.,  472.,  482.,  492.,  502.,
        512.,  522.,  532.,  542.,  552.,  562.,  572.,  582.,  592.,
        602.,  612.,  622.,  632.,  642.,  652.,  662.,  672.,  682.,
        692.,  702.,  712.,  722.,  732.,  742.,  752.,  762.,  772.,
        782.,  792.,  802.,  812.,  822.,  832.,  842.,  852.,  862.,
        872.,  882.,  892.,  902.,  912.,  922.,  932.,  942.,  952.,
        962.,  972.,  982.,  992., 1002., 1012., 1022., 1032., 1042.,
       1052., 1062., 1072., 1082., 1092., 1102., 1112., 1122., 1132.,
       1142., 1152., 1162., 1172., 1182., 1192., 1202., 1212., 1222.,
       1232., 1242., 1252., 1262., 1272., 1282., 1292., 1302., 1312.,
       1322., 1332., 1342., 1352., 1362., 1372., 1382., 1392., 1402.,
       1412., 1422., 1432., 1442., 1452., 1462., 1472., 1482., 1492.,
       1502., 1512., 1522., 1532., 1542., 1552., 1562., 1572., 1582.,
       1592., 1602., 1612., 1622., 1632., 1642., 1652., 1662., 1672.,
       1682., 1692., 1702., 1712., 1722., 1732., 1742., 1752., 1762.,
       1772., 1782., 1792., 1802., 1812., 1822., 1832., 1842., 1852.,
       1862., 1872., 1882., 1892., 1902., 1912., 1922., 1932., 1942.,
       1952., 1962., 1972., 1982., 1992.};
   Double_t yresolution[34];
*/
   chdir("/p/scratch/cvsk10/li8/analysis_v7.0/");
   TFile file("RigidityResolution_new.root");
   TH2D *RigidityResolution = (TH2D*)file.Get("RigidityResolutioni_Choutko");

   for (int i =1;i<35;i=i+1){
   yresolution[i-1] = RigidityResolution->ProjectionY("",i,i)->GetStdDev(); // true x, reconstructed y.
   yresolution_ratio[i-1] = RigidityResolution->ProjectionY("",i,i)->GetStdDev()/x[i-1];  // true x, reconstructed y.
   cout<< x[i-1] <<endl;
   cout<< yresolution[i-1] <<endl;
   cout<< i << endl;
   }

/*
   Double_t x2[n-1];
   strncpy(x2,x,n-1,1);  
   Double_t y2[n-1];
   strncpy(y2,yresolution,n-1,1);
*/

   TCanvas *c = new TCanvas;
   TGraph *gr1 = new TGraph (34, x, yresolution);
   gr1->Draw("AC*");
   gPad->SetLogx(1);
   gr1->GetXaxis()->SetRangeUser(10, 1000);
   gPad->Print("RigidityResolution.png");

/*
 * TImage *img = TImage::Create();
   img->FromPad(c);
   img->WriteImage("RigidityResolution.png");
*/
   TGraph *gr2 = new TGraph (34, x, yresolution_ratio);
   gr2->Draw("AC*");
   gPad->SetLogx(1);
   gr2->GetXaxis()->SetRangeUser(10, 1000);
   gPad->Print("RigidityResolution_2.png");

}



