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

/**
  * \brief A gaussian filter on a 2D histogram. 
  * \param h The histogram to be smeared. It is modified, so copy first the original one if you do not want to lose it!
  * \param sX the sigma of the x-Gaussian
  * \param sY the sigma of the y-Gaussian
  * \note x and y are treated independently
  */

void Gaus2DFilter_test()
{
   chdir("/home/bo791269");
   TFile file("histo_93.0_108.0cccut0.20CCN20TRDN12.root");
   TH2D * h = (TH2D*)file.Get("template_electron");
   const Double_t sX = 20;
   const Double_t sY = 12;

    if(!h) return;
    
    TH2F* h2 = (TH2F*)h->Clone();
    h->Reset();
    double content;
    double x, y, w, z;
    const UInt_t nX = h->GetNbinsX();
    const UInt_t nY = h->GetNbinsY();
    UInt_t kmin, kmax, lmin, lmax;

    for(UInt_t i=1;i<=nX;++i)
    {
        x = h->GetXaxis()->GetBinCenter(i);
        for(UInt_t j=1;j<=nY;++j)
        {
            y = h->GetYaxis()->GetBinCenter(j);
            content = 0;
            //For bin i,j, we look in the old histogram at all bins in a nearby region and calculate the probability of being smeared in the current one
            kmin = h2->GetXaxis()->FindBin(x-3*sX);//For being faster, we reduce the range to 3-sigma. Use kmin = 1 if you want the whole range
            kmax = h2->GetXaxis()->FindBin(x+3*sX);//For being faster, we reduce the range to 3-sigma. Use kmax = nX if you want the whole range
            
            for(UInt_t k=kmin;k<=kmax;++k)
            {
                w = h2->GetXaxis()->GetBinCenter(k);
                lmin = h2->GetYaxis()->FindBin(y-3*sY);//For being faster, we reduce the range to 3-sigma. Use lmin = 1 if you want the whole range
                lmax = h2->GetYaxis()->FindBin(y+3*sY);//For being faster, we reduce the range to 3-sigma. Use lmax = nY if you want the whole range
                
                for(UInt_t l=lmin;l<=lmax;++l)
                {
                    z = h2->GetYaxis()->GetBinCenter(l);
                    //We approximate the integral of the gaussian function between binLowEdge and binUpEdge by the integrand times the bin width. This bin width and the normalization is divided after the end of the loop. The Error Function could be used instead.
                    content+=h2->GetBinContent(k,l)*TMath::Gaus(x,w,sX,false)*TMath::Gaus(y,z,sY,false);
                }
            }
            content=content*h->GetXaxis()->GetBinWidth(i)*h->GetYaxis()->GetBinWidth(j)/2./TMath::Pi()/sX/sY;
            h->SetBinContent(i,j,content);
        }
    }   

   TCanvas *c = new TCanvas;
   h->Draw("");
   gPad->Print("Gaus2DFilter_test.png");

}



