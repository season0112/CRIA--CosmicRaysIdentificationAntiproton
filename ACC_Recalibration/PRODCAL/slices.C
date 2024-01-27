#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TCut.h>
#include <iostream>
#include <TF1.h>

using namespace std;

double Dividemean(TH1D* h){
  // each bin is divided to the mean-Y                                                                                          
  // assume that bin at y=0 are empty (so not in mean)                                                                           
  double mean = 0;
  double pesum = 0;
  for (int ii=0; ii<=h->GetNbinsX(); ii++) {
    if(h->GetBinContent(ii)!=0 && h->GetBinError(ii)!=0){
      double pes = fabs(1./h->GetBinError(ii));
      mean=mean+pes*h->GetBinContent(ii);
      pesum=pesum+pes;
    }
  }
  mean = mean/pesum;
  for (int ii=0; ii<=h->GetNbinsX(); ii++) {
    h->SetBinError(ii,h->GetBinError(ii)/mean);
    if(h->GetBinContent(ii)!=0){
      if (h->GetBinError(ii)==0) {
        h->SetBinContent(ii,0);}
      else {
        h->SetBinContent(ii,h->GetBinContent(ii)/mean);
      }
    }
  }
  return mean;
}

double Removemean(TH1D* h){
  // remove the Y-mean from an histo                                                                                             
  // assume that bin at y=0 are empty (so not in mean)                                                                           
  double mean = 0;
  double pesum = 0;
  for (int ii=0; ii<=h->GetNbinsX(); ii++) {
    if(h->GetBinContent(ii)!=0 && h->GetBinError(ii)!=0){
      double pes = fabs(1./h->GetBinError(ii));
      mean=mean+pes*h->GetBinContent(ii);
      pesum=pesum+pes;
    }
  }
  mean = mean/pesum;
  for (int ii=0; ii<=h->GetNbinsX(); ii++) {
    if(h->GetBinContent(ii)!=0){
      if (h->GetBinError(ii)==0) {
        h->SetBinContent(ii,0);}
      else {
        h->SetBinContent(ii,h->GetBinContent(ii)-mean);}
    }
  }
  return mean;
}

void SetXaxisTimeHist(TH2D *h){
  // nice time display
  h->GetXaxis()->SetTimeDisplay(1);
  h->GetXaxis()->SetNdivisions(-512);
  // gr->GetXaxis()->SetTimeFormat("%d-%b");                                                                                     
  h->GetXaxis()->SetTimeOffset(0,"gmt");
  h->GetXaxis()->SetLabelOffset(0.04);
  h->GetXaxis()->SetTimeFormat("#splitline{%Y}{%b}");
}


int SliceNormalizeY(TH2* h, double tail, TGraph* lower, TGraph* upper, int iofl, double ecur){
  // ecur = N will fill error of empty bins with N/norm                                                                                                          
  // iofl = 1 consider overflow  iofl=0 only what you see                                                                                                        
  // if tails not necessary call is as:                                                                                                                          
  // SliceNormalizeY(histo,-1.,(TGraph*) 1,(TGraph*) 1,0,0);                                                                                                       
  // will return 2                                                                                                                                               
    double integralslice=0;
  if (iofl != 0 && iofl != 1) return 5;  // wronf iofl                                                                                                           
  if (!h) return 1;  // no histogram                                                                                                                             
  for (int ii=0; ii<=h->GetNbinsY()+1; ii++) {
    integralslice=0;
    integralslice=h->Integral(1-iofl, h->GetNbinsX()+iofl, ii, ii);
    if (integralslice>0) {
      for (int jj=1-iofl; jj<=h->GetNbinsX()+iofl; jj++) {
        if (integralslice) {
          h->SetBinContent(jj, ii, h->GetBinContent(jj, ii)/integralslice);
          h->SetBinError(jj, ii, h->GetBinError(jj, ii)/integralslice);
          if ((h->GetBinContent(jj, ii)) == 0) h->SetBinError(jj, ii, ecur/integralslice);
        } } } }
  if (tail>0.5 || tail<= 0) return 2;  // not reasonable tail                                                                                                    
  if(lower->GetN() != h->GetNbinsY()) return 3; // lower tgraph not correct                                                                                      
  if(upper->GetN() != h->GetNbinsY()) return 4; // upper tgraph not correct                                                                                      
  for (int ii=1; ii<=h->GetNbinsY(); ii++) {
    double upp=h->GetXaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    double low=h->GetXaxis()->GetBinLowEdge(1);
    double cent = h->GetYaxis()->GetBinCenter(ii);
    for (int ix=1-iofl; ix<=h->GetNbinsX(); ix++) {
      integralslice=h->Integral(1-iofl, ix, ii, ii);
      if (integralslice<tail) low = h->GetXaxis()->GetBinLowEdge(ix+1);
      if (integralslice<1.-tail) upp = h->GetXaxis()->GetBinLowEdge(ix+1);
    }
    lower->SetPoint(ii-1,low,cent);
    upper->SetPoint(ii-1,upp,cent); }
  return 0;  // regular                                                                                                                                          
}



int SliceNormalizeX(TH2* h, double tail, TGraph* lower, TGraph* upper, int iofl, double ecur){
  // ecur = N will fill error of empty bins with N/norm                                                                          
  // iofl = 1 consider overflow  iofl=0 only what you see      
  // if tails not necessary call is as:                        
  // SliceNormalizeX(histo,-1.,(TGraph*) 1,(TGraph*) 1,0,0);   
  // will return 2              
  double integralslice=0;
  if (iofl != 0 && iofl != 1) return 5;  // wronf iofl         
  if (!h) return 1;  // no histogram                                                             

  for (int ii=0; ii<=h->GetNbinsX()+1; ii++) {
    integralslice=0;
    integralslice=h->Integral( ii, ii, 1-iofl, h->GetNbinsY()+iofl);
    if (integralslice>0) {
      for (int jj=1-iofl; jj<=h->GetNbinsY()+iofl; jj++) {
        if (integralslice) {
          h->SetBinContent(ii, jj, h->GetBinContent(ii, jj)/integralslice);
          h->SetBinError(ii, jj, h->GetBinError(ii, jj)/integralslice);
          if ((h->GetBinContent(ii, jj)) == 0) h->SetBinError(ii, jj, ecur/integralslice);
        } } } }
  if (tail>0.5 || tail<= 0) return 2;  // not reasonable tail                                         
  if(lower->GetN() != h->GetNbinsX()) return 3; // lower tgraph not correct                           
  if(upper->GetN() != h->GetNbinsX()) return 4; // upper tgraph not correct                           
  for (int ii=1; ii<=h->GetNbinsX(); ii++) {
    double upp=h->GetYaxis()->GetBinLowEdge(h->GetNbinsX()+1);
    double low=h->GetYaxis()->GetBinLowEdge(1);
    double cent = h->GetXaxis()->GetBinCenter(ii);
    for (int iy=1-iofl; iy<=h->GetNbinsY(); iy++) {
      integralslice=h->Integral(ii,ii,1-iofl, iy);
      if (integralslice<tail) low = h->GetYaxis()->GetBinLowEdge(iy+1);
      if (integralslice<1.-tail) upp = h->GetYaxis()->GetBinLowEdge(iy+1);
    }
    lower->SetPoint(ii-1,cent,low);
    upper->SetPoint(ii-1,cent,upp); }
  return 0;  // regular                                                                                                          
}


int FindBinOfMax(TH1* h,int iofl){
  // find bin of maximum with request that is not isolated
  if (!h) return -2;  // no histogram                                                           
    int imax=-1;
    int inext=-1;
    double vmax=-1.e-30;
    double vnext=-1.e-30;
  for (int ii=0-iofl; ii<=h->GetNbinsX()+iofl; ii++) {
    //      cout << ii << " " << vmax << " " << h->GetBinContent(ii) << endl;
    if (vmax < h->GetBinContent(ii)){
      vnext=vmax;
      inext=imax;
      vmax = h->GetBinContent(ii);
      imax = ii;}
    else{
    if (vnext < h->GetBinContent(ii)){
      vnext = h->GetBinContent(ii);
      inext = ii;}
    }
  }
  // check sidebands
  double sidemax=0;
  double sidenext=0;
  int ii=(imax-1);
  if (ii>0) sidemax=sidemax+h->GetBinContent(ii);
  ii=(imax+1);
  if (ii<h->GetNbinsX()) sidemax=sidemax+h->GetBinContent(ii);
  ii=(inext-1);
  if (ii>0) sidenext=sidenext+h->GetBinContent(ii);
  ii=(inext+1);
  if (ii<h->GetNbinsX()) sidenext=sidenext+h->GetBinContent(ii);

  if (sidenext>sidemax && abs(inext-imax)>2) imax = inext;
  return imax;}

double lagrint(double x1,double x2,double x3,double y1,double y2,double y3,double x, double& xextr, double& a, double& b, double& c ){ 
// return 3 points quadratic lagrange interpolation
  // y = ax^2+bx+c

  //  cout << xextr << " prima" << endl;
  if (x1==x2 || x1==x3 || x2==x3) return (x2+x1+x3)/3.;
  a = y1/((x1-x2)*(x1-x3))+y2/((x2-x1)*(x2-x3))+y3/((x3-x1)*(x3-x2));
  b = -(x2+x3)*y1/((x1-x2)*(x1-x3))-(x1+x3)*y2/((x2-x1)*(x2-x3))-(x1+x2)*y3/((x3-x1)*(x3-x2));
  c = (x2*x3)*y1/((x1-x2)*(x1-x3))+(x1*x3)*y2/((x2-x1)*(x2-x3))+(x1*x2)*y3/((x3-x1)*(x3-x2));
  if (a==0) {xextr=(x2+x1+x3)/3.;}
  else {xextr=-b/(2.*a);}
  //  cout << xextr << " dopo" << endl;
  return a*x*x+b*x+c;
}


int FindMaxY(TH2* h, TH1* hx, int iofl, int irms = 0, int i3 = 3, double escale = 1){
  // For each Y slice fill the hx with the position of maximum and error is RMS/sqrt(N) of the slice
  // iofl = 1 consider overflow  iofl=0 only what you see                                      
  // if irms = 0 the error is the rms/sqrt(integral)
  // if irms = 1 the error is the rms
  // if irms = 2 the error is the sigma extracted from curvature
  // i3 is the distance of the nearest bin used in lagrange interpolation
  // escale = error multiplicative scaling
  if (i3<1) i3 = 1;
  if (iofl != 0 && iofl != 1) return 5;  // wronf iofl                                         

  if (!h) return 1;  // no histogram                                                           

  if(hx->GetNbinsX() != h->GetNbinsY()) return 3; // lower tgraph not correct                  

  for (int ii=0; ii<=h->GetNbinsY()+1; ii++) {
    double integralslice=h->Integral(1-iofl, h->GetNbinsX()+iofl,ii,ii);
    TH1D* htmp = (TH1D*) h->ProjectionX("htmp",ii,ii);
    int iy = FindBinOfMax(htmp,iofl);
    double maxy = 0;
    double a = 0;
    double ymax = 0;
    double emaxy = fabs((htmp->GetBinCenter(1))-htmp->GetBinCenter(htmp->GetNbinsX()));
    if (iy>0) {maxy=htmp->GetBinCenter(iy);
      if (iy>i3 && iy<htmp->GetNbinsX()+1-i3) {
	double b,c,xx=0;
	lagrint(htmp->GetBinCenter(iy-i3),htmp->GetBinCenter(iy),htmp->GetBinCenter(iy+i3),htmp->GetBinContent(iy-i3),htmp->GetBinContent(iy),htmp->GetBinContent(iy+i3),xx,maxy,a,b,c);
	ymax = -b*b/(4.*a)+c;
      }
      emaxy = htmp->GetRMS();
      if (emaxy == 0 || integralslice <= 2) {maxy =0; emaxy=fabs((htmp->GetBinCenter(1))-htmp->GetBinCenter(htmp->GetNbinsX()));}
      else {
	if (irms==2 && a<0) {emaxy = sqrt(-ymax/(2.*a));}
      if (irms==1) {emaxy = emaxy;}
      if (irms==0) {emaxy = emaxy/sqrt(fabs(integralslice));}

      }}
    //       cout << ii << " " << maxy << " " << emaxy << " " << integralslice << " " << iy << " " << htmp->GetRMS() << endl;
    htmp->Delete();
    hx->SetBinContent(ii,maxy);
    hx->SetBinError(ii,escale*emaxy);
  }
  return 0;
}



int FindMaxX(TH2* h, TH1* hx, int iofl, int irms = 0, int i3 = 3, double escale = 1){
  // For each X slice fill the hx with the position of maximum from fit and error is RMS/sqrt(N) of the slice
  // iofl = 1 consider overflow  iofl=0 only what you see                                      
  // if irms = 0 the error is the rms/sqrt(integral)
  // if irms = 1 the error is the rms
  // if irms = 2 the error is the sigma extracted from curvature
  // i3 is the distance of the nearest bin used in lagrange interpolation
  // escale = error multiplicative scaling
  if (i3<1) i3 = 1;
  if (iofl != 0 && iofl != 1) return 5;  // wronf iofl                                         

  if (!h) return 1;  // no histogram                                                           

  if(hx->GetNbinsX() != h->GetNbinsX()) return 3; // lower tgraph not correct                  

  for (int ii=0; ii<=h->GetNbinsX()+1; ii++) {
    double integralslice=h->Integral( ii, ii, 1-iofl, h->GetNbinsY()+iofl);
    TH1D* htmp = (TH1D*) h->ProjectionY("htmp",ii,ii);
    int iy = FindBinOfMax(htmp,iofl);
    double maxy = 0;
    double a = 0;
    double ymax = 0;
    double emaxy = fabs((htmp->GetBinCenter(1))-htmp->GetBinCenter(htmp->GetNbinsX()));
    //        cout << "*****************************************************" << endl;
    //      cout << ii << " " << integralslice << " " << iy << endl;

    if (iy>0) {maxy=htmp->GetBinCenter(iy);
      if (iy>i3 && iy<htmp->GetNbinsX()+1-i3) {
	double b,c,xx=0;
	lagrint(htmp->GetBinCenter(iy-i3),htmp->GetBinCenter(iy),htmp->GetBinCenter(iy+i3),htmp->GetBinContent(iy-i3),htmp->GetBinContent(iy),htmp->GetBinContent(iy+i3),xx,maxy,a,b,c);
	ymax = -b*b/(4.*a)+c;
      }

      emaxy = htmp->GetRMS();
      if (emaxy == 0 || integralslice <= 2) {maxy =0; emaxy=fabs((htmp->GetBinCenter(1))-htmp->GetBinCenter(htmp->GetNbinsX()));}
      else {
      if (irms==0) {emaxy = emaxy/sqrt(fabs(integralslice));}
      if (irms==1) {emaxy = emaxy;}
      if (irms==2 && a<0) {emaxy = sqrt(-ymax/(2.*a));}
      }}
    //       cout << ii << " " << maxy << " " << emaxy << " " << integralslice << " " << iy << " " << htmp->GetRMS() << endl;
    htmp->Delete();
    hx->SetBinContent(ii,maxy);
    hx->SetBinError(ii,escale*emaxy);
  }
  return 0;
}

int FindSep(TH1* sx, TH1* dx, double & ratio, double & thr, int iofl = 0, double approx = 0.005){
  // return the bin of separation threshold form sx and dx
  // return also the ratio of the confused tail and threshold value
  // IOFL = 1 IF YOU WANT OVERFLOW
  // approx is the approximation you have sx and dx with same area
  if (!sx || !dx) return -2;  // no histogram                                                                                   
  if(dx->GetNbinsX() != sx->GetNbinsX()) return -3; // not same bin number

  int nbin = sx->GetNbinsX();
  double  intsx=sx->Integral(1-iofl, nbin+iofl);
  double  intdx=dx->Integral(1-iofl, nbin+iofl);
  double tot = (intsx+intdx)/2.;
  ratio = fabs((intsx - intdx)/(2.*tot)); 
  if (ratio>approx) {
    cerr << " slices.C error; intsx = " << intsx << " differ from intdx " << intdx << endl;
    return -4; // not exactly same integral
  }
  double ms = sx->GetMean();
  double md = dx->GetMean();
  if (ms > md) return -5; // sx and dx switched
  int iret = 0;
  for (int ii=1-iofl; ii<=nbin+iofl; ii++) {
    intsx=sx->Integral(1-iofl, ii);
    intdx=tot-dx->Integral(1-iofl, ii);
    if (intdx<=intsx){
      iret = ii;
      break;}
  }

  ratio = (intsx+intdx)/(2.*tot);
  thr = sx->GetBinLowEdge(iret);
  return iret;
}


int FitRiseTimeX(TH2* h, TH1* hx, double ymin, double ymax){
  // For each X slice fill the hx with the position of risetime fitted with a gaussian from ymin to ymax

  if (!h) return 1;  // no histogram                                                           
  if(hx->GetNbinsX() != h->GetNbinsX()) return 3; // hx not correct                  
  if(ymax<=ymin) return 2; // limit not correct                  

  for (int ii=0; ii<=h->GetNbinsX()+1; ii++) {
    double risetime=0;
    double risetimerror=ymax-ymin;

    TH1D* htmp = (TH1D*) h->ProjectionY("htmp",ii,ii);
    int bmin = htmp->FindBin(ymin);
    int bmax = htmp->FindBin(ymax);
    if(bmax>bmin) {
      if(htmp->Integral(bmin,bmax)>0){
      htmp->GetXaxis()->SetRangeUser(ymin,ymax);
      htmp->Fit("gaus","Q");
      risetime = htmp->GetFunction("gaus")->GetParameter(1)-fabs(htmp->GetFunction("gaus")->GetParameter(2));
      risetimerror = sqrt(pow(htmp->GetFunction("gaus")->GetParError(1),2)+pow(htmp->GetFunction("gaus")->GetParError(2),2));
      }
    }
    htmp->Delete();
    hx->SetBinContent(ii,risetime);
    hx->SetBinError(ii,risetimerror);
  }
  return 0;
}
