#include "Tofrec02_ihep.h"
#include "Tofcharge_ihep.h"

#include <EcalChi2CY.h>

#include "root_RVSP.h"
#include "amschain.h"
#include "root.h"
#include "TStopwatch.h"

#include "TrExtAlignDB.h"
#include "TrCharge.h"
#include "TofTrack.h"
#include "TofRecon.h"
#include "TkCoo.h"
#include "TrRecon.h"
#include "TrTrackSelection.h"
#include "TRandom.h"
#include "tkdcards.h"
#ifdef _PGTRACK_
#include <TrdSCalib.h>
#include <TrdKCluster.h>
#endif

#include <TH3D.h>
#include <signal.h>

#include "TApplication.h"




  bool TPreSelect(AMSSetupR::RTI &a){ 
    bool cut[10]={0};
    cut[0]=(a.ntrig/a.nev>0.98);
    cut[1]=(a.npart/a.ntrig>0.07/1600*a.ntrig&&a.npart/a.ntrig<0.25);
    cut[2]=(a.lf>0.5);
    cut[3]=(a.zenith<40);
    cut[4]=(a.nerr>=0&&a.nerr/a.nev<0.1);
    cut[5]=(a.npart>0&&a.nev<1800);
    cut[6]=(a.good==0);
    //  cut[7]=(!(a.IsInSAA()));
    //  bool tcut=(cut[0]&&cut[1]&&cut[2]&&cut[3]&&cut[4]&&cut[5]&&cut[6]&&cut[7]);
    bool tcut=(cut[0]&&cut[1]&&cut[2]&&cut[3]&&cut[4]&&cut[5]&&cut[6]);
    return tcut;
  }

static Int_t SigTERM = 0;


void handler(int sig)
{
  switch(sig){
  case SIGTERM: case SIGINT: 
    cerr << " SIGTERM intercepted" << endl;
    SigTERM = 1;
    break;
  }
}


int rck(const char *fname, char *fnout)
{



  double oldtime = 1.e9;
  double oldtimemusec = 1.e9;

  double pig = 3.141592654;
  Double_t emass = 510.998928e-6;
  Double_t pmass = 0.938272297;  // proton
  Double_t hemass = 3.727379;  // helium 
  AMSChain ch;
  ch.Add(fname);
  cout << " ******************* " << endl;

  Int_t ntr  = ch.GetNtrees();
  Int_t nent = ch.GetEntries();
  if (ntr <= 0 || nent <= 0) return 0;

  cout << "Ntr,Nent= " << ntr << " " << nent << endl;

  for (Int_t i = 0; i < ntr; i++)
    cout << ch.GetListOfFiles()->At(i)->GetTitle() << endl;

  static const int ICHAN = 50; //31;
  int ICHANT;
  Int_t nc[ICHAN];
  for (int ikk =0;ikk<ICHAN;ikk++)
    {nc[ikk]=1;}

  char reasons[ICHAN][300];
  TString oname = fnout;
  oname.ReplaceAll("root","log");
  ofstream logout;
  logout.open(oname.Data(),ios::out);
  

//  TFile of(fnout, "recreate");
  TFile *of = new TFile(fnout, "recreate");

  TTree *tree = new TTree("tree", "trcls");

  // recognize charge and momentum limit from file name of MC
  int namecharge = 0;
  float emin = 0;
  float emax = 8000;
  TString signtst = fname;
  if (signtst.Contains("pos") || signtst.Contains("prot")) namecharge = 1;
  if (signtst.Contains("el") || signtst.Contains("antiprot")) namecharge = -1;
  if (signtst.Contains("he")) namecharge = 2;
  if (signtst.Contains("0_255/")) {emin=0.25; emax=5.;}
  if (signtst.Contains("5100/")) {emin=5.; emax=100.;}
  if (signtst.Contains("22100/")) {emin=22.; emax=100.;}
  if (signtst.Contains("221000.")) {emin=22.; emax=1000.;}
  if (signtst.Contains("22800/")) {emin=22.; emax=800.;}
  if (signtst.Contains("1002000/")) {emin=100.; emax=2000.;}
  if (signtst.Contains("4472000/")) {emin=447.; emax=2000.;}
  if (signtst.Contains("0_510/")) {emin=0.5; emax=10.;}
  if (signtst.Contains("10200/")) {emin=10; emax=200.;}
  if (signtst.Contains("2004000/")) {emin=200; emax=4000.;}
  if (signtst.Contains("20400/")) {emin=20; emax=400.;}
  if (signtst.Contains("4008000/")) {emin=400; emax=8000.;}

  if (emin>0) {
    cout << "inferred from filename: MONTECARLO charge " << namecharge << " MONTECARLO RANGE " << emin << " - " << emax << endl;  }
  else { cout << "inferred from filename: DATA (NO MC) run " << endl;  }

  cout << "*************************************" << endl;
  static const int ILENGTH = 12;  
  Int_t   idata[ILENGTH];  
  tree->Branch("idata", idata, "nfil/I:run/I:event/I:ient/I:nrawside/I:ntofh/I:nanti/I:ntrdhits/I:"
	       "npairs/I:sector/I:adc0/I:adc1/I");

  static const int FLENGTH = 4+4+3+3+4+4+4+2+4+3+4+4;
  Float_t fdata[FLENGTH];
  tree->Branch("fdata", fdata, "timev/F:timusec/F:lat/F:lon/F:"
               "phimiss/F:AntiCrossX/F:AntiCrossY/F:AntiCrossZ/F:"
               "TrdTof0X/F:TrdTof0Y/F:TrdTof0Z/F:"
               "TrdTof1X/F:TrdTof1Y/F:TrdTof1Z/F:"
               "QTRD/F:Likep/F:LikeHe/F:LikHep/F:"
	       "Tusedtof0/F:Tusedtof1/F:Xusedtof0/F:Xusedtof1/F:"
	       "Yusedtof0/F:Yusedtof1/F:Zusedtof0/F:Zusedtof1/F:"
	       "Qusedtof0/F:Qusedtof1/F:"
	       "incl/F:dista/F:beta/F:phibar/F:"
	       "acchi2/F:zadc/F:unf_zadc/F:"
	       "uncal_time/F:uncal_zeta/F:uncal_qside0/F:uncal_qside1/F:"
	       "cal_time/F:cal_zeta/F:cal_qside0/F:cal_qside1/F"); 


  signal(SIGTERM, handler);
  signal(SIGINT,  handler);


  Int_t nevt = 0;
  Int_t nfil = 0;

  Int_t intv = 100000;

  TStopwatch timer;
  timer.Start();

  gRandom->SetSeed(0);



AMSSetupR::RTI::Version=AMSSetupR::RTI::UseLatest();
AMSSetupR::RTI rti;

 AntiRecoPG* Acci = AntiRecoPG::gethead();


Int_t iend=0;
for ( ULong64_t ient = 0; ient < nent && !SigTERM; ient++) {
Int_t indx=0;
// cout << ient << " inizio" << endl;
    AMSEventR *evt = ch.GetEvent(ient);

    int runtype = Acci->SelectRun();

    //  MC ONLY
   MCEventgR *pmc = evt->pMCEventg(0);
   nc[indx]++;  strcpy(reasons[indx],"MC without pMCE pointer"); indx++;   
   if (runtype == 99 && !pmc) continue;
   double MCM;
   double MCH;
   if(!pmc) {MCM = 1; MCH=0;} 
   else { MCM = pmc->Momentum; MCH = pmc->Charge;}
   nc[indx]++;  strcpy(reasons[indx],"MC momentum outside range"); indx++;   
   if (MCM<emin || MCM>emax) continue;

   nc[indx]++;  strcpy(reasons[indx],"MC charge different from filenamecharge"); indx++;   
   if (runtype == 99 && fabs(MCH-namecharge)>0.2) continue;
   double rgen = MCM;
   if (MCH!=0) rgen = rgen/MCH;

    if (nevt%intv == 0 || nevt == nent) {
      Double_t tm = timer.RealTime();
      timer.Continue();
      cout << Form("%6d %6d %7d (%5.1f%%) %4.0f sec (%4.1f kHz)",
		   ient, nfil, nevt,
		   100.*nevt/nent, tm, nevt/tm*1e-3)
	   << endl;
      logout <<  Form("%6d %6d %6d %7d (%5.1f%%) %4.0f sec (%4.1f kHz)",
	nc[0], nc[15], nfil, nevt,
		   100.*nevt/nent, tm, nevt/tm*1e-3)
	   << endl;
      logout << flush;
    }
   
    // done first raw compatibility selection for MC
  
   nevt++;

   bool oldRTIpass = true;
    double ttime = evt->fHeader.Time[0];
    if(runtype == 1) {  // RTI only for DATA
     if(ttime!=oldtime){
       if (AMSEventR::GetRTI(rti,ttime)!=0) continue; //Time no information
       oldRTIpass=TPreSelect(rti);
     }}

    double ttimemusec = evt->fHeader.Time[1];
    //    double dtime = ttime-oldtime;
    //    dtime = dtime + 1.e-6*(ttimemusec-oldtimemusec);
    oldtime = ttime;
    oldtimemusec = ttimemusec;
   nc[indx]++;  strcpy(reasons[indx],"RTI BAD SECOND"); indx++;
   if(!oldRTIpass) continue;


   bool badruntag = false;
    if(runtype == 1) {  // RTI only for DATA
      //      if ((evt->fHeader.RunType)<0xF016) badruntag=true;
      //      cout << std::hex << evt->fHeader.RunType << std::dec << endl;
    }
   nc[indx]++;  strcpy(reasons[indx],"BAD RUNTAG"); indx++;
   if(badruntag) continue;


   nc[indx]++;  strcpy(reasons[indx],"no trdhtrk"); indx++;
   if (evt->nTrdHTrack() < 1) continue;

   nc[indx]++;  strcpy(reasons[indx],"! anti found"); indx++;
   if (evt->nAntiRawSide()==0) continue;

   TofRecH::ReBuild();
   nc[indx]++;  strcpy(reasons[indx],"! tof found"); indx++;
   if (evt->nTofClusterH()==0) continue;

   AntiClusterR* pcl;
   TrdHTrackR* htrd=evt->pTrdHTrack(0);
   AMSPoint Point1;
   AMSPoint Point2;
   AMSPoint DTRD;  // TRD h trk direction 
   AMSPoint AntiCross;
   int iposec;
   float phimiss;

   double dusedtof[2] ={1000,1000};
   double tusedtof[2] ={0,0};
   double xusedtof[2] ={0,0};
   double yusedtof[2] ={0,0};
   double zusedtof[2] ={0,0};
   double qusedtof[2] ={0,0};
   AMSPoint trdtof[2];

   Double_t Likep = -1000;
   Double_t LikHep = -1000;
   Double_t LikeHe = -1000;
   double QTRD=0.1;
   if (runtype==99) QTRD=fabs(MCH);
   int trdkhit = 0;
   bool goodTrdKCalib = true;
   double kLike[3]={-1,-1,-1}; //To be filled with 3 LikelihoodRatio : e/P, e/H, P/H        
   int NHits=0;                //To be filled with number of hits taken into account in Likelihood Calculation
   float threshold=15;         //min ADC taken into account in Likelihood Calculation, 15 ADC is the recommended value for the moment.     
   double phibar;

   int hitsec=0;


   for(int i=0; i<evt->nTrdHTrack(); i++){ // scan trdtrk if hitting acc       
     // FOUND THE TRD TRK
     htrd=evt->pTrdHTrack(i);
     if( !htrd ) continue;
     Point1.setp(htrd->Coo[0],htrd->Coo[1],htrd->Coo[2]);
     Point2.setp(Point1.x()+htrd->Dir[0],Point1.y()+htrd->Dir[1],Point1.z()+htrd->Dir[2]);
     DTRD.setp(htrd->Dir[0],htrd->Dir[1],htrd->Dir[2]);  // TRD trk direction

     //TEST IF TRD TRK IS CROSSING AN ANTICOUNTER                                                
     iposec=Acci->GetCrossing(& Point1, & Point2, & AntiCross, phimiss);
     if (iposec<=0) continue; // not acc crossing    

     // SCAN TOF CLUSTERS searching for matching with the TRD trk 
     dusedtof[0] = 1000;
     dusedtof[1] = 1000;
     for(int i=0; i<evt->nTofClusterH(); i++){
       TofClusterHR* tclu = evt->pTofClusterH(i);
       if(!tclu) continue;
       if(tclu->IsGoodTime()==false) continue; //if true have a signal in both sides  
       int tlayer = tclu->Layer;
       if (tlayer > 1) continue;  // only upper TOF layers 
       double cttime= tclu->Time;
       double cx = tclu->Coo[0];
       double cy = tclu->Coo[1];
       double cz = tclu->Coo[2];
       double tstar = (cz-Point1.z())/DTRD.z();
       double tx = Point1.x()+DTRD.x()*tstar;  // estrapolation X TRDTRK @ZTOF              
       double ty = Point1.y()+DTRD.y()*tstar;  // estrapolation Y TRDTRK @ZTOF              
       double dd = sqrt((tx-cx)*(tx-cx)+(ty-cy)*(ty-cy));  // distance of extrapolation
       if (dd<dusedtof[tlayer]){  //find the nearest cluster to trk                       
	 dusedtof[tlayer]=dd;
	 tusedtof[tlayer]=cttime;
	 xusedtof[tlayer]=cx;
	 yusedtof[tlayer]=cy;
	 zusedtof[tlayer]=cz;
	 qusedtof[tlayer]=tclu->GetQSignal(2,TofClusterHR::DefaultQOpt,fabs(DTRD.z()));  //TOF charge corrected only for inclination
	 trdtof[tlayer].setp(tx,ty,cz);
       }
     }

     int imatch=1;
     for (int it=0;it<2;it++) {if (dusedtof[it]>10.) imatch=0;}  //10 cm window
     if (imatch==0) continue;  // not tof crossed pass to next trk     

     // find the charge for TRD in case of DATA/BT
     if(evt->NTrdRawHit()>0 &&  (runtype<3) )
       {
         TrdKCluster* trdkcluster = new TrdKCluster(evt, htrd );
         // Get the status of the Calculation                                        
         int trdKIsReadAlignmentOK=trdkcluster->IsReadAlignmentOK;
         int trdKIsReadCalibOK=trdkcluster->IsReadCalibOK;
         int isvalid = trdkcluster->GetLikelihoodRatio_TRDRefit(threshold,kLike,NHits);
         if( NHits>0 && trdKIsReadAlignmentOK>0 && trdKIsReadCalibOK>0 && isvalid )
           {
             goodTrdKCalib=true;
             Likep = kLike[0];
             LikeHe = kLike[1];
             LikHep = kLike[2];
             int QTRDStatus = trdkcluster->CalculateTRDCharge(0, 1);
             if (QTRDStatus>=0) QTRD = trdkcluster->GetTRDCharge();
	     if (! (fabs(QTRD)<1.e6)) QTRD=0.1; // to manage null/inf
           }
         else
           {
             goodTrdKCalib=false;
           }
         delete trdkcluster;
       }
     else
       {
         goodTrdKCalib=false;
       }
     trdkhit = NHits;


     //PRESELECTION MATCHING OF TRDTRK and UPPER TOF
     // MAYBE FOR STD STREAM YOU SHOULD STUDY THIS NUMBERS
     if (!(fabs(QTRD)<30)) continue;  // charge nonsense
     if (fabs(qusedtof[0]-qusedtof[1])>0.5) continue;  // QTOF self-consistent   
     if (fabs(tusedtof[0]-tusedtof[1])>1.5) continue;  // timing TOF self-consistent
     if (fabs(1.15*QTRD-0.5*(qusedtof[0]+qusedtof[1]))>0.3) continue; //QTRD and QTOF matching (check)
     if(fabs(trdtof[0].x()-xusedtof[0])>5)continue;    // toftrk and tofposition matching
     if(fabs(trdtof[0].y()-yusedtof[0])>7)continue;
     if(fabs(trdtof[1].x()-xusedtof[1])>7)continue;
     if(fabs(trdtof[1].y()-yusedtof[1])>5)continue;


     phibar = (iposec*45.-22.5); // center of the hitted acc bar

     // now see if the BAR have signal
     evt->RebuildAntiClusters(iposec,AntiCross.z());
     for (int icl = 0; icl<evt->nAntiCluster();icl++){ //scan acc for hitted from trd                                         
       pcl = evt->pAntiCluster(icl);
       if (!pcl) continue;
       int bar = pcl->Sector;
       if (bar != iposec) continue;
       if(pcl->adc[0]==0 && pcl->adc[1]==0) continue;  // if no signal on both side not interesting
       if(pcl->chi>100) continue;  // if chisquare is too large probably is casual, then not interesting
       hitsec=iposec;
       break;}
     if (hitsec != 0) break; // found an hitted acc                                                         
   }

   nc[indx]++;  strcpy(reasons[indx],"! anti hitted and fired"); indx++;
   if (hitsec==0) continue;

   // GET QUANTITY WITH and WITHOUT CORRECTIONS
   int sector = pcl->Sector;
   double accchi2=pcl->chi;
   int np=pcl->Npairs;
   int adc0=pcl->adc[0];
   int adc1=pcl->adc[1];

   double zadc=pcl->zeta_adc;  // z from ADC
   double unf_zadc=pcl->unfzeta_adc; // unfolded z from ADC

   double uncal_time = Acci->GetTime(sector,0,1,0,phimiss);
   double cal_time = Acci->GetTime(sector,0,1,1,phimiss);
   double uncal_zeta = Acci->GetZeta(sector,0,0); 
   double cal_zeta = Acci->GetZeta(sector,0,1); 

   // UPPER TOF avg time
   float exttime = 0.5*(tusedtof[0]+tusedtof[1]);
   // UPPER TOF extrapolated position
   AMSPoint trdtofpoint;
   trdtofpoint.setp(0.5*(trdtof[0].x()+trdtof[1].x()),0.5*(trdtof[0].y()+trdtof[1].y()),0.5*(trdtof[0].z()+trdtof[1].z()));
   //find BETA and inclination from tracking and from timing                   
   float beta;
   float incl;
   float dista = pcl->DoBeta(&trdtofpoint,exttime,beta,incl,&AntiCross);

   float beta_corr = beta;

   // of course if time is badly calibrated or not calibrated you could force to +-1
   //   if (beta_corr>0) beta_corr = 1;
   //   if (beta_corr<0) beta_corr = -1;


   // without corrections
   double uncal_qside0=pcl->Charge(incl,beta_corr,0,AntiCross[2],0);
   double uncal_qside1=pcl->Charge(incl,beta_corr,1,AntiCross[2],0);

   // with corrections
   double cal_qside0=pcl->Charge(incl,beta_corr,0,AntiCross[2]);
   double cal_qside1=pcl->Charge(incl,beta_corr,1,AntiCross[2]);


   float thetas = evt->fHeader.ThetaS;  // [rad]       
   float phis = evt->fHeader.PhiS;      // [rad]         
   double lat = thetas*180/TMath::Pi();         // [ -90,  90]             
   double lon = phis*180/TMath::Pi();           // [   0, 360]              
   if (lon>=180) lon = (lon-360);      // [-180, 180]  

    int iind = 0;
    idata[iind] = nfil; iind++;
    idata[iind] = evt->Run(); iind++;
    idata[iind] = evt->Event(); iind++;
    idata[iind] = ient; iind++;
    idata[iind] = evt->nAntiRawSide(); iind++;
    idata[iind] = evt->nTofClusterH(); iind++;
    idata[iind] = evt->nAntiCluster(); iind++;
    idata[iind] = NHits; iind++;

    idata[iind] = np; iind++;
    idata[iind] = sector; iind++;
    idata[iind] = adc0; iind++;
    idata[iind] = adc1; iind++;



    if (iind != ILENGTH) {
      cout << ILENGTH <<" != " << iind << endl;
      logout << ILENGTH <<" != " << iind << endl;
      return 0;
    }
    

    iind = 0;
    fdata[iind] = ttime;    iind++;
    fdata[iind] = ttimemusec;    iind++;
    fdata[iind] = lat;    iind++;
    fdata[iind] = lon;    iind++;

    fdata[iind] = phimiss; iind++;
    fdata[iind] = AntiCross.x(); iind++;
    fdata[iind] = AntiCross.y(); iind++;
    fdata[iind] = AntiCross.z(); iind++;

    fdata[iind] = trdtof[0].x(); iind++;
    fdata[iind] = trdtof[0].y(); iind++;
    fdata[iind] = trdtof[0].z(); iind++;
    fdata[iind] = trdtof[1].x(); iind++;
    fdata[iind] = trdtof[1].y(); iind++;
    fdata[iind] = trdtof[1].z(); iind++;
    fdata[iind] = QTRD; iind++;
    fdata[iind] = Likep; iind++;
    fdata[iind] = LikeHe; iind++;
    fdata[iind] = LikHep; iind++;

    for (int iii=0;iii<2;iii++){
      fdata[iind] = tusedtof[iii]; iind++;}
    for (int iii=0;iii<2;iii++){
      fdata[iind] = xusedtof[iii]; iind++;}
    for (int iii=0;iii<2;iii++){
      fdata[iind] = yusedtof[iii]; iind++;}
    for (int iii=0;iii<2;iii++){
      fdata[iind] = zusedtof[iii]; iind++;}
    for (int iii=0;iii<2;iii++){
      fdata[iind] = qusedtof[iii]; iind++;}

    fdata[iind] = incl;    iind++;
    fdata[iind] = dista;    iind++;
    fdata[iind] = beta;    iind++;
    fdata[iind] = phibar;    iind++;

    fdata[iind] = accchi2; iind++;
    fdata[iind] = zadc; iind++;
    fdata[iind] = unf_zadc; iind++;
    
    fdata[iind] = uncal_time; iind++;
    fdata[iind] = uncal_zeta; iind++;
    fdata[iind] = uncal_qside0; iind++;
    fdata[iind] = uncal_qside1; iind++;

    fdata[iind] = cal_time; iind++;
    fdata[iind] = cal_zeta; iind++;
    fdata[iind] = cal_qside0; iind++;
    fdata[iind] = cal_qside1; iind++;




    if (iind != FLENGTH) {
      cout << FLENGTH <<" != " << iind << endl;
      logout << FLENGTH <<" != " << iind << endl;
      return 0;
    }

    
      tree->Fill();
    nfil++;
    ICHANT = indx;
 }
  of->cd();
      
  tree->Write();
      
  //      logout << 0 << " " <<  nc[0] << " " << (float) nc[0]/nevt << " " << nevt << " " << nent << endl;
  //logout << 0 << " " << reasons[0] << endl;
  for (int ikk =1;ikk<=ICHANT;ikk++)
    {
logout << ikk << " " << reasons[ikk-1] << endl;
      logout << ikk << " " <<  nc[ikk] << " " << (float) nc[ikk]/nc[ikk-1] << " " << nevt << " " << nent << endl;

    }
  logout << nfil << endl;
      logout << flush;
      //  cout << " here i am " << iend << endl;
      //      TCC->~TrkCC();
if (iend !=0) return 2;
  return 1;
}


int main(int argc,char ** argv){

  char fname[255]; 
  char fnout[255]; 

  // int iout=atoi(argv[2]);


 if((argc==3))
{sprintf(fname,"%s",argv[1]);
sprintf(fnout,"%s",argv[2]);
  }
  else{
     printf("Usage: %s  <file_to_analyze> <out filename> \n",argv[0]);
     exit(-1);
   } 
 signal(SIGTERM, handler);
 signal(SIGINT, handler);
 printf("doing: %s %s %s \n",argv[0],argv[1],argv[2]);

 TApplication theApp("App", &argc, argv);

 ofstream fout;
  fout.open("bad.dat",ios::app);
  fout << fnout << endl;
  fout << flush;
  fout.close();


 int iret = 0;
iret = rck(fname,fnout);
// cout << " ecco iret " << iret << endl;
 if (iret == 1){
 ofstream fout;
 fout.open("finish.dat",ios::app);
 //  cout << " opened " << iret << endl;
  fout << fnout << endl;
  //  cout << " writtened " << iret << endl;
  fout << flush;
  //  cout << " flushed " << iret << endl;
  fout.close();
  //  cout << " closed " << iret << endl;
  system("rm -rf bad.dat");
  system("touch bad.dat");
 }

  return !iret; }
