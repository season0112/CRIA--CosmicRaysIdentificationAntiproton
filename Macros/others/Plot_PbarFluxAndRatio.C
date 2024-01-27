
#include "Plot_PbarFluxAndRatio.hh"

void Plot_PbarFluxAndRatio(){

    //// Files
    // PbarFlux
    TFile *f1 = new TFile("/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonFlux/AntiprotonFluxMulR_Compare_CRDB_above10.root", "READ"); // Flux * R^2.7, In R or E?
    TFile *f5 = new TFile("/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonFlux/AntiprotonFlux_Compare_CRDB.root","READ");      // Original Flux, In R
    TFile *f6 = new TFile("/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonFlux/AntiprotonFlux_Compare_CRDB_More.root","READ"); // Original Flux, In R
    TFile *f7 = new TFile("/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonFlux/AntiprotonFluxMulR_Compare_CRDB_More.root","READ"); // Flux * R^2.7 

    // PbarOverProton
    TFile *f2 = new TFile("/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonToProtonRatio/PbarOverProton_Compare_CRDB_above10.root", "READ"); 
    TFile *f3 = new TFile("/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/AntiprotonToProtonRatio/PbarOverProton_Compare_CRDB.root", "READ");  //In Rigidity
    // This analysis
    TFile *f4 = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultspass7.8.root", "READ");

    //// Load
    // Load PbarFlux Result (Flux * R^2.7, Above 10 GV)
    TGraphAsymmErrors *AMSPbarFluxMulR       = (TGraphAsymmErrors*)f1->Get("gr_exp1");
    TGraphAsymmErrors *CAPRICE98PbarFluxMulR = (TGraphAsymmErrors*)f1->Get("gr_exp2");
    TGraphAsymmErrors *MASS91PbarFluxMulR    = (TGraphAsymmErrors*)f1->Get("gr_exp3");
    TGraphAsymmErrors *PAMELAPbarFluxMulR    = (TGraphAsymmErrors*)f1->Get("gr_exp4");
    TGraphAsymmErrors *PAMELA2PbarFluxMulR   = (TGraphAsymmErrors*)f1->Get("gr_exp5");
    // Load PbarFlux Result (Flux, Full Range)
    TGraphAsymmErrors *AMS02Flux       = (TGraphAsymmErrors*)f5->Get("gr_exp1_errtot");
    TGraphAsymmErrors *BESSPolarIIFlux = (TGraphAsymmErrors*)f5->Get("gr_exp2_errtot");
    TGraphAsymmErrors *PAMELAFlux      = (TGraphAsymmErrors*)f5->Get("gr_exp3_errtot");

    RemoveXaxis(AMS02Flux);
    RemoveXaxis(BESSPolarIIFlux);
    RemoveXaxis(PAMELAFlux);

    TGraphAsymmErrors *BESSPolarIFlux   = (TGraphAsymmErrors*)f6->Get("gr_exp3");
    TGraphAsymmErrors *BESSTeVFlux      = (TGraphAsymmErrors*)f6->Get("gr_exp5");
    TGraphAsymmErrors *BESS00Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp6");
    TGraphAsymmErrors *BESS93Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp7");
    TGraphAsymmErrors *BESS95Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp8");
    TGraphAsymmErrors *BESS97Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp9");
    TGraphAsymmErrors *BESS98Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp10");
    TGraphAsymmErrors *BESS99Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp11");
    TGraphAsymmErrors *CAPRICE94Flux     = (TGraphAsymmErrors*)f6->Get("gr_exp12");
    TGraphAsymmErrors *CAPRICE98Flux     = (TGraphAsymmErrors*)f6->Get("gr_exp13");
    TGraphAsymmErrors *IMAX92Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp14");
    TGraphAsymmErrors *MASS91Flux        = (TGraphAsymmErrors*)f6->Get("gr_exp15");

    RemoveXaxis(BESSPolarIFlux);
    RemoveXaxis(BESSTeVFlux);
    RemoveXaxis(BESS00Flux);
    RemoveXaxis(BESS93Flux);
    RemoveXaxis(BESS95Flux);
    RemoveXaxis(BESS97Flux);
    RemoveXaxis(BESS98Flux);
    RemoveXaxis(BESS99Flux);
    RemoveXaxis(CAPRICE94Flux);
    RemoveXaxis(CAPRICE98Flux);
    RemoveXaxis(IMAX92Flux);
    RemoveXaxis(MASS91Flux);


    // Load PbarOverProton Result (Above 10 GV)
    TGraphAsymmErrors *PbarOverProton_CAPRICE98 = (TGraphAsymmErrors*)f2->Get("gr_exp1");
    TGraphAsymmErrors *PbarOverProton_HEATpbar  = (TGraphAsymmErrors*)f2->Get("gr_exp2");
    TGraphAsymmErrors *PbarOverProton_PAMELA1   = (TGraphAsymmErrors*)f2->Get("gr_exp3");
    TGraphAsymmErrors *PbarOverProton_PAMELA2   = (TGraphAsymmErrors*)f2->Get("gr_exp4");
    TGraphAsymmErrors *PbarOverProton_PAMELA3   = (TGraphAsymmErrors*)f2->Get("gr_exp5");
    // Load PbarOverProton Result (Full Range)
    TGraphAsymmErrors *ratioPAMELA    = (TGraphAsymmErrors*)f3->Get("graph3");

    // Load This analysis
    TGraphErrors      *ratio_pass78   = (TGraphErrors*)f4     ->Get("gRatio_unfolded");


    //// Plot
    //Plot_PbarOverProton_CompareInLowEnergy(ratioPAMELA             , ratio_pass78);
    //Plot_PbarOverProton_Compare_above10   (PbarOverProton_CAPRICE98, PbarOverProton_HEATpbar, PbarOverProton_PAMELA3);

    //Plot_PbarFluxMulR_Compare_above10     (AMSPbarFluxMulR         , CAPRICE98PbarFluxMulR  , MASS91PbarFluxMulR, PAMELA2PbarFluxMulR);
    Plot_PbarFlux_CRDB_Compare            (AMS02Flux, BESSPolarIIFlux, PAMELAFlux, BESSPolarIFlux, BESSTeVFlux, BESS00Flux, BESS93Flux, BESS95Flux, BESS97Flux, BESS98Flux, BESS99Flux, CAPRICE94Flux, CAPRICE98Flux, IMAX92Flux, MASS91Flux);

}


