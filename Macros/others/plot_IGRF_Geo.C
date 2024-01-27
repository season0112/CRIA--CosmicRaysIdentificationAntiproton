#include <stdlib.h>
#include <stdio.h> 
#include <time.h>
#include "Utilities.hh"
#include "plot_IGRF_Geo.hh"

int plot_IGRF_Geo() {

    //// Open File and get histograms
    TFile *GEOMETRIC     = new TFile("/hpcwork/jara0052/sichen/Measuringtime/MeasuringTime_pass7.8_05_2021_GEOMETRIC35_1.2.root");
    TH1D *GEOMETRIC_time                         = (TH1D*)GEOMETRIC->Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff");
    TH1D *GEOMETRIC_fMeasuringTimeVsLiveTime     = (TH1D*)GEOMETRIC->Get("MeasuringTime/fMeasuringTimeVsLiveTime");
    TH2D *GEOMETRIC_fParticlesVsTriggers         = (TH2D*)GEOMETRIC->Get("MeasuringTime/fParticlesVsTriggers");
    TH3F *GEOMETRIC_fTriggerRateVsPosition       = (TH3F*)GEOMETRIC->Get("MeasuringTime/fTriggerRateVsISSPosition");
    TH3F *GEOMETRIC_fCutOffRigidityVsISSPosition = (TH3F*)GEOMETRIC->Get("MeasuringTime/fCutOffRigidityVsISSPosition");
    TH3F *GEOMETRIC_fLiveTimeVsISSPosition       = (TH3F*)GEOMETRIC->Get("MeasuringTime/fLiveTimeVsISSPosition");
    //GEOMETRIC->Close();

    TFile *IGRF          = new TFile("/hpcwork/jara0052/sichen/Measuringtime/MeasuringTime_pass7.8_05_2021_IGRF35_1.2.root");
    TH1D *IGRF_time      = (TH1D*)IGRF->Get("MeasuringTime/fIntegratedMeasuringTimeOverCutOff");
    //IGRF->Close();

    //// Plot 
    //Plot_MeasuringTime_OnlyStormer    (GEOMETRIC_time);
    Plot_MeasuringTime              (GEOMETRIC_time, IGRF_time); // Stomer Vs IGRF 
    //Plot_MeasuringTimeVsLiveTime    (GEOMETRIC_fMeasuringTimeVsLiveTime);
    //Plot_ParticlesVsTriggers        (GEOMETRIC_fParticlesVsTriggers);
    //Plot_TriggerRateVsPosition      (GEOMETRIC_fTriggerRateVsPosition);
    //Plot_CutOffRigidityVsISSPosition(GEOMETRIC_fCutOffRigidityVsISSPosition);
    //Plot_LiveTimeVsISSPosition      (GEOMETRIC_fLiveTimeVsISSPosition); // projection to 2D ???

    std::cout << "\u0394V" << '\n';
    std::cout << "\u00F8" << '\n';
    std::cout << "\u0444" << '\n';

    return EXIT_SUCCESS;
}


