//Usage: root -L /home/bo791269/Software/AntiprotonAnalysis/Libraries/AntiprotonBinning/AntiprotonBinning.C -x Plot_FullRangeRatio.C
//Note: Sys error(ACC) is relative error from acceptance calculation. then times ratio; Sys error(CC) is load with unfolding root result with absolute value, then divide by ratio into relative one.

#include "Utilities.hh"
#include "Plot_FullRangeRatio.hh"

void Plot_FullRangeRatio(){

    //std::string issversion = "2016paper";
    //std::string issversion = "PhyRep2021";
    std::string issversion = "pass7.8";

    std::string issversionname;
    if (issversion == "pass7.8"){
        issversionname = "pass78";        
    }
    else if (issversion == "2016paper"){
        issversionname = "2016paper";
    }
    else if (issversion == "PhyRep2021"){
        issversionname = "PhyRep2021";
    }

    double Scaler = 1.0;


    //// Load PhysicsReport result
    TFile *f_phyreport = new TFile( (string(getenv("MY_ANALYSIS")) + string("/ReferenceFiles/AntiprotonToProtonRatio/AntiprotonToProtonRatio_PhysicsReport/ssdc_canvas.root")).c_str() );
    TGraphErrors *gPhyrReortRatio = (TGraphErrors*)f_phyreport->Get("graph1");
    

    //// Load Low range and PRL2016 result
    string LowFilename;
    if (issversion == "pass7.8"){
        LowFilename     = string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultspass7.8.root");
        //LowFilename     = string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots_old/unfolded_resultspass7.8.root");
        //LowFilename     = string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v10.0") + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultspass7.8.root");
    }
    else if (issversion == "2016paper"){
        LowFilename     = string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_results2016paper.root");
        //LowFilename     = string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots_old/unfolded_results2016paper.root");
        //LowFilename     = string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v10.0") + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_results2016paper.root");
    }
    else if (issversion == "PhyRep2021"){
        LowFilename = string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultsPhysicsReport.root");
        //LowFilename = string(getenv("HPCLOWENERGYDIR")) + string("/totalall/Time_Averaged_ratio_Low/plots_old/unfolded_resultsPhysicsReport.root");
        //LowFilename     = string("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v10.0") + string("/totalall/Time_Averaged_ratio_Low/plots/unfolded_resultsPhysicsReport.root");
    }
    TFile *LowFilenameFile     = new TFile(LowFilename.c_str());
    // Ratio
    TGraphErrors *g_RawLowResult  = (TGraphErrors*)LowFilenameFile->Get("gRatio_Raw");
    TGraphErrors *g_LowResult     = (TGraphErrors*)LowFilenameFile->Get("gRatio_unfolded");
    TGraphErrors *gPublishedRatio = (TGraphErrors*)LowFilenameFile->Get("gPublishedRatio");
    // Number
    TGraph *g_antiproton_number_raw_Low      = (TGraph *)LowFilenameFile->Get("g_antiproton_number_raw");
    TGraph *g_Proton_number_raw_Low          = (TGraph *)LowFilenameFile->Get("g_proton_number_raw");
    TH1D *h_Antiproton_number_unfolded_Low   = (TH1D*)LowFilenameFile   ->Get("hAntiproton_number_unfolded");
    TH1D *h_Proton_number_unfolded_Low       = (TH1D*)LowFilenameFile   ->Get("hProton_number_unfolded");
    TGraph *g_Antiproton_number_unfolded_Low = new TGraph(h_Antiproton_number_unfolded_Low);
    TGraph *g_Proton_number_unfolded_Low     = new TGraph(h_Proton_number_unfolded_Low);
    // EffectiveAcceptance
    TGraphErrors *Effective_Acceptance_Low_B1042              = (TGraphErrors*)LowFilenameFile->Get("g_Effective_Acceptance");
    // Error
    TH1D *h_SystematicError_Low                         = (TH1D*)LowFilenameFile->Get("h_SystematicError");
    TGraphErrors *g_SystematicError_Shape_Low           = (TGraphErrors*)LowFilenameFile->Get("g_SystematicError");
    TGraphErrors *g_StatisticalError_Low                = (TGraphErrors*)LowFilenameFile->Get("g_StatisticalError");
    TGraphErrors *g_StatisticalRelError_Low             = (TGraphErrors*)LowFilenameFile->Get("g_StatisticalRelError");
    TGraph *g_PublishedRatioToralError                  = (TGraph*)LowFilenameFile->Get("g_PublishedError");
    TGraph *g_PublishedRatioStatisticErrorPRL           = (TGraph*)LowFilenameFile->Get("g_PublishedRatioStatisticErrorPRL");
    TGraph *g_PublishedRatioSystematicErrorPRL          = (TGraph*)LowFilenameFile->Get("g_PublishedRatioSystematicErrorPRL");
    TGraph *g_PublishedRatioRelativeErrorPRL            = (TGraph*)LowFilenameFile->Get("g_PublishedRatioRelativeErrorPRL");
    TGraph *g_PublishedRatioStatisticRelativeErrorPRL   = (TGraph*)LowFilenameFile->Get("g_PublishedRatioStatisticRelativeErrorPRL");
    TGraph *g_PublishedRatioSystematicRelativeErrorPRL  = (TGraph*)LowFilenameFile->Get("g_PublishedRatioSystematicRelativeErrorPRL");
    TGraph *g_PhysicsReportRatioError                   = (TGraph*)LowFilenameFile->Get("g_PhysicsReportRatioError");
    TGraph *g_PhysicsReportRatioStatisticError          = (TGraph*)LowFilenameFile->Get("g_PhysicsReportRatioStatisticError");
    TGraph *g_PhysicsReportRatioSystematicError         = (TGraph*)LowFilenameFile->Get("g_PhysicsReportRatioSystematicError");
    TGraph *g_PhysicsReportRatioRelativeError           = (TGraph*)LowFilenameFile->Get("g_PhysicsReportRatioRelativeError");
    TGraph *g_PhysicsReportStatisticRatioRelativeError  = (TGraph*)LowFilenameFile->Get("g_PhysicsReportStatisticRatioRelativeError");
    TGraph *g_PhysicsReportSystematicRatioRelativeError = (TGraph*)LowFilenameFile->Get("g_PhysicsReportSystematicRatioRelativeError");
    TH1D *h_StatisticalRelError_Low                     = (TH1D*)LowFilenameFile->Get("h_StatisticalRelError");

    // Binning
    std::vector<double> PublishedPRLBinEdge ( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.end());            // since first point (0.8-1.0) don't have result, therefore it should be removed.
    std::vector<double> PhysicsReportBinEdge ( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_zhili525.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_zhili525.end()); // since first point (0.8-1.0) don't have result, therefore it should be removed.
    // Convert Graph to Histogram
    TH1D h_PublishedRatioStatisticRelativeErrorPRL = TH1D("", "", g_PublishedRatioStatisticRelativeErrorPRL->GetN(), PublishedPRLBinEdge.data());
    Utilities::ConvertToHistogram ( g_PublishedRatioStatisticRelativeErrorPRL, h_PublishedRatioStatisticRelativeErrorPRL);
    TH1D h_PhysicsReportStatisticRatioRelativeError = TH1D("", "", g_PhysicsReportStatisticRatioRelativeError->GetN(), PhysicsReportBinEdge.data());
    Utilities::ConvertToHistogram ( g_PhysicsReportStatisticRatioRelativeError, h_PhysicsReportStatisticRatioRelativeError);

    TH1D h_PublishedRatioSystematicRelativeErrorPRL = TH1D("", "", g_PublishedRatioSystematicRelativeErrorPRL->GetN(), PublishedPRLBinEdge.data());
    Utilities::ConvertToHistogram ( g_PublishedRatioSystematicRelativeErrorPRL, h_PublishedRatioSystematicRelativeErrorPRL);
    TH1D h_PhysicsReportSystematicRatioRelativeError = TH1D("", "", g_PhysicsReportSystematicRatioRelativeError->GetN(), PhysicsReportBinEdge.data());
    Utilities::ConvertToHistogram ( g_PhysicsReportSystematicRatioRelativeError, h_PhysicsReportSystematicRatioRelativeError);

    TH1D h_PublishedRatioRelativeErrorPRL = TH1D("", "", g_PublishedRatioRelativeErrorPRL->GetN(), PublishedPRLBinEdge.data());
    Utilities::ConvertToHistogram ( g_PublishedRatioRelativeErrorPRL, h_PublishedRatioRelativeErrorPRL);
    TH1D h_PhysicsReportRatioRelativeError = TH1D("", "", g_PhysicsReportRatioRelativeError->GetN(), PhysicsReportBinEdge.data());
    Utilities::ConvertToHistogram ( g_PhysicsReportRatioRelativeError, h_PhysicsReportRatioRelativeError);

    // StatisticalError in Low range change back to absolute value.
    for (int q = 0; q <g_StatisticalError_Low->GetN(); q++){
        g_StatisticalError_Low->SetPoint(q, g_StatisticalError_Low->GetX()[q], g_StatisticalError_Low->GetY()[q]/100000);
    }

    //// Load Intermediate range
    string IntermediateFilename;
    string IntermediateFilename_UnfoldingIterationMore = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_UnfoldingIteration_6pass7.8.root"); // Here for testing:  give different Unfolding Iteration Numbers.
    if (issversion == "pass7.8"){
        IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_pass7.8.root");
        //IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding_old/unfolded_results_pass7.8.root");
        //IntermediateFilename = string("/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v5.0") + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_pass7.8.root");
    }
    else if (issversion == "2016paper"){
        IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_2016paper.root");
        //IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding_old/unfolded_results_2016paper.root");
        //IntermediateFilename = string("/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v5.0") + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_2016paper.root");
    }
    else if (issversion == "PhyRep2021"){
        IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_PhysicsReport.root");
        //IntermediateFilename = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding_old/unfolded_results_PhysicsReport.root");
        //IntermediateFilename = string("/hpcwork/jara0052/sichen/AntiprotonIntermediateEnergy_v5.0") + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_PhysicsReport.root");
    }
    TFile *IntermediateFile                        = new TFile(IntermediateFilename.c_str());
    TFile *IntermediateFile_UnfoldingIterationMore = new TFile(IntermediateFilename_UnfoldingIterationMore.c_str());
    // Ratio
    TGraphErrors *g_RawIntermediateResult  = (TGraphErrors*)IntermediateFile->Get("gRatio_Raw");
    TGraphErrors *g_IntermediateResult     = (TGraphErrors*)IntermediateFile->Get("gRatio_unfolded");
    TGraphErrors *g_RawIntermediateResult_UnfoldingIterationMore = (TGraphErrors*)IntermediateFile_UnfoldingIterationMore->Get("gRatio_Raw");
    TGraphErrors *g_IntermediateResult_UnfoldingIterationMore    = (TGraphErrors*)IntermediateFile_UnfoldingIterationMore->Get("gRatio_unfolded");
    // Numbers
    TH1D *h_Antiproton_number_unfolded_Intermediate   = (TH1D*)IntermediateFile->Get("hAntiproton_number_unfolded");
    TH1D *h_Proton_number_unfolded_Intermediate       = (TH1D*)IntermediateFile->Get("hProton_number_unfolded");
    TH1D *h_Antiproton_number_raw_Intermediate        = (TH1D*)IntermediateFile->Get("hantiproton_number_raw");
    TH1D *h_Proton_number_raw_Intermediate            = (TH1D*)IntermediateFile->Get("hproton_number_raw");
    TGraph *g_Antiproton_number_unfolded_Intermediate = new TGraph(h_Antiproton_number_unfolded_Intermediate);
    TGraph *g_Proton_number_unfolded_Intermediate     = new TGraph(h_Proton_number_unfolded_Intermediate);
    TGraph *g_Antiproton_number_raw_Intermediate      = new TGraph(h_Antiproton_number_raw_Intermediate);
    TGraph *g_Proton_number_raw_Intermediate          = new TGraph(h_Proton_number_raw_Intermediate);
    // Acceptance
    TGraphErrors *Effective_Acceptance_Intermediate_B1042        = (TGraphErrors*)IntermediateFile->Get("Effective_Acceptance_IntermediateOnly");
    TGraphErrors *Effective_Acceptance_Intermediate_all_B1042    = (TGraphErrors*)IntermediateFile->Get("Effective_Acceptance");
    TGraphErrors *Effective_Acceptance_Intermediate_inHigh_B1042 = (TGraphErrors*)IntermediateFile->Get("Effective_Acceptance_Intermediate_inHigh");
    // Statistical Error
    TGraph *g_StatisticError_Intermediate            = (TGraphErrors*)IntermediateFile->Get("g_StatisticError");
    TGraph *g_StatisticalRelError_Intermediate       = (TGraphErrors*)IntermediateFile->Get("g_StatisticalRelError");
    TH1D *h_StatisticalRelError_Intermediate         = (TH1D*)IntermediateFile->Get("h_StatisticalRelError");
    // Systematic Error (Taken from pass7.8, FIXME)
    string IntermediateFilename_Pass78 = string(getenv("HPCINTERMEDIATEDIR")) + string("/total/Time_Averaged_ratio/unfolding/unfolded_results_pass7.8.root");     
    TFile *IntermediateFile_Pass78 = new TFile(IntermediateFilename_Pass78.c_str());
    TGraphErrors *g_SystematicError_TRD_Intermediate = (TGraphErrors*)IntermediateFile_Pass78->Get("g_SystematicError_TRD");
    TH1D         *h_SystematicRelativeError_TRD      = (TH1D*)IntermediateFile_Pass78->Get("h_SystematicError_TRD");  // To be refilled with relative error.

    //// Load High range
    string HighFilename_MixedPattern;
    if (issversion == "pass7.8"){
        HighFilename_MixedPattern = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_MixedPattern.root");
    }
    else if(issversion == "2016paper"){
        HighFilename_MixedPattern = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_MixedPattern_May2015.root");
    }
    else if (issversion == "PhyRep2021"){
        HighFilename_MixedPattern = string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/unfolding/unfolded_results_MixedPattern_Nov2017.root");
    }
    TFile *highfile_MixedPattern = new TFile(HighFilename_MixedPattern.c_str());
    // Ratio
    TGraphErrors *g_HighResult         = (TGraphErrors*)highfile_MixedPattern->Get("g_Ratio_unfolded_with_AcceptanceCorrection_mixedPattern");
    TGraphErrors *g_HighResult_P0      = (TGraphErrors*)highfile_MixedPattern->Get("g_HighResult_P0");
    TGraphErrors *g_HighResult_P0VGGNN = (TGraphErrors*)highfile_MixedPattern->Get("g_HighResult_P0VGGNN");
    TGraphErrors *g_HighResult_P1      = (TGraphErrors*)highfile_MixedPattern->Get("g_HighResult_P1");
    TGraphErrors *g_HighResult_P2      = (TGraphErrors*)highfile_MixedPattern->Get("g_HighResult_P2");
    TGraphErrors *g_HighResult_P3      = (TGraphErrors*)highfile_MixedPattern->Get("g_HighResult_P3");
    TGraphErrors *g_HighResult_P4      = (TGraphErrors*)highfile_MixedPattern->Get("g_HighResult_P4");
    TGraphErrors *g_HighResult_P5      = (TGraphErrors*)highfile_MixedPattern->Get("g_HighResult_P5");
    TGraphErrors *g_Raw_Unfold_Compare_P0VGGNN = (TGraphErrors*)highfile_MixedPattern->Get("g_Raw_Unfold_Compare_P0VGGNN");  
    TGraphErrors *g_Raw_Unfold_Compare_P1      = (TGraphErrors*)highfile_MixedPattern->Get("g_Raw_Unfold_Compare_P1");
    TGraphErrors *g_Raw_Unfold_Compare_P2      = (TGraphErrors*)highfile_MixedPattern->Get("g_Raw_Unfold_Compare_P2");
    TGraphErrors *g_Raw_Unfold_Compare_P4      = (TGraphErrors*)highfile_MixedPattern->Get("g_Raw_Unfold_Compare_P4");

    // Number
    TGraph *g_Antiproton_number_unfolded_High         = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_mixedPattern") );
    TGraph *g_Proton_number_unfolded_High             = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_mixedPattern") );
    TGraph *g_Antiproton_number_raw_High              = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_mixedPattern") );
    TGraph *g_Proton_number_raw_High                  = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_mixedPattern") );

    TGraph *g_Antiproton_number_unfolded_High_P0      = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P0") );
    TGraph *g_Antiproton_number_unfolded_High_P0VGGNN = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P0VGGNN") );
    TGraph *g_Antiproton_number_unfolded_High_P1      = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P1") );
    TGraph *g_Antiproton_number_unfolded_High_P2      = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P2") );
    TGraph *g_Antiproton_number_unfolded_High_P3      = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P3") );
    TGraph *g_Antiproton_number_unfolded_High_P4      = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P4") );
    TGraph *g_Antiproton_number_unfolded_High_P5      = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P5") );

    TGraph *g_Antiproton_number_raw_High_P0           = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_P0") );
    TGraph *g_Antiproton_number_raw_High_P0VGGNN      = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_P0VGGNN") );
    TGraph *g_Antiproton_number_raw_High_P1           = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_P1") );
    TGraph *g_Antiproton_number_raw_High_P2           = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_P2") );
    TGraph *g_Antiproton_number_raw_High_P3           = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_P3") );
    TGraph *g_Antiproton_number_raw_High_P4           = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_P4") );
    TGraph *g_Antiproton_number_raw_High_P5           = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_Raw_High_P5") );

    TGraph *g_Proton_number_unfolded_High_P0          = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_P0") );
    TGraph *g_Proton_number_unfolded_High_P0VGGNN     = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_P0VGGNN") );
    TGraph *g_Proton_number_unfolded_High_P1          = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_P1") );
    TGraph *g_Proton_number_unfolded_High_P2          = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_P2") );
    TGraph *g_Proton_number_unfolded_High_P3          = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_P3") );
    TGraph *g_Proton_number_unfolded_High_P4          = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_P4") );
    TGraph *g_Proton_number_unfolded_High_P5          = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_unfolded_High_P5") );

    TGraph *g_Proton_number_raw_High_P0               = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_P0") );
    TGraph *g_Proton_number_raw_High_P0VGGNN          = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_P0VGGNN") );
    TGraph *g_Proton_number_raw_High_P1               = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_P1") );
    TGraph *g_Proton_number_raw_High_P2               = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_P2") );
    TGraph *g_Proton_number_raw_High_P3               = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_P3") );
    TGraph *g_Proton_number_raw_High_P4               = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_P4") );
    TGraph *g_Proton_number_raw_High_P5               = new TGraph( (TH1D*)highfile_MixedPattern->Get("h_Proton_number_Raw_High_P5") );

    TH1D *h_Antiproton_number_unfolded_High_P0        = (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P0");
    TH1D *h_Antiproton_number_unfolded_High_P0VGGNN   = (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P0VGGNN");
    TH1D *h_Antiproton_number_unfolded_High_P1        = (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P1");
    TH1D *h_Antiproton_number_unfolded_High_P2        = (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P2");
    TH1D *h_Antiproton_number_unfolded_High_P3        = (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P3");
    TH1D *h_Antiproton_number_unfolded_High_P4        = (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P4");
    TH1D *h_Antiproton_number_unfolded_High_P5        = (TH1D*)highfile_MixedPattern->Get("h_Antiproton_number_unfolded_High_P5");

    // EffectiveAcceptance
    TGraphErrors *Effective_Acceptance_High_B1042_P0 = (TGraphErrors*)highfile_MixedPattern->Get("Effective_Acceptance_P0");
    TGraphErrors *Effective_Acceptance_High_B1042_P1 = (TGraphErrors*)highfile_MixedPattern->Get("Effective_Acceptance_P1");
    TGraphErrors *Effective_Acceptance_High_B1042_P2 = (TGraphErrors*)highfile_MixedPattern->Get("Effective_Acceptance_P2");
    TGraphErrors *Effective_Acceptance_High_B1042_P4 = (TGraphErrors*)highfile_MixedPattern->Get("Effective_Acceptance_P4");
    TGraphErrors *Effective_Acceptance_High_B1042    = (TGraphErrors*)highfile_MixedPattern->Get("Effective_Acceptance_mixedPattern");
    // Error
    TGraph *g_StatisticalError_High = (TGraph*)highfile_MixedPattern->Get("g_StatisticalError_mixedPattern");
    TH1D *h_Statistic_Error_Relative_High = (TH1D*)highfile_MixedPattern->Get("h_Statistic_Error_Relative_mixedPattern");
    TGraph *g_System_CC_High = (TGraph*)highfile_MixedPattern->Get("g_System_CC_High_mixedPattern");
    // Fix and Smooth g_System_CC_High (Due statistics in first few bins, fix trival values for plots)
    g_System_CC_High->SetPoint(0, g_System_CC_High->GetX()[0], g_System_CC_High->GetY()[7]);
    g_System_CC_High->SetPoint(1, g_System_CC_High->GetX()[1], g_System_CC_High->GetY()[7]);
    g_System_CC_High->SetPoint(2, g_System_CC_High->GetX()[2], g_System_CC_High->GetY()[7]);
    g_System_CC_High->SetPoint(3, g_System_CC_High->GetX()[3], g_System_CC_High->GetY()[7]);
    g_System_CC_High->SetPoint(4, g_System_CC_High->GetX()[4], g_System_CC_High->GetY()[7]);
    g_System_CC_High->SetPoint(5, g_System_CC_High->GetX()[5], g_System_CC_High->GetY()[7]);
    g_System_CC_High->SetPoint(6, g_System_CC_High->GetX()[6], g_System_CC_High->GetY()[7]);
    g_System_CC_High->SetPoint(9, g_System_CC_High->GetX()[9], g_System_CC_High->GetY()[8]);
    g_System_CC_High->SetPoint(11, g_System_CC_High->GetX()[11], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(12, g_System_CC_High->GetX()[12], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(13, g_System_CC_High->GetX()[13], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(14, g_System_CC_High->GetX()[14], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(15, g_System_CC_High->GetX()[15], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(16, g_System_CC_High->GetX()[16], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(17, g_System_CC_High->GetX()[17], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(18, g_System_CC_High->GetX()[18], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(19, g_System_CC_High->GetX()[19], g_System_CC_High->GetY()[10]);
    g_System_CC_High->SetPoint(20, g_System_CC_High->GetX()[20], g_System_CC_High->GetY()[10]);

    //// Load Reference pbar Number
    TGraph *g_Published_pbarNumber;
    TGraph *g_PhysicsReport_pbarNumber;
    tie(g_Published_pbarNumber, g_PhysicsReport_pbarNumber) = Load_Reference_pbar_Number();


    //// Load EffectiveAcceptance correction
    TGraphErrors *Effective_Acceptance_ratio_B1220;
    TGraphErrors *Effective_Acceptance_ratio_B1220_minus10;
    TGraphErrors *Effective_Acceptance_ratio_B1220_plus10;
    TGraph *g_SysUncertaintyRel_ACCratio_B1220;
    TH1D h_SysUncertaintyRel_ACCratio;
    tie(Effective_Acceptance_ratio_B1220, Effective_Acceptance_ratio_B1220_minus10, Effective_Acceptance_ratio_B1220_plus10, g_SysUncertaintyRel_ACCratio_B1220, h_SysUncertaintyRel_ACCratio) = Load_EffectiveAcceptance_Correction(PhysicsReportBinEdge);


    //// General define the overlap range.
    int LowRangeRemoveAtEnd            = 3;
    int IntermediateRangeRemoveAtBegin = 4;
    int IntermediateRangeRemoveAtEnd   = 2;
    int HighRangeRemovedAtBegin        = 1;


    //// Calculate "Total error" and "Total Systematic error" in This analysis
    // Prepare: (A Simple Fix to avoid small y axis limit)
    h_StatisticalRelError_Low         ->SetMaximum(100);
    h_StatisticalRelError_Intermediate->SetMaximum(100);
    h_Statistic_Error_Relative_High   ->SetMaximum(100);

    // Total Systematic Error (absolute) (TGraph)
    TGraph *g_TotalSysError_Low          = new TGraph(g_StatisticalError_Low->GetN());           //(To be finished...)
    TGraph *g_TotalSysError_Intermediate = new TGraph(g_StatisticError_Intermediate->GetN()); 
    TGraph *g_TotalSysError_High         = new TGraph(g_StatisticalError_High->GetN());
    TGraph *g_System_ACC_Low             = new TGraph(g_StatisticalError_Low->GetN());
    TGraph *g_System_ACC_Intermediate    = new TGraph(g_StatisticError_Intermediate->GetN());
    TGraph *g_System_ACC_High            = new TGraph(g_StatisticalError_High->GetN());

    // Total Systematic Error (relative) (TGraph)
    TGraph *g_TotalSysRelError_Low          = new TGraph(g_StatisticalError_Low->GetN());        //(To be finished...)
    TGraph *g_TotalSysRelError_Intermediate = new TGraph(g_StatisticError_Intermediate->GetN()); 
    TGraph *g_TotalSysRelError_High         = new TGraph(g_StatisticalError_High->GetN());
    TGraph *g_System_CC_relative_High       = new TGraph(g_StatisticalError_High->GetN());
    // Total Systematic Error (relative) (TH1D)
    TH1D *h_TotalSysRelError_Low          = new TH1D(*h_StatisticalRelError_Low);                //(To be finished...)
    TH1D *h_TotalSysRelError_Intermediate = new TH1D(*h_StatisticalRelError_Intermediate); 
    TH1D *h_TotalSysRelError_High         = new TH1D(*h_Statistic_Error_Relative_High);    
    TH1D *h_System_CC_relative_High       = new TH1D(*h_Statistic_Error_Relative_High);
    TH1D *h_SystematicError_Shape_Low     = new TH1D(*h_SystematicError_Low);

    // Total Error (absolute)
    TGraph *g_TotalError_Low          = new TGraph(g_StatisticalError_Low->GetN());
    TGraph *g_TotalError_Intermediate = new TGraph(g_StatisticError_Intermediate->GetN());
    TGraph *g_TotalError_High         = new TGraph(g_StatisticalError_High->GetN());

    // Total Error (relative) (TGraph)
    TGraph *g_TotalRelError_Low          = new TGraph(g_StatisticalError_Low->GetN());
    TGraph *g_TotalRelError_Intermediate = new TGraph(g_StatisticError_Intermediate->GetN());
    TGraph *g_TotalRelError_High         = new TGraph(g_StatisticalError_High->GetN());
    // Total Error (relative) (TH1D)
    TH1D *h_TotalRelError_Low          = new TH1D(*h_StatisticalRelError_Low);
    TH1D *h_TotalRelError_Intermediate = new TH1D(*h_StatisticalRelError_Intermediate);
    TH1D *h_TotalRelError_High         = new TH1D(*h_Statistic_Error_Relative_High);


    //// Calculate "Total error" and "Total Systematic error" in This analysis
    Calculate_TotalError_TotalSystematicError_LowRange(g_StatisticalError_Low, g_SysUncertaintyRel_ACCratio_B1220, g_System_ACC_Low, g_SystematicError_Shape_Low, g_LowResult, g_TotalSysError_Low, g_TotalSysRelError_Low, h_TotalSysRelError_Low, g_TotalError_Low, g_TotalRelError_Low, h_TotalRelError_Low, h_SystematicError_Shape_Low);

    Calculate_TotalError_TotalSystematicError_IntermediateRange(g_StatisticError_Intermediate, g_SysUncertaintyRel_ACCratio_B1220, g_StatisticalError_Low, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, g_System_ACC_Intermediate, g_IntermediateResult, g_TotalSysError_Intermediate, g_TotalSysRelError_Intermediate, h_TotalSysRelError_Intermediate, g_TotalError_Intermediate, g_TotalRelError_Intermediate, h_TotalRelError_Intermediate, g_SystematicError_TRD_Intermediate, h_SystematicRelativeError_TRD);

    Calculate_TotalError_TotalSystematicError_HighRange(g_StatisticalError_High, g_SysUncertaintyRel_ACCratio_B1220, g_StatisticalError_Low, g_StatisticError_Intermediate, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin, g_System_ACC_High, g_TotalSysError_High, g_HighResult, g_System_CC_High, g_TotalSysRelError_High, g_System_CC_relative_High, h_TotalSysRelError_High, h_System_CC_relative_High, g_TotalError_High, g_TotalRelError_High, h_TotalRelError_High);


    //// Binning shift.
    double binshift = 1.0; //0.95
    BinningShift(binshift, g_LowResult, g_IntermediateResult, g_HighResult);


    //// Plot: Part1 Plot the Pbar Ratio with Statistics Error Only
    Plot_Ratio_OnlyStaErr(g_HighResult, g_LowResult, g_IntermediateResult, issversionname, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin);
    Plot_RawRatioAndUnfoldedRatio_Compare(issversion, g_Raw_Unfold_Compare_P0VGGNN, g_Raw_Unfold_Compare_P1, g_Raw_Unfold_Compare_P2, g_Raw_Unfold_Compare_P4, g_RawLowResult, g_LowResult, g_RawIntermediateResult_UnfoldingIterationMore, g_IntermediateResult_UnfoldingIterationMore, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin);
    

    //// Make high rigidity result with scaled acceptance error
    TGraphErrors *g_LowResult_ScaledAcceptanceError          = new TGraphErrors(g_LowResult->GetN());
    TGraphErrors *g_IntermediateResult_ScaledAcceptanceError = new TGraphErrors(g_IntermediateResult->GetN());
    TGraphErrors *g_HighResult_ScaledAcceptanceError         = new TGraphErrors(g_HighResult->GetN());
    MakeResultWithScaledAcceptanceError(g_LowResult_ScaledAcceptanceError, g_IntermediateResult_ScaledAcceptanceError, g_HighResult_ScaledAcceptanceError, g_StatisticalError_Low, g_StatisticError_Intermediate, g_StatisticalError_High, g_SysUncertaintyRel_ACCratio_B1220, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin, g_LowResult, g_IntermediateResult, g_HighResult, g_System_CC_High, g_System_ACC_Low, g_System_ACC_Intermediate, g_System_ACC_High, g_SystematicError_Shape_Low, g_SystematicError_TRD_Intermediate, Scaler);


    //// Reset Error from Statistical Error to total error: Low, Intermediate, High.
    Reset_Error_From_Statistical_To_Total(g_LowResult   , g_IntermediateResult, g_HighResult, g_TotalError_Low, g_TotalError_Intermediate, g_TotalError_High);

    //// Deal with overlap range for: Fill the Relative Error Only for Historams.
    Reset_Error_For_OverlapedRange(h_StatisticalRelError_Low, h_StatisticalRelError_Intermediate, h_TotalRelError_Low, h_TotalRelError_Intermediate, 
    h_TotalSysRelError_Low, h_TotalSysRelError_Intermediate, h_Statistic_Error_Relative_High, h_TotalRelError_High, h_TotalSysRelError_High, h_System_CC_relative_High, h_SystematicError_Shape_Low, h_SystematicRelativeError_TRD, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin);


    //// Plot: Part2
    
    Plot_Antiproton_and_Proton_Numbers_LowRange(g_PhysicsReport_pbarNumber, g_Published_pbarNumber, g_antiproton_number_raw_Low, g_Proton_number_raw_Low, issversionname); // Input can be raw or unfolded numbers.
    Plot_Antiproton_and_Proton_Numbers_IntermediateRange(g_PhysicsReport_pbarNumber, g_Published_pbarNumber, g_Antiproton_number_raw_Intermediate, g_Proton_number_raw_Intermediate, issversionname); // Input can be raw or unfolded numbers.
    Plot_Antiproton_and_Proton_Numbers_HighRange(g_PhysicsReport_pbarNumber, g_Published_pbarNumber, g_Antiproton_number_unfolded_High,issversionname); // Input can be raw or unfolded numbers.
    
    
    // Remove overlaped range for: Antiproton and Proton Numbers and Acceptance.
    Remove_OverlapedRange_Antiproton_Proton_Numbers_and_Acceptance(
    g_Antiproton_number_unfolded_Low         , g_Proton_number_unfolded_Low         , Effective_Acceptance_Low_B1042, 
    g_Antiproton_number_unfolded_Intermediate, g_Proton_number_unfolded_Intermediate, Effective_Acceptance_Intermediate_B1042, 
    g_Antiproton_number_unfolded_High        , g_Proton_number_unfolded_High        , Effective_Acceptance_Intermediate_inHigh_B1042, 
    g_Antiproton_number_unfolded_High_P0VGGNN, g_Antiproton_number_unfolded_High_P0 , g_Antiproton_number_unfolded_High_P1, g_Antiproton_number_unfolded_High_P2, g_Antiproton_number_unfolded_High_P4,
    g_Proton_number_unfolded_High_P0VGGNN    , g_Proton_number_unfolded_High_P0     , g_Proton_number_unfolded_High_P1    , g_Proton_number_unfolded_High_P2    , g_Proton_number_unfolded_High_P4, 
    g_Antiproton_number_raw_High_P0VGGNN     , g_Antiproton_number_raw_High_P0      , g_Antiproton_number_raw_High_P1     , g_Antiproton_number_raw_High_P2     , g_Antiproton_number_raw_High_P4, 
    g_Proton_number_raw_High_P0VGGNN         , g_Proton_number_raw_High_P0          , g_Proton_number_raw_High_P1         , g_Proton_number_raw_High_P2         , g_Proton_number_raw_High_P4,
    LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin, Effective_Acceptance_High_B1042, 
    g_antiproton_number_raw_Low              , g_Antiproton_number_raw_Intermediate , g_Antiproton_number_raw_High,
    g_Proton_number_raw_Low                  , g_Proton_number_raw_Intermediate     , g_Proton_number_raw_High);

    //// Apply Relative Uncertainty From B1220 version To B1042 version 
    TGraphErrors *Effective_Acceptance_ratio_B1042_minus10_AllRange;
    TGraphErrors *Effective_Acceptance_ratio_B1042_plus10_AllRange;
    TGraphAsymmErrors *Effective_Acceptance_ratio_B1042_AllRange;
    tie(Effective_Acceptance_ratio_B1042_minus10_AllRange, Effective_Acceptance_ratio_B1042_plus10_AllRange, Effective_Acceptance_ratio_B1042_AllRange) = ApplyRelativeUncertaintyToB1042(Effective_Acceptance_ratio_B1220, Effective_Acceptance_ratio_B1220_minus10, Effective_Acceptance_ratio_B1220_plus10, Effective_Acceptance_High_B1042, Effective_Acceptance_Low_B1042, Effective_Acceptance_Intermediate_B1042, Effective_Acceptance_Intermediate_inHigh_B1042); 
   
    
    //// Plot: Part3
    Plot_Antiproton_and_Proton_Numbers_OverRigidityBinWidth(g_PhysicsReport_pbarNumber, g_Published_pbarNumber, 
    g_Antiproton_number_unfolded_Low         , g_Proton_number_unfolded_Low, 
    g_Antiproton_number_unfolded_Intermediate, g_Proton_number_unfolded_Intermediate, 
    g_Antiproton_number_unfolded_High        , g_Antiproton_number_unfolded_High_P0, g_Antiproton_number_unfolded_High_P0VGGNN, g_Antiproton_number_unfolded_High_P1, g_Antiproton_number_unfolded_High_P2, g_Antiproton_number_unfolded_High_P3, g_Antiproton_number_unfolded_High_P4, g_Antiproton_number_unfolded_High_P5, 
    g_Proton_number_unfolded_High            , g_Proton_number_unfolded_High_P0         , g_Proton_number_unfolded_High_P0VGGNN, g_Proton_number_unfolded_High_P1   , g_Proton_number_unfolded_High_P2    , g_Proton_number_unfolded_High_P3    , g_Proton_number_unfolded_High_P4    , g_Proton_number_unfolded_High_P5,
    g_antiproton_number_raw_Low         , g_Proton_number_raw_Low, 
    g_Antiproton_number_raw_Intermediate, g_Proton_number_raw_Intermediate, 
    g_Antiproton_number_raw_High        , g_Antiproton_number_raw_High_P0, g_Antiproton_number_raw_High_P0VGGNN, g_Antiproton_number_raw_High_P1, g_Antiproton_number_raw_High_P2, g_Antiproton_number_raw_High_P3, g_Antiproton_number_raw_High_P4, g_Antiproton_number_raw_High_P5, 
    g_Proton_number_raw_High            , g_Proton_number_raw_High_P0    , g_Proton_number_raw_High_P0VGGNN    , g_Proton_number_raw_High_P1    , g_Proton_number_raw_High_P2    , g_Proton_number_raw_High_P3    , g_Proton_number_raw_High_P4    , g_Proton_number_raw_High_P5,
    h_Antiproton_number_unfolded_Low, h_Antiproton_number_unfolded_Intermediate, h_Antiproton_number_unfolded_High_P0, issversionname, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin); //Also print numbers.    
    
    Plot_Antiproton_and_Proton_Numbers_Ratio(g_PhysicsReport_pbarNumber, g_Published_pbarNumber, g_Antiproton_number_unfolded_Low, g_Proton_number_unfolded_Low, g_Antiproton_number_unfolded_Intermediate, g_Proton_number_unfolded_Intermediate, g_Antiproton_number_unfolded_High, g_Antiproton_number_unfolded_High_P0VGGNN, g_Antiproton_number_unfolded_High_P1, g_Antiproton_number_unfolded_High_P2, g_Antiproton_number_unfolded_High_P3, g_Antiproton_number_unfolded_High_P4, g_Antiproton_number_unfolded_High_P5, g_Proton_number_unfolded_High, issversionname, LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin);
    
    Plot_AntiprotonNumbersAllPattern(g_PhysicsReport_pbarNumber, g_Published_pbarNumber, g_Antiproton_number_unfolded_Low, g_Proton_number_unfolded_Low, g_Antiproton_number_unfolded_Intermediate, g_Proton_number_unfolded_Intermediate, g_Antiproton_number_unfolded_High, g_Antiproton_number_unfolded_High_P0VGGNN, g_Antiproton_number_unfolded_High_P1, g_Antiproton_number_unfolded_High_P2, g_Antiproton_number_unfolded_High_P3, g_Antiproton_number_unfolded_High_P4, g_Antiproton_number_unfolded_High_P5, g_Proton_number_unfolded_High, issversionname);

    Plot_AntiprotonNumbersOverRigidityBinWidth_High(g_Antiproton_number_raw_High_P0VGGNN, g_Antiproton_number_raw_High_P1, g_Antiproton_number_raw_High_P2, g_Antiproton_number_raw_High_P3, g_Antiproton_number_raw_High_P4, g_Antiproton_number_raw_High_P5, h_Antiproton_number_unfolded_High_P0VGGNN, h_Antiproton_number_unfolded_High_P1, h_Antiproton_number_unfolded_High_P2, h_Antiproton_number_unfolded_High_P3, h_Antiproton_number_unfolded_High_P4, h_Antiproton_number_unfolded_High_P5, g_Proton_number_raw_High_P0VGGNN, g_Proton_number_raw_High_P1, g_Proton_number_raw_High_P2, g_Proton_number_raw_High_P3, g_Proton_number_raw_High_P4, g_Proton_number_raw_High_P5, issversionname, HighRangeRemovedAtBegin); // Input can be raw or unfolded numbers. 

    
    Plot_AntiprotonNumbersToRrferencesRatio_InHighRange(g_PhysicsReport_pbarNumber, g_Published_pbarNumber, g_Antiproton_number_unfolded_Low, g_Proton_number_unfolded_Low, g_Antiproton_number_unfolded_Intermediate, g_Proton_number_unfolded_Intermediate, g_Antiproton_number_unfolded_High, g_Antiproton_number_unfolded_High_P0VGGNN, g_Antiproton_number_unfolded_High_P1, g_Antiproton_number_unfolded_High_P2, g_Antiproton_number_unfolded_High_P4, g_Proton_number_unfolded_High, issversionname);    

    Plot_EffectiveAcceptanceInThreeRanges(Effective_Acceptance_High_B1042_P0, Effective_Acceptance_High_B1042_P1, Effective_Acceptance_High_B1042_P2, Effective_Acceptance_High_B1042_P4, Effective_Acceptance_High_B1042, Effective_Acceptance_Low_B1042, Effective_Acceptance_Intermediate_B1042, Effective_Acceptance_Intermediate_inHigh_B1042);
    
    //// Plot: Part4. Overlapped Ratio And Remove the first two points in interrmediate range for illustration. 
    Plot_Ratio_ThisAnalysisOnly_Overlapped(gPublishedRatio, g_LowResult, g_IntermediateResult, g_HighResult, issversionname); // (Overlapped range Included)
     
    //// Remove overlaped range For: 1. Ratio (Except the first two points in the intermediate range, which have been removed in "Plot_Ratio_ThisAnalysisOnly_Overlapped") 2. Error
    Remove_OverlapedRange_Ratio(g_LowResult, g_IntermediateResult, g_HighResult, g_LowResult_ScaledAcceptanceError, g_IntermediateResult_ScaledAcceptanceError, g_HighResult_ScaledAcceptanceError, 
    g_System_ACC_Low, g_System_ACC_Intermediate, g_System_ACC_High, g_TotalSysError_Low, g_TotalSysError_Intermediate, g_TotalSysError_High, 
    g_StatisticalError_Low, g_StatisticError_Intermediate, g_StatisticalError_High, 
    LowRangeRemoveAtEnd, IntermediateRangeRemoveAtBegin, IntermediateRangeRemoveAtEnd, HighRangeRemovedAtBegin);

    //// Plot: Part5 Continue to Plot 
    Plot_Ratio_CompareWithReference_NotOverlapped(gPublishedRatio, gPhyrReortRatio, g_HighResult, g_LowResult, g_IntermediateResult, issversionname, issversion); // (Overlapped range NOT Included)
    
    Plot_Ratio_AllTrackerPatterns_CompareWithReference_NotOverlapped(gPublishedRatio, gPhyrReortRatio, g_HighResult, g_LowResult, g_IntermediateResult, g_HighResult_P0VGGNN, g_HighResult_P1, g_HighResult_P2, g_HighResult_P3, g_HighResult_P4, g_HighResult_P5, issversionname); 
    Plot_Ratio_CompareBetweenReferences(gPhyrReortRatio, gPublishedRatio);
    //Plot_Ratio_ThisAnalysisOnly(g_HighResult, g_LowResult, g_IntermediateResult, issversionname, Scaler);  // Also do the linear fit for the high range result.
    Plot_Ratio_ThisAnalysisOnly(g_HighResult_ScaledAcceptanceError, g_LowResult_ScaledAcceptanceError, g_IntermediateResult_ScaledAcceptanceError, issversionname, Scaler);  // Also do the linear fit for the high range result.

    Plot_Statictical_Error_Absolute(g_PublishedRatioStatisticErrorPRL, g_StatisticalError_Low, g_StatisticError_Intermediate, g_StatisticalError_High, g_PhysicsReportRatioStatisticError, issversionname);
    Plot_Statictical_Error_Relative(h_PhysicsReportStatisticRatioRelativeError, h_PublishedRatioStatisticRelativeErrorPRL, h_StatisticalRelError_Low, h_StatisticalRelError_Intermediate, h_Statistic_Error_Relative_High, PublishedPRLBinEdge, issversion, issversionname);
    
    Plot_Systematical_Error_Absolute_Acc     (g_System_ACC_Low, g_System_ACC_Intermediate, g_System_ACC_High);
    Plot_Systematical_Error_Absolute_FitRange(g_SystematicError_TRD_Intermediate, g_SystematicError_Shape_Low, issversionname);
    Plot_Systematical_Error_Absolute_CC      (g_System_CC_High, issversionname);
    Plot_Systematical_Error_Absolute         (g_TotalSysError_High, g_System_CC_High, g_TotalSysError_Intermediate, g_SystematicError_TRD_Intermediate, g_TotalSysError_Low, g_SystematicError_Shape_Low, g_System_ACC_Low, g_System_ACC_Intermediate, g_System_ACC_High, g_PublishedRatioSystematicErrorPRL, g_PhysicsReportRatioSystematicError, issversionname);

    Plot_Systematical_Error_Relative_CC      (h_System_CC_relative_High   , issversionname);
    Plot_Systematical_Error_Relative_ACC     (h_SysUncertaintyRel_ACCratio, issversionname);
    Plot_Systematical_Error_Relative_FitRange(h_SystematicError_Shape_Low, h_SystematicRelativeError_TRD, g_SystematicError_Shape_Low, issversionname);
 
    Plot_Relatvie_Systematical_Error(h_SysUncertaintyRel_ACCratio, h_PublishedRatioSystematicRelativeErrorPRL, h_PhysicsReportSystematicRatioRelativeError, h_System_CC_relative_High, h_TotalSysRelError_High, h_TotalSysRelError_Intermediate, h_TotalSysRelError_Low, PublishedPRLBinEdge, issversionname);

    Plot_StaSysRelErrCompare        (h_PublishedRatioSystematicRelativeErrorPRL, h_PhysicsReportSystematicRatioRelativeError, h_TotalSysRelError_High, h_TotalSysRelError_Intermediate, h_TotalSysRelError_Low, h_PhysicsReportStatisticRatioRelativeError, h_PublishedRatioStatisticRelativeErrorPRL, h_StatisticalRelError_Low, h_StatisticalRelError_Intermediate, h_Statistic_Error_Relative_High, PublishedPRLBinEdge, issversionname, issversion);

    Plot_Total_Error_Absolute                 (g_PublishedRatioToralError, g_TotalError_Low, g_TotalError_Intermediate, g_TotalError_High, g_PhysicsReportRatioError, issversionname); //(FIXME: Deal with overlap range for g_TotalError_Low, Intermediate, High )
    Plot_Total_Relative_Error_CompareReference(h_PhysicsReportRatioRelativeError, h_PublishedRatioRelativeErrorPRL, h_TotalRelError_Low, h_TotalRelError_Intermediate, h_TotalRelError_High, PublishedPRLBinEdge, issversion, issversionname);
    Plot_Total_Relative_Error_ThisAnalysis    (h_PhysicsReportRatioRelativeError, h_TotalRelError_High, h_TotalRelError_Intermediate, h_TotalRelError_Low, h_TotalSysRelError_Low, h_TotalSysRelError_Intermediate, h_TotalSysRelError_High, h_StatisticalRelError_Low, h_StatisticalRelError_Intermediate, h_Statistic_Error_Relative_High, issversionname);

    Plot_Effective_Acceptance_Ratios_With_Uncertainty_B1220(Effective_Acceptance_ratio_B1220, Effective_Acceptance_ratio_B1220_minus10, Effective_Acceptance_ratio_B1220_plus10); 
    Plot_Effective_Acceptance_Ratios_With_Uncertainty_B1042(Effective_Acceptance_High_B1042 , Effective_Acceptance_Low_B1042, Effective_Acceptance_Intermediate_B1042, Effective_Acceptance_Intermediate_inHigh_B1042, Effective_Acceptance_ratio_B1042_minus10_AllRange, Effective_Acceptance_ratio_B1042_plus10_AllRange, Effective_Acceptance_ratio_B1042_AllRange);  // Also remove the y error bars for Inllustrtation
    

    /*    
    //// Save in root file to check
    TFile checkrootfile( (string("CheckResult.root")).c_str(),"RECREATE");

    g_LowResult         ->Write("g_LowResul");
    g_IntermediateResult->Write("g_Intermediate");
    g_HighResult        ->Write("g_HighResult");

    g_Antiproton_number_unfolded_Low         ->Write("g_Antiproton_number_unfolded_Low");
    g_Antiproton_number_unfolded_Intermediate->Write("g_Antiproton_number_unfolded_Intermediate");
    g_Antiproton_number_unfolded_High        ->Write("g_Antiproton_number_unfolded_High");

    g_TotalError_Low          ->Write("g_TotalError_Low");
    g_TotalError_Intermediate ->Write("g_TotalError_Intermediate");
    g_TotalError_High         ->Write("g_TotalError_High");

    checkrootfile.Close();
    */ 
}

