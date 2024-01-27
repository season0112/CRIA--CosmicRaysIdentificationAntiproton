#include "Plot_Acceptance.hh"

void Plot_Acceptance(){

    // Low & Intrmediate
    TFile *f_eff_antiproton_low = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root");
    TFile *f_eff_proton_low     = new TFile("/hpcwork/jara0052/sichen/AntiprotonLowEnergy_v11.0/totalall/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root");
    TGraphAsymmErrors *Acceptance_antiproton_low = (TGraphAsymmErrors*)f_eff_antiproton_low->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_proton_low     = (TGraphAsymmErrors*)f_eff_proton_low->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");

    // High
    TFile *f_Anti_P0 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_Pattern") + std::string("0") + std::string("_525version.root")).c_str() );
    TFile *f_Pro_P0  = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_Pattern") + std::string("0") + std::string("_525version.root")).c_str() );
    TFile *f_Anti_P1 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_Pattern") + std::string("1") + std::string("_525version.root")).c_str() );
    TFile *f_Pro_P1  = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_Pattern") + std::string("1") + std::string("_525version.root")).c_str() );
    TFile *f_Anti_P2 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_Pattern") + std::string("2") + std::string("_525version.root")).c_str() );
    TFile *f_Pro_P2  = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_Pattern") + std::string("2") + std::string("_525version.root")).c_str() );
    TFile *f_Anti_P4 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_Pattern") + std::string("4") + std::string("_525version.root")).c_str() );
    TFile *f_Pro_P4  = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_Pattern") + std::string("4") + std::string("_525version.root")).c_str() );
    TFile *f_Anti_PMinus1 = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all_Pattern") + std::string("-1") + std::string("_525version.root")).c_str() );
    TFile *f_Pro_PMinus1  = new TFile( (std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all_Pattern") + std::string("-1") + std::string("_525version.root")).c_str() );
    TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton_high_P0      = (TGraphAsymmErrors*)f_Anti_P0->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Proton_high_P0          = (TGraphAsymmErrors*)f_Pro_P0 ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton_high_P1      = (TGraphAsymmErrors*)f_Anti_P1->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Proton_high_P1          = (TGraphAsymmErrors*)f_Pro_P1 ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton_high_P2      = (TGraphAsymmErrors*)f_Anti_P2->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Proton_high_P2          = (TGraphAsymmErrors*)f_Pro_P2 ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton_high_P4      = (TGraphAsymmErrors*)f_Anti_P4->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Proton_high_P4          = (TGraphAsymmErrors*)f_Pro_P4 ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Antiproton_high_PMinus1 = (TGraphAsymmErrors*)f_Anti_PMinus1->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *g_EffectiveAcceptance_Proton_high_PMinus1     = (TGraphAsymmErrors*)f_Pro_PMinus1 ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");

    cout<< Acceptance_proton_low->GetN()     << endl;
    cout<< Acceptance_antiproton_low->GetN() << endl;
    cout<< g_EffectiveAcceptance_Antiproton_high_P0->GetN() << endl;
    cout<< g_EffectiveAcceptance_Proton_high_P0->GetN()     << endl;

    for (int i=0; i<35; i++){
        Acceptance_proton_low    ->RemovePoint(28);
        Acceptance_antiproton_low->RemovePoint(28);
    }
    for (int i=0; i<28; i++){
        g_EffectiveAcceptance_Antiproton_high_P0->RemovePoint(0);
        g_EffectiveAcceptance_Proton_high_P0    ->RemovePoint(0);
        g_EffectiveAcceptance_Antiproton_high_P1->RemovePoint(0);
        g_EffectiveAcceptance_Proton_high_P1    ->RemovePoint(0);
        g_EffectiveAcceptance_Antiproton_high_P2->RemovePoint(0);
        g_EffectiveAcceptance_Proton_high_P2    ->RemovePoint(0);
        g_EffectiveAcceptance_Antiproton_high_P4->RemovePoint(0);
        g_EffectiveAcceptance_Proton_high_P4    ->RemovePoint(0);
    }

    // FIXME: for First or Last point in plot, due to rigidity edge acceptance goes down.
    Acceptance_proton_low                   ->SetPoint(26, Acceptance_proton_low->GetX()[26]                   , Acceptance_proton_low->GetY()[26]*1.02);
    Acceptance_antiproton_low               ->SetPoint(26, Acceptance_antiproton_low->GetX()[26]               , Acceptance_antiproton_low->GetY()[26]*1.02);    
    Acceptance_proton_low                   ->SetPoint(27, Acceptance_proton_low->GetX()[27]                   , Acceptance_proton_low->GetY()[27]*1.05);
    Acceptance_antiproton_low               ->SetPoint(27, Acceptance_antiproton_low->GetX()[27]               , Acceptance_antiproton_low->GetY()[27]*1.05);

    g_EffectiveAcceptance_Antiproton_high_P0->SetPoint(1, g_EffectiveAcceptance_Antiproton_high_P0->GetX()[1], g_EffectiveAcceptance_Antiproton_high_P0->GetY()[2]*1.00);
    g_EffectiveAcceptance_Proton_high_P0    ->SetPoint(1, g_EffectiveAcceptance_Proton_high_P0->GetX()[1]    , g_EffectiveAcceptance_Proton_high_P0->GetY()[2]*1.00);
    g_EffectiveAcceptance_Antiproton_high_P1->SetPoint(1, g_EffectiveAcceptance_Antiproton_high_P1->GetX()[1], g_EffectiveAcceptance_Antiproton_high_P1->GetY()[2]*1.00);
    g_EffectiveAcceptance_Proton_high_P1    ->SetPoint(1, g_EffectiveAcceptance_Proton_high_P1->GetX()[1]    , g_EffectiveAcceptance_Proton_high_P1->GetY()[2]*1.00);
    g_EffectiveAcceptance_Antiproton_high_P2->SetPoint(1, g_EffectiveAcceptance_Antiproton_high_P2->GetX()[1], g_EffectiveAcceptance_Antiproton_high_P2->GetY()[2]*1.00);
    g_EffectiveAcceptance_Proton_high_P2    ->SetPoint(1, g_EffectiveAcceptance_Proton_high_P2->GetX()[1]    , g_EffectiveAcceptance_Proton_high_P2->GetY()[2]*1.00);
    g_EffectiveAcceptance_Antiproton_high_P4->SetPoint(1, g_EffectiveAcceptance_Antiproton_high_P4->GetX()[1], g_EffectiveAcceptance_Antiproton_high_P4->GetY()[2]*1.00);
    g_EffectiveAcceptance_Proton_high_P4    ->SetPoint(1, g_EffectiveAcceptance_Proton_high_P4->GetX()[1]    , g_EffectiveAcceptance_Proton_high_P4->GetY()[2]*1.00);

    g_EffectiveAcceptance_Antiproton_high_P0->SetPoint(0, g_EffectiveAcceptance_Antiproton_high_P0->GetX()[0], g_EffectiveAcceptance_Antiproton_high_P0->GetY()[1]*1.00);
    g_EffectiveAcceptance_Proton_high_P0    ->SetPoint(0, g_EffectiveAcceptance_Proton_high_P0->GetX()[0]    , g_EffectiveAcceptance_Proton_high_P0->GetY()[1]*1.00);
    g_EffectiveAcceptance_Antiproton_high_P1->SetPoint(0, g_EffectiveAcceptance_Antiproton_high_P1->GetX()[0], g_EffectiveAcceptance_Antiproton_high_P1->GetY()[1]*1.00);
    g_EffectiveAcceptance_Proton_high_P1    ->SetPoint(0, g_EffectiveAcceptance_Proton_high_P1->GetX()[0]    , g_EffectiveAcceptance_Proton_high_P1->GetY()[1]*1.00);
    g_EffectiveAcceptance_Antiproton_high_P2->SetPoint(0, g_EffectiveAcceptance_Antiproton_high_P2->GetX()[0], g_EffectiveAcceptance_Antiproton_high_P2->GetY()[1]*1.00);
    g_EffectiveAcceptance_Proton_high_P2    ->SetPoint(0, g_EffectiveAcceptance_Proton_high_P2->GetX()[0]    , g_EffectiveAcceptance_Proton_high_P2->GetY()[1]*1.00);
    g_EffectiveAcceptance_Antiproton_high_P4->SetPoint(0, g_EffectiveAcceptance_Antiproton_high_P4->GetX()[0], g_EffectiveAcceptance_Antiproton_high_P4->GetY()[1]*1.00);
    g_EffectiveAcceptance_Proton_high_P4    ->SetPoint(0, g_EffectiveAcceptance_Proton_high_P4->GetX()[0]    , g_EffectiveAcceptance_Proton_high_P4->GetY()[1]*1.00);


    //Plot
    Plot_EffectiveAcceptance(Acceptance_proton_low, g_EffectiveAcceptance_Proton_high_P0, g_EffectiveAcceptance_Proton_high_P1, g_EffectiveAcceptance_Proton_high_P2, g_EffectiveAcceptance_Proton_high_P4, std::string("Proton"));
    Plot_EffectiveAcceptance(Acceptance_antiproton_low, g_EffectiveAcceptance_Antiproton_high_P0, g_EffectiveAcceptance_Antiproton_high_P1, g_EffectiveAcceptance_Antiproton_high_P2, g_EffectiveAcceptance_Antiproton_high_P4, std::string("Antiproton"));



}
