
void Plot_Trigger(){

    string PatternName = "_Pattern0";
    string MCName = "B1042_antipr.pl1.1800_7.6_all";
    //string MCName = "B1042_pr.pl1.1800_7.6_all";

    //// Load Trigger
    TFile *f_Trigger   = new TFile( (getenv("HPCHIGHENERGYDATADIR") + string("/TriggerEff_") + MCName + string("_525version") + PatternName + string(".root")) .c_str()) ;

    TH1F *h_TriggerEff              = (TH1F*)f_Trigger->Get("TriggerEff");
    TH1F *h_TriggerEff_noprescaling = (TH1F*)f_Trigger->Get("TriggerEff_noprescaling");

    h_TriggerEff              ->SetTitle("");
    h_TriggerEff_noprescaling ->SetTitle("");

    /*
    TCanvas c1("c1","c1",1000,500);
    TPad *p1 = new TPad("p1", "p1", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p1->Draw();
    p1->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetOptStat(0);
    //gPad->SetLogx();

    h_TriggerEff->Draw("");

    TAxis * xaxis1 = h_TriggerEff->GetXaxis();
    TAxis * yaxis1 = h_TriggerEff->GetYaxis();
    yaxis1->SetRangeUser(0.0, 1.1);
    xaxis1->SetTitle("Rigidity (GV)");
    yaxis1->SetTitle("Trigger Efficiency");    
    xaxis1->SetTitleFont(62);
    yaxis1->SetTitleFont(62);
    xaxis1->SetTitleSize(0.045);
    yaxis1->SetTitleSize(0.045);
    xaxis1->SetLabelFont(62);
    xaxis1->SetLabelSize(0.05);
    yaxis1->SetLabelFont(62);
    yaxis1->SetLabelSize(0.05);

    c1.SaveAs( (string("TriggerEff_") + MCName + PatternName + string(".pdf")).c_str() );
    */       


    TCanvas c2("c2","c2",1000,500);
    
    TPad *p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 1.0, 0);  //(name, title, xlow, ylow, xup, yup, color, bordersize, bordermode)
    p2->Draw();
    p2->cd();
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetOptStat(0);
    gPad->SetLogx();
    
    h_TriggerEff_noprescaling->Draw("HIST");

    TAxis * xaxis2 = h_TriggerEff_noprescaling->GetXaxis();
    TAxis * yaxis2 = h_TriggerEff_noprescaling->GetYaxis();

    xaxis2->SetLimits   (1.0, 800);
    xaxis2->SetRangeUser(1.0, 800);
    yaxis2->SetRangeUser(0.5, 0.85);

    xaxis2->SetTitle("Rigidity (GV)");
    yaxis2->SetTitle("Trigger Efficiency");

    xaxis2->SetTitleFont(62);
    yaxis2->SetTitleFont(62);
    xaxis2->SetTitleSize(0.06);
    yaxis2->SetTitleSize(0.06);
    xaxis2->SetLabelFont(62);
    xaxis2->SetLabelSize(0.05);
    yaxis2->SetLabelFont(62);
    yaxis2->SetLabelSize(0.05);

    xaxis2->SetTitleOffset(1.3);
    //yaxis2->SetTitleOffset();

    gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.15); 
    
    c2.SaveAs( (string("TriggerEff_noprescaling_") + MCName + PatternName + string(".pdf")).c_str() );

              

}



