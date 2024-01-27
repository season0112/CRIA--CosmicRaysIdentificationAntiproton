
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


void  Plot_1D_TrdSegmentsXZNumber(vector<double> binning, int i, TH1F *AntiprotonMC_TrdSegmentsXZNumber, TH1F *IssNegativeData_TrdSegmentsXZNumber, TH1F *PionMC_TrdSegmentsXZNumber, TH1F *ElectronMC_TrdSegmentsXZNumber){
        TCanvas cplot1D_TrdSegmentsXZNumber("cplot1D_TrdSegmentsXZNumber","cplot1D_TrdSegmentsXZNumber",1000,500);

        AntiprotonMC_TrdSegmentsXZNumber->Draw("HIST");
        //IssNegativeData_TrdSegmentsXZNumber->Draw("HIST same");
        PionMC_TrdSegmentsXZNumber->Draw("HIST same");
        ElectronMC_TrdSegmentsXZNumber->Draw("HIST same");

        AntiprotonMC_TrdSegmentsXZNumber->SetMarkerColor(3);
        AntiprotonMC_TrdSegmentsXZNumber->SetLineColor(3);
        AntiprotonMC_TrdSegmentsXZNumber->SetLineWidth(2.0);
        IssNegativeData_TrdSegmentsXZNumber->SetMarkerColor(4);
        IssNegativeData_TrdSegmentsXZNumber->SetLineColor(4);
        IssNegativeData_TrdSegmentsXZNumber->SetLineWidth(2.0);
        PionMC_TrdSegmentsXZNumber->SetMarkerColor(1);
        PionMC_TrdSegmentsXZNumber->SetLineColor(1);
        PionMC_TrdSegmentsXZNumber->SetLineWidth(2.0);
        ElectronMC_TrdSegmentsXZNumber->SetMarkerColor(2);
        ElectronMC_TrdSegmentsXZNumber->SetLineColor(2);
        ElectronMC_TrdSegmentsXZNumber->SetLineWidth(2.0);

        TLegend * leg_TrdSegmentsXZNumber = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
        //leg_TrdSegmentsXZNumber->AddEntry(IssNegativeData_TrdSegmentsXZNumber,"IssNegativeData","lp");
        leg_TrdSegmentsXZNumber->AddEntry(PionMC_TrdSegmentsXZNumber,"PionMC","lp");  // New way to select Pion MC.
        leg_TrdSegmentsXZNumber->AddEntry(ElectronMC_TrdSegmentsXZNumber,"ElectronMC","lp");
        leg_TrdSegmentsXZNumber->AddEntry(AntiprotonMC_TrdSegmentsXZNumber,"AntiprotonMC","lp");
        gStyle->SetLegendTextSize(0.03);
        leg_TrdSegmentsXZNumber->Draw();

        gStyle->SetOptStat(0);

        gStyle->SetErrorX(0);
        cplot1D_TrdSegmentsXZNumber.SaveAs( (std::string("cplot1D_TrdSegmentsXZNumber") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );

        cplot1D_TrdSegmentsXZNumber.SetLogy();
        cplot1D_TrdSegmentsXZNumber.SaveAs( (std::string("cplot1D_TrdSegmentsXZNumber") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string("_LogY.pdf")).c_str() );

        cplot1D_TrdSegmentsXZNumber.Close();
}


void Plot_1D_TrdSegmentsYZNumber(vector<double> binning, int i, TH1F *AntiprotonMC_TrdSegmentsYZNumber, TH1F *IssNegativeData_TrdSegmentsYZNumber, TH1F *PionMC_TrdSegmentsYZNumber, TH1F *ElectronMC_TrdSegmentsYZNumber){
        TCanvas cplot1D_TrdSegmentsYZNumber("cplot1D_TrdSegmentsYZNumber","cplot1D_TrdSegmentsYZNumber",1000,500);

        AntiprotonMC_TrdSegmentsYZNumber->Draw("HIST");
        //IssNegativeData_TrdSegmentsYZNumber->Draw("HIST same");
        PionMC_TrdSegmentsYZNumber->Draw("HIST same");
        ElectronMC_TrdSegmentsYZNumber->Draw("HIST same");

        AntiprotonMC_TrdSegmentsYZNumber->SetMarkerColor(3);
        AntiprotonMC_TrdSegmentsYZNumber->SetLineColor(3);
        AntiprotonMC_TrdSegmentsYZNumber->SetLineWidth(2.0);
        IssNegativeData_TrdSegmentsYZNumber->SetMarkerColor(4);
        IssNegativeData_TrdSegmentsYZNumber->SetLineColor(4);
        IssNegativeData_TrdSegmentsYZNumber->SetLineWidth(2.0);
        PionMC_TrdSegmentsYZNumber->SetMarkerColor(1);
        PionMC_TrdSegmentsYZNumber->SetLineColor(1);
        PionMC_TrdSegmentsYZNumber->SetLineWidth(2.0);
        ElectronMC_TrdSegmentsYZNumber->SetMarkerColor(2);
        ElectronMC_TrdSegmentsYZNumber->SetLineColor(2);
        ElectronMC_TrdSegmentsYZNumber->SetLineWidth(2.0);

        TLegend * leg_TrdSegmentsYZNumber = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
        //leg_TrdSegmentsYZNumber->AddEntry(IssNegativeData_TrdSegmentsYZNumber,"IssNegativeData","lp");
        leg_TrdSegmentsYZNumber->AddEntry(PionMC_TrdSegmentsYZNumber,"PionMC","lp");  // New way to select Pion MC.
        leg_TrdSegmentsYZNumber->AddEntry(ElectronMC_TrdSegmentsYZNumber,"ElectronMC","lp");
        leg_TrdSegmentsYZNumber->AddEntry(AntiprotonMC_TrdSegmentsYZNumber,"AntiprotonMC","lp");
        gStyle->SetLegendTextSize(0.03);
        leg_TrdSegmentsYZNumber->Draw();

        gStyle->SetOptStat(0);

        gStyle->SetErrorX(0);
        cplot1D_TrdSegmentsYZNumber.SaveAs( (std::string("cplot1D_TrdSegmentsYZNumber") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );

        cplot1D_TrdSegmentsYZNumber.SetLogy();
        cplot1D_TrdSegmentsYZNumber.SaveAs( (std::string("cplot1D_TrdSegmentsYZNumber") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string("_LogY.pdf")).c_str() );

        cplot1D_TrdSegmentsYZNumber.Close();                
}


void Plot_1D_TrdNumberOfHits(vector<double> binning, int i, TH1F *AntiprotonMC_TrdNumberOfHits, TH1F *IssNegativeData_TrdNumberOfHits, TH1F *PionMC_TrdNumberOfHits, TH1F *ElectronMC_TrdNumberOfHits){
        TCanvas cplot1D_TrdNumberOfHits("cplot1D_TrdNumberOfHits","cplot1D_TrdNumberOfHits",1000,500);

        AntiprotonMC_TrdNumberOfHits->Draw("HIST");
        //IssNegativeData_TrdNumberOfHits->Draw("HIST same");
        PionMC_TrdNumberOfHits->Draw("HIST same");
        ElectronMC_TrdNumberOfHits->Draw("HIST same");

        AntiprotonMC_TrdNumberOfHits->SetMarkerColor(3);
        AntiprotonMC_TrdNumberOfHits->SetLineColor(3);
        AntiprotonMC_TrdNumberOfHits->SetLineWidth(2.0);
        IssNegativeData_TrdNumberOfHits->SetMarkerColor(4);
        IssNegativeData_TrdNumberOfHits->SetLineColor(4);
        IssNegativeData_TrdNumberOfHits->SetLineWidth(2.0);
        PionMC_TrdNumberOfHits->SetMarkerColor(1);
        PionMC_TrdNumberOfHits->SetLineColor(1);
        PionMC_TrdNumberOfHits->SetLineWidth(2.0);
        ElectronMC_TrdNumberOfHits->SetMarkerColor(2);
        ElectronMC_TrdNumberOfHits->SetLineColor(2);
        ElectronMC_TrdNumberOfHits->SetLineWidth(2.0);

        TLegend * leg_TrdNumberOfHits = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
        //leg_TrdNumberOfHits->AddEntry(IssNegativeData_TrdNumberOfHits,"IssNegativeData","lp");
        leg_TrdNumberOfHits->AddEntry(PionMC_TrdNumberOfHits,"PionMC","lp");  // New way to select Pion MC.
        leg_TrdNumberOfHits->AddEntry(ElectronMC_TrdNumberOfHits,"ElectronMC","lp");
        leg_TrdNumberOfHits->AddEntry(AntiprotonMC_TrdNumberOfHits,"AntiprotonMC","lp");
        gStyle->SetLegendTextSize(0.03);
        leg_TrdNumberOfHits->Draw();

        gStyle->SetOptStat(0);

        gStyle->SetErrorX(0);
        cplot1D_TrdNumberOfHits.SaveAs( (std::string("cplot1D_TrdNumberOfHits") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );

        cplot1D_TrdNumberOfHits.SetLogy();
        cplot1D_TrdNumberOfHits.SaveAs( (std::string("cplot1D_TrdNumberOfHits") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string("_LogY.pdf")).c_str() );

        cplot1D_TrdNumberOfHits.Close();                
}

void Plot_1D_TRDVTracksSize(vector<double> binning, int i, TH1F *AntiprotonMC_TRDVTracksSize, TH1F *IssNegativeData_TRDVTracksSize, TH1F *PionMC_TRDVTracksSize, TH1F *ElectronMC_TRDVTracksSize){
        TCanvas cplot1D_TRDVTracksSize("cplot1D_TRDVTracksSize","cplot1D_TRDVTracksSize",1000,500);

        AntiprotonMC_TRDVTracksSize->Draw("HIST");
        //IssNegativeData_TRDVTracksSize->Draw("HIST same");
        PionMC_TRDVTracksSize->Draw("HIST same");
        ElectronMC_TRDVTracksSize->Draw("HIST same");

        AntiprotonMC_TRDVTracksSize->SetMarkerColor(3);
        AntiprotonMC_TRDVTracksSize->SetLineColor(3);
        AntiprotonMC_TRDVTracksSize->SetLineWidth(2.0);
        IssNegativeData_TRDVTracksSize->SetMarkerColor(4);
        IssNegativeData_TRDVTracksSize->SetLineColor(4);
        IssNegativeData_TRDVTracksSize->SetLineWidth(2.0);
        PionMC_TRDVTracksSize->SetMarkerColor(1);
        PionMC_TRDVTracksSize->SetLineColor(1);
        PionMC_TRDVTracksSize->SetLineWidth(2.0);
        ElectronMC_TRDVTracksSize->SetMarkerColor(2);
        ElectronMC_TRDVTracksSize->SetLineColor(2);
        ElectronMC_TRDVTracksSize->SetLineWidth(2.0);

        TLegend * leg_TRDVTracksSize = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
        //leg_TRDVTracksSize->AddEntry(IssNegativeData_TRDVTracksSize,"IssNegativeData","lp");
        leg_TRDVTracksSize->AddEntry(PionMC_TRDVTracksSize,"PionMC","lp");  // New way to select Pion MC.
        leg_TRDVTracksSize->AddEntry(ElectronMC_TRDVTracksSize,"ElectronMC","lp");
        leg_TRDVTracksSize->AddEntry(AntiprotonMC_TRDVTracksSize,"AntiprotonMC","lp");
        gStyle->SetLegendTextSize(0.03);
        leg_TRDVTracksSize->Draw();

        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        cplot1D_TRDVTracksSize.SaveAs( (std::string("cplot1D_TRDVTracksSize") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );

        cplot1D_TRDVTracksSize.SetLogy();
        cplot1D_TRDVTracksSize.SaveAs( (std::string("cplot1D_TRDVTracksSize") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string("_LogY.pdf")).c_str() );

        cplot1D_TRDVTracksSize.Close();
}


void Plot_1D_ACCHits(vector<double> binning, int i, TH1F *AntiprotonMC_ACCHits, TH1F *IssNegativeData_ACCHits, TH1F *PionMC_ACCHits, TH1F *ElectronMC_ACCHits){
        TCanvas cplot1D_ACCHits("cplot1D_ACCHits","cplot1D_ACCHits",1000,500);

        AntiprotonMC_ACCHits->Draw("HIST");
        //IssNegativeData_ACCHits->Draw("HIST same");
        PionMC_ACCHits->Draw("HIST same");
        ElectronMC_ACCHits->Draw("HIST same");

        AntiprotonMC_ACCHits->SetMarkerColor(3);
        AntiprotonMC_ACCHits->SetLineColor(3);
        AntiprotonMC_ACCHits->SetLineWidth(2.0);
        IssNegativeData_ACCHits->SetMarkerColor(4);
        IssNegativeData_ACCHits->SetLineColor(4);
        IssNegativeData_ACCHits->SetLineWidth(2.0);
        PionMC_ACCHits->SetMarkerColor(1);
        PionMC_ACCHits->SetLineColor(1);
        PionMC_ACCHits->SetLineWidth(2.0);
        ElectronMC_ACCHits->SetMarkerColor(2);
        ElectronMC_ACCHits->SetLineColor(2);
        ElectronMC_ACCHits->SetLineWidth(2.0);

        TLegend * leg_ACCHits = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
        //leg_ACCHits->AddEntry(IssNegativeData_ACCHits,"IssNegativeData","lp");
        leg_ACCHits->AddEntry(PionMC_ACCHits,"PionMC","lp");  // New way to select Pion MC.
        leg_ACCHits->AddEntry(ElectronMC_ACCHits,"ElectronMC","lp");
        leg_ACCHits->AddEntry(AntiprotonMC_ACCHits,"AntiprotonMC","lp");
        gStyle->SetLegendTextSize(0.03);
        leg_ACCHits->Draw();

        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        cplot1D_ACCHits.SaveAs( (std::string("cplot1D_ACCHits") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );

        cplot1D_ACCHits.SetLogy();
        cplot1D_ACCHits.SaveAs( (std::string("cplot1D_ACCHits") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string("_LogY.pdf")).c_str() );

        cplot1D_ACCHits.Close();       
}

void Plot_1D_EcalBDT_EnergyD(vector<double> binning, int i, TH1F *AntiprotonMC_EcalBDT_EnergyD, TH1F *IssNegativeData_EcalBDT_EnergyD, TH1F *PionMC_EcalBDT_EnergyD, TH1F *ElectronMC_EcalBDT_EnergyD){
        TCanvas cplot1D_EcalBDT_EnergyD("cplot1D_EcalBDT_EnergyD","cplot1D_EcalBDT_EnergyD",1000,500);

        AntiprotonMC_EcalBDT_EnergyD->Draw("HIST");
        //IssNegativeData_EcalBDT_EnergyD->Draw("HIST same");
        PionMC_EcalBDT_EnergyD->Draw("HIST same");
        ElectronMC_EcalBDT_EnergyD->Draw("HIST same");

        AntiprotonMC_EcalBDT_EnergyD->SetMarkerColor(3);
        AntiprotonMC_EcalBDT_EnergyD->SetLineColor(3);
        AntiprotonMC_EcalBDT_EnergyD->SetLineWidth(2.0);
        IssNegativeData_EcalBDT_EnergyD->SetMarkerColor(4);
        IssNegativeData_EcalBDT_EnergyD->SetLineColor(4);
        IssNegativeData_EcalBDT_EnergyD->SetLineWidth(2.0);
        PionMC_EcalBDT_EnergyD->SetMarkerColor(1);
        PionMC_EcalBDT_EnergyD->SetLineColor(1);
        PionMC_EcalBDT_EnergyD->SetLineWidth(2.0);
        ElectronMC_EcalBDT_EnergyD->SetMarkerColor(2);
        ElectronMC_EcalBDT_EnergyD->SetLineColor(2);
        ElectronMC_EcalBDT_EnergyD->SetLineWidth(2.0);

        TLegend * leg_EcalBDT_EnergyD = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
        //leg_EcalBDT_EnergyD->AddEntry(IssNegativeData_EcalBDT_EnergyD,"IssNegativeData","lp");
        leg_EcalBDT_EnergyD->AddEntry(PionMC_EcalBDT_EnergyD,"PionMC","lp");  // New way to select Pion MC.
        leg_EcalBDT_EnergyD->AddEntry(ElectronMC_EcalBDT_EnergyD,"ElectronMC","lp");
        leg_EcalBDT_EnergyD->AddEntry(AntiprotonMC_EcalBDT_EnergyD,"AntiprotonMC","lp");
        gStyle->SetLegendTextSize(0.03);
        leg_EcalBDT_EnergyD->Draw();

        gStyle->SetOptStat(0);
        gStyle->SetErrorX(0);

        cplot1D_EcalBDT_EnergyD.SaveAs( (std::string("cplot1D_EcalBDT_EnergyD") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );

        cplot1D_EcalBDT_EnergyD.SetLogy();
        cplot1D_EcalBDT_EnergyD.SaveAs( (std::string("cplot1D_EcalBDT_EnergyD") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string("_LogY.pdf")).c_str() );

        cplot1D_EcalBDT_EnergyD.Close();                
}


void Plot_1D_TrdLogLikelihoodRatioProtonHeliumTracker(vector<double> binning, int i, TH1F *AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker, TH1F *IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker, TH1F *PionMC_TrdLogLikelihoodRatioProtonHeliumTracker, TH1F *ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker){
    TCanvas cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker("cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker","cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker",1000,500);

    AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker->Draw("HIST");
    //IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker->Draw("HIST same");
    PionMC_TrdLogLikelihoodRatioProtonHeliumTracker->Draw("HIST same");
    ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker->Draw("HIST same");

    AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetMarkerColor(3);
    AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineColor(3);
    AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineWidth(2.0);
    IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker->SetMarkerColor(4);
    IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineColor(4);
    IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineWidth(2.0);
    PionMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetMarkerColor(1);
    PionMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineColor(1);
    PionMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineWidth(2.0);
    ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetMarkerColor(2);
    ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineColor(2);
    ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker->SetLineWidth(2.0);

    TLegend * leg_TrdLogLikelihoodRatioProtonHeliumTracker = new TLegend(0.7,0.7,0.9,0.85); //(xmin, ymin, xmax, ymax)
    //leg_TrdLogLikelihoodRatioProtonHeliumTracker->AddEntry(IssNegativeData_TrdLogLikelihoodRatioProtonHeliumTracker,"IssNegativeData","lp");
    leg_TrdLogLikelihoodRatioProtonHeliumTracker->AddEntry(PionMC_TrdLogLikelihoodRatioProtonHeliumTracker,"PionMC","lp");  // New way to select Pion MC.
    leg_TrdLogLikelihoodRatioProtonHeliumTracker->AddEntry(ElectronMC_TrdLogLikelihoodRatioProtonHeliumTracker,"ElectronMC","lp");
    leg_TrdLogLikelihoodRatioProtonHeliumTracker->AddEntry(AntiprotonMC_TrdLogLikelihoodRatioProtonHeliumTracker,"AntiprotonMC","lp");
    gStyle->SetLegendTextSize(0.03);
    leg_TrdLogLikelihoodRatioProtonHeliumTracker->Draw();

    gStyle->SetOptStat(0);
    gStyle->SetErrorX(0);

    cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker.SaveAs( (std::string("cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string(".pdf")).c_str() );

    cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker.SetLogy();
    cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker.SaveAs( (std::string("cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker") + std::string("_") + to_string_with_precision(binning.at(i), 2) + std::string("_") + to_string_with_precision(binning.at(i+1), 2) + std::string("_LogY.pdf")).c_str() );

    cplot1D_TrdLogLikelihoodRatioProtonHeliumTracker.Close();
}



