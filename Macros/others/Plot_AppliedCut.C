// root -L /home/bo791269/Software/AntiprotonAnalysis/Libraries/AntiprotonBinning/AntiprotonBinning.C Plot_AppliedCut.C

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

void Plot_AppliedCut(){

    std::list<std::string> variableAll;
    //variableAll.push_back("InnerTrackerCharge"); // Only high range
    variableAll.push_back("LowerTofCharge");
    //variableAll.push_back("UpperTofCharge");
    //variableAll.push_back("TrdLogLikelihoodRatioProtonHeliumTracker");

    double lowedge   = 0;
    double highedge  = 5;
    double binnumber = 500;

    std::list<std::string> MCNameAll;

    /*
    //// High range
    std::string path = std::string(getenv("HPCHIGHENERGYDATADIR")) + std::string("/");
    std::string treename = "ExampleAnalysisTree";

    MCNameAll.push_back("B1042_antipr.pl1.1800_7.6_all_Tree_");
    //MCNameAll.push_back("B1042_pr.pl1.1800_7.6_all_Tree_positive_");
    //MCNameAll.push_back("B1042_pr.pl1.flux.l1a9.2016000_7.6_all_Tree_positive_"); //big
    //MCNameAll.push_back("B1220_pr.pl1ph.021000_7.8_all_Tree_positive_");
    //MCNameAll.push_back("B1220_pr.pl1phpsa.l19.5016000.4_00_7.8_all_Tree_positive_"); //big

    std::string issReduced = std::string("");

    int index_start = 50;
    int index_end   = 54;

    int digit = 0;

    */

    //// Low range
    std::string path = std::string(getenv("HPCLOWENERGYDIR")) + std::string("/totalall/");
    std::string treename = "AntiprotonLowEnergyTree";

    //MCNameAll.push_back("B1042_antipr.pl1.1800_7.6_all_Tree_");
    //MCNameAll.push_back("B1042_pr.pl1.1800_7.6_all_Tree_positive_");
    //MCNameAll.push_back("B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_");
    MCNameAll.push_back("B1042_antipr.pl1.1800_7.6_all_Tree_negative_RemoveTOFChargeCut");

    std::string issReduced = std::string("_test");
     
    int index_start = 17; //10:3.29 14:4.88
    int index_end   = 18;

    int digit = 2;

    //TCut PionTemplateDataCut = (std::string("RichIsNaF==0 && abs(RichBeta-(-Rigidity)/sqrt(0.139^2+Rigidity^2))<0.002 && TrdLogLikelihoodRatioElectronProtonTracker>") + std::to_string(TrdLOW) + std::string("&& TrdSegmentsXZNumber>1 && TrdSegmentsYZNumber>1 && TrdNumberOfHits>50")).c_str();
    //TCut ProtonSelection = ( std::string("TrdLogLikelihoodRatioElectronProtonTracker>0.6")  ).c_str();
    //TCut ProtonSelection = ( std::string("TrdLogLikelihoodRatioProtonHeliumTracker>-0.05 && TrdLogLikelihoodRatioProtonHeliumTracker<0.3")  ).c_str();
    TCut ProtonSelection = ( std::string("TrdLogLikelihoodRatioElectronProtonTracker>0.6 && TrdLogLikelihoodRatioProtonHeliumTracker>-0.05 && TrdLogLikelihoodRatioProtonHeliumTracker<0.3")  ).c_str();
    //TCut ProtonSelection;

    /*
    Pattern
    UpperTofBeta
    LowerTofBeta
    TofMassonecharge
    EcalBDT_EnergyD
    TrdLogLikelihoodRatioElectronProtonTracker
    TrdLogLikelihoodRatioProtonHeliumTracker
    TrdKCharge
    */

    ///////////////////////////////////////////////////////

    std::vector<double> PublishedPRLBinEdge ( AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.begin()+1, AntiprotonNewBinning::AntiprotonAllBinning::AntiprotonBinning_450.end()); //size=58;  index=50

    //// Loop
    for (std::string variable:variableAll){
        cout<< "variable is " << variable <<endl;

        for (std::string MCName:MCNameAll){
            cout<< "MCName is " << MCName <<endl;

            for (int i=index_start; i<index_end; i++ ){
                cout << '\n' << endl;
                cout << i << endl;
                cout<< PublishedPRLBinEdge.at(i);

                TChain *f_MCPbar = new TChain((treename).c_str());
                //f_MCPbar->AddFile(  (path +  MCName + to_string_with_precision(PublishedPRLBinEdge.at(i),digit) + std::string("_") + to_string_with_precision(PublishedPRLBinEdge.at(i+1), digit) + std::string(".root")).c_str() );
                f_MCPbar->AddFile(  (path +  MCName + std::string(".root")).c_str() );
                f_MCPbar->Draw( (variable + std::string(">>h_MCPbar(") + to_string_with_precision(binnumber, 2) + std::string(",") + to_string_with_precision(lowedge, 2) + std::string(",") + to_string_with_precision(highedge,2) + std::string(")") ).c_str(), ProtonSelection);
                TH1F *h_MCPbar = (TH1F*)gDirectory->Get("h_MCPbar");


                TChain *f_ISSdata = new TChain((treename).c_str());
                //f_ISSdata->AddFile(  (path + std::string("B1130_pass7_7.8_all_Tree_positive_May2015_") + to_string_with_precision(PublishedPRLBinEdge.at(i), digit) + std::string("_") + to_string_with_precision(PublishedPRLBinEdge.at(i+1), digit) + issReduced + std::string(".root")).c_str() );
                f_ISSdata->AddFile(  (path + std::string("B1130_pass7_7.8_all_Tree_positive_RemoveTOFChargeCut") + std::string(".root")).c_str() );
                f_ISSdata->Draw( (variable + std::string(">>h_ISSdata(") + to_string_with_precision(binnumber,2) + std::string(",") + to_string_with_precision(lowedge, 2) + std::string(",") + to_string_with_precision(highedge,2) + std::string(")") ).c_str(), ProtonSelection );
                TH1F *h_ISSdata = (TH1F*)gDirectory->Get("h_ISSdata");

                double scale_MCPbar = 1.0/h_MCPbar->Integral();
                h_MCPbar->Scale(scale_MCPbar);
                double scale_ISSdata = 1.0/h_ISSdata->Integral();
                h_ISSdata->Scale(scale_ISSdata);

                h_MCPbar->Sumw2();
                h_ISSdata->Sumw2();

                TGraphErrors *g_ISSdata = new TGraphErrors(h_ISSdata);


                //// Plot
                TCanvas c1("c1","c1",1000, 800);

                h_MCPbar->Draw("HIST");
                //h_ISSdata->Draw("same HIST");
                g_ISSdata->Draw("E1 P");

                h_MCPbar->SetTitle("");
                gStyle->SetOptStat(0);

                TAxis * xaxis = h_MCPbar->GetXaxis();
                TAxis * yaxis = h_MCPbar->GetYaxis();
                xaxis->SetTitle("Q_{LowTOF}");
                yaxis->SetTitle("Normalized events");
                yaxis->SetRangeUser(0.00001, 0.1);  // (0.000000001, 0.1)
                xaxis->SetTitleFont(62);
                yaxis->SetTitleFont(62);
                xaxis->SetTitleSize(0.045);
                yaxis->SetTitleSize(0.045);
                xaxis->SetLabelFont(62);
                yaxis->SetLabelFont(62);
                xaxis->SetLabelSize(0.05);
                yaxis->SetLabelSize(0.05);
                //xaxis->SetTitleOffset(0.8);
                //yaxis->SetTitleOffset(0.8);

                h_MCPbar->SetLineColor(4);
                h_MCPbar->SetLineWidth(2);
                h_MCPbar->SetMarkerSize(0.9);
                h_MCPbar->SetMarkerColor(4);
                h_MCPbar->SetMarkerStyle(15);
                h_ISSdata->SetLineColor(2);
                h_ISSdata->SetLineWidth(3);
                h_ISSdata->SetMarkerSize(0.9);
                g_ISSdata->SetMarkerColor(2);
                g_ISSdata->SetMarkerStyle(15);
                g_ISSdata->SetMarkerSize(0.9);

                TLine *line_Cut = new TLine(2, 0, 2, 0.1);
                line_Cut->SetLineColor(kBlack);
                line_Cut->SetLineWidth(3);
                //line_Cut->SetNDC(1);
                line_Cut->SetLineStyle(9);
                line_Cut->Draw();

                TLegend *legend = new TLegend(0.70, 0.7, 0.89, 0.85);
                legend->AddEntry(h_MCPbar , "MC" , "l");
                legend->AddEntry(g_ISSdata, "ISS", "p");
                legend->SetTextSize(0.04);
                legend->SetTextFont(62);
                legend->SetBorderSize(0);
                legend->Draw();

                xaxis->SetTitleOffset(1.3);
                gPad->SetBottomMargin(0.2);
                gPad->SetLeftMargin(0.16);

                gPad->SetLogy();


                c1.Update();
                //c1.SaveAs( (string("Compare_") + variable + string("_") + MCName + string("_") + to_string_with_precision(PublishedPRLBinEdge.at(i), digit) + std::string("_") + to_string_with_precision(PublishedPRLBinEdge.at(i+1), digit) + string(".pdf")).c_str());
                c1.SaveAs( (string("Compare_") + variable + string("_") + MCName + string(".pdf")).c_str());

            }

        }
    }
}


