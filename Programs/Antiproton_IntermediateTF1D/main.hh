//#include "CutDefinition.hh"
#include "BinningTools.hh"

//// For Time Averaged AND Time Dependent.
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


std::tuple<std::vector<double>> LoadRichbetaCut(std::string lowpath, std::string signalefficiency, std::string binmerge){
    std::vector<double> v_richbetacut;
    std::ifstream richbeta;
    richbeta.open( lowpath + std::string("/RichBetaCutValue_eff_") + signalefficiency + std::string("_"+ binmerge +".txt") );
    std::string cutvalue;
    assert(v_richbetacut.empty());
    while (getline(richbeta, cutvalue)) {
        v_richbetacut.push_back(stod(cutvalue));
        //std::cout << stod(cutvalue) <<std::endl;
        //std::cout<< v_richbetacut.size()  << std::endl;
    }
    return {v_richbetacut};
}


std::tuple<std::vector<double>> LoadTRDLogLikelihoodCut(std::string TRDefficiency, std::string binmerge, std::string lowpath){
    std::vector<double> v_TRDcut;
    std::ifstream TRDfile;
    TRDfile.open( lowpath + std::string("TRDLogLikelihood_CutValue_eff_") + TRDefficiency + std::string("_") + binmerge + std::string(".txt") );
    std::string trdcutvalue;
    assert(v_TRDcut.empty());
    while (getline(TRDfile, trdcutvalue)) {
        v_TRDcut.push_back(stod(trdcutvalue));
    }
    TRDfile.close();
    return {v_TRDcut};
}


//// For Time Averaged ONLY
std::tuple<TH1D, TH1D> LoadPublishedResultError(std::vector<double> subbincenter, std::vector<double> v_StatisticError_relative, std::vector<double> subbinedge, std::vector<double> v_SystematicError_relative){
    TGraph *g_Published_Statistic_Error_Relative = new TGraph(subbincenter.size(), subbincenter.data(), v_StatisticError_relative.data());
    TH1D h_Published_Statistic_Error_Relative = TH1D("", "", subbincenter.size(), subbinedge.data());
    Utilities::ConvertToHistogram(g_Published_Statistic_Error_Relative, h_Published_Statistic_Error_Relative);

    TGraph *g_Published_Systematic_Error_Relative = new TGraph(subbincenter.size(), subbincenter.data(), v_SystematicError_relative.data());
    TH1D h_Published_Systematic_Error_Relative = TH1D("", "", subbincenter.size(), subbinedge.data());
    Utilities::ConvertToHistogram(g_Published_Systematic_Error_Relative, h_Published_Systematic_Error_Relative);

    return {h_Published_Statistic_Error_Relative, h_Published_Systematic_Error_Relative};
}


std::tuple<TGraphErrors *> LoadEffectiveAcceptance(std::string lowpath){
    TFile *f_eff_antiproton = new TFile( (lowpath + std::string("/EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root")).c_str() );
    TFile *f_eff_proton     = new TFile( (lowpath + std::string("/EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root")).c_str() );
    TGraphAsymmErrors *Acceptance_antiproton = (TGraphAsymmErrors*)f_eff_antiproton->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_proton     = (TGraphAsymmErrors*)f_eff_proton    ->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    int number = Acceptance_antiproton->GetN();
    Double_t xa[number],Aa[number],ea[number],xp[number],Ap[number],ep[number];

    TGraphErrors *Effective_Acceptance = new TGraphErrors(number);
    for (int p = 0; p < number; p=p+1) {       //from 0, 0 is 0.8-1.0, last is from 822 to 1130. so for B1042MC there two are 0 (generated momentum from 1 to 800).
        Acceptance_antiproton->GetPoint(p,xa[p],Aa[p]);
        Acceptance_proton    ->GetPoint(p,xp[p],Ap[p]);
        ea[p] = Acceptance_antiproton->GetErrorY(p);
        ep[p] = Acceptance_proton    ->GetErrorY(p);
        Effective_Acceptance->SetPoint     (p, xa[p], Ap[p]/Aa[p]);
        Effective_Acceptance->SetPointError(p, 0.0  , sqrt(pow(ep[p],2)/pow(Aa[p],2)+pow(Ap[p],2)/pow(Aa[p],4)*pow(ea[p],2)));
    }
    f_eff_antiproton->Close();
    f_eff_proton->Close();

    return {Effective_Acceptance};
}


std::tuple<TH1F *> LoadAntiprotonMCTemplate(std::string lowpath, std::string binmerge, std::vector<double> binning, int index, TCut AllCut){
    TChain *fantiproton = new TChain("AntiprotonLowEnergyTree");
    for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
        fantiproton->AddFile(( lowpath + std::string("/B1042_antipr.pl1.1800_7.6_all_Tree_negative_") + doubleToString(binning[index+binindex]) + std::string("_") + doubleToString(binning[index+binindex+1]) + std::string(".root")).c_str());}
    fantiproton->Draw("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_antiproton(50,-2.0, 0)", AllCut);
    TH1F *template_antiproton = (TH1F*)gDirectory->Get("th1f_antiproton");
    template_antiproton->SetTitle("");
    delete fantiproton;
    return {template_antiproton};
}


std::tuple<TH1F *> LoadElectronMCTemplate(std::string lowpath, std::string binmerge, std::vector<double> binning, int index, TCut AllCut){
    TChain *felectron = new TChain("AntiprotonLowEnergyTree");
    for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
        felectron->AddFile(( lowpath + std::string("/B1091_el.pl1.0_25200_7.6_all_Tree_negative_") + doubleToString(binning[index+binindex]) + std::string("_") + doubleToString(binning[index+binindex+1]) + std::string(".root")).c_str());}
    felectron->Draw("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_electron(50,-2.0, 0)", AllCut);
    TH1F *template_electron = (TH1F*)gDirectory->Get("th1f_electron");
    template_electron->SetTitle("Electron");
    delete felectron;
    return {template_electron};
}


std::tuple<TH1F *> LoadPionMCTemplate(std::string lowpath, std::string binmerge, std::vector<double> binning, int index, TCut AllCut){
    TChain *fpion = new TChain("AntiprotonLowEnergyTree");
    for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
        fpion->AddFile(( lowpath + std::string("/B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + doubleToString(binning[index+binindex]) + std::string("_") + doubleToString(binning[index+binindex+1]) + std::string(".root")).c_str());}
    fpion->Draw("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_pion(50, -2.0, 0)", AllCut);
    TH1F *template_pion =  (TH1F*)gDirectory->Get("th1f_pion");
    template_pion->SetTitle("Pion");
    delete fpion;
    return {template_pion};
}


std::tuple<TH1F *, TChain *>LoadISSPositiveData(std::string lowpath, std::string binmerge, std::vector<double> binning, int index, TCut PositiveCut, std::string issversion, std::string TestMode, double TrdLOW, double TrdHIGH, double TrdBINNUMBER){

    TChain *fpass7_positive = new TChain("AntiprotonLowEnergyTree");
    if (issversion == "pass7.8"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            std::cout<< "left:" << left << std::endl;
            std::cout<< "right:" << right << std::endl;
            if (TestMode == "Yes"){
                fpass7_positive->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_positive_") + left + std::string("_") + right + std::string("_test.root")).c_str());}
            else{
                fpass7_positive->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_positive_") + left + std::string("_") + right + std::string(".root")).c_str());}}
    }
    else if (issversion == "PhysicsReport"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13)){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            std::cout<< "left:" << left << std::endl;
            std::cout<< "right:" << right << std::endl;
            if (TestMode == "Yes"){
                 fpass7_positive->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_positive_Nov2017_")+ left + std::string("_") + right + std::string("_test.root")).c_str());}
            else{
                fpass7_positive->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_positive_Nov2017_")+ left + std::string("_") + right + std::string(".root")).c_str());}}
    }
    else if (issversion == "2016paper"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13)){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            std::cout<< "left:" << left << std::endl;
            std::cout<< "right:" << right << std::endl;
            if (TestMode == "Yes"){
                fpass7_positive->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_positive_May2015_")+ left + std::string("_") + right + std::string("_test.root")).c_str());}
            else{
                fpass7_positive->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_positive_May2015_")+ left + std::string("_") + right + std::string(".root")).c_str());}}
    }

    fpass7_positive->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_positive(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")") ).c_str(), PositiveCut);
    TH1F *data_pass7_positive =  (TH1F*)gDirectory->Get("th1f_positive");
    data_pass7_positive->SetTitle("Antiproton");

    return {data_pass7_positive, fpass7_positive};
}



std::tuple<TH1F *, TH1F *, TH1F *, TChain *>LoadISSNegativeData(std::string lowpath, std::string binmerge, std::vector<double> binning, int index, TCut NegativeCut, TCut ElectronTemplateDataCut, TCut PionTemplateDataCut, std::string issversion, std::string TestMode, double TrdLOW, double TrdHIGH, double TrdBINNUMBER){

    TChain *fpass7_negative = new TChain("AntiprotonLowEnergyTree");
    if (issversion == "pass7.8"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13)){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            std::cout<< "left:" << left << std::endl;
            std::cout<< "right:" << right << std::endl;
            fpass7_negative->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_negative_") + left + std::string("_") + right + std::string(".root")).c_str());}
    }
    else if (issversion == "PhysicsReport"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            std::cout<< "left:" << left << std::endl;
            std::cout<< "right:" << right << std::endl;
            fpass7_negative->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_negative_Nov2017_")+ left + std::string("_") + right + std::string(".root")).c_str());}
    }
    else if (issversion == "2016paper"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            std::cout<< "left:" << left << std::endl;
            std::cout<< "right:" << right << std::endl;
            fpass7_negative->AddFile(( lowpath + std::string("/B1130_pass7_7.8_all_Tree_negative_May2015_")+ left + std::string("_") + right + std::string(".root")).c_str());}
    }


    fpass7_negative->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_negative(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")") ).c_str(), NegativeCut);
    TH1F *data_pass7_negative =  (TH1F*)gDirectory->Get("th1f_negative");
    data_pass7_negative->SetTitle("");

    fpass7_negative->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_ElectronDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")") ).c_str(), ElectronTemplateDataCut);
    TH1F *template_electron_Data =  (TH1F*)gDirectory->Get("th1f_ElectronDataTemplate");
    template_electron_Data->SetTitle("Electron");

    fpass7_negative->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_PionDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")") ).c_str(), PionTemplateDataCut);
    TH1F *template_pion_Data =  (TH1F*)gDirectory->Get("th1f_PionDataTemplate");
    template_pion_Data->SetTitle("Pion");

    return {data_pass7_negative, template_electron_Data, template_pion_Data, fpass7_negative};
}


std::tuple<TH1F *, TH1F *, TH1F *, TH1F *, double>LoadTemplates(std::string lowpath, std::string binmerge, std::vector<double> binning, int index, std::string issversion, std::string ParametrilizedMode, double RescaleFactor, int numberindex, std::string trdeff, std::string FullRange){

    double ProtonNumber;
    std::string left  = doubleToString(binning[index]);
    std::string right = doubleToString(binning[index + stoi(binmerge)]);
    if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
        left = to_string_with_precision(binning[index], 1);}
    if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
        right = to_string_with_precision(binning[index + stoi(binmerge)], 1);}

    std::string nameindex; 
    std::string templateindex;
    if (ParametrilizedMode == "Yes"){
        nameindex     = "_new";
        //templateindex = "_fit";
        templateindex = std::string("_RandomIndex_") + std::to_string(numberindex);
    }
    else if (ParametrilizedMode == "No"){
        nameindex     = "";
        templateindex = "";
    }

    std::string TemplateRootName;
    if (FullRange == "Yes"){
        TemplateRootName = "averaged_ratio_fullRange_";
    }
    else{
        TemplateRootName = "averaged_ratio_";
    }

    TFile *f_histogram = new TFile( (std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/RigidityRootFiles/") + TemplateRootName + left + std::string("_") + right + std::string("_") + issversion + nameindex + std::string(".root")).c_str(), "OPEN" );
    TH1F *data_pass7_positive    = (TH1F*) f_histogram->Get( (std::string("data_pass7_positive")   + std::string("_TRDeff_")    + trdeff + templateindex).c_str() );
    TH1F *template_electron_Data = (TH1F*) f_histogram->Get( (std::string("template_electron_Data")+ std::string("_TRDeff_")    + trdeff + templateindex).c_str() );
    TH1F *template_pion_Data     = (TH1F*) f_histogram->Get( (std::string("template_pion_Data")    + std::string("_TRDeff_")    + trdeff + templateindex).c_str() );
    TH1F *data_pass7_negative    = (TH1F*) f_histogram->Get( (std::string("data_pass7_negative")   + std::string("_TRDeff_")    + trdeff                ).c_str() );
    //f_histogram->Close();
    //TH1F *data_pass7_positive    = (TH1F*) f_histogram->Get( (std::string("data_pass7_positive")    + templateindex).c_str() );
    //TH1F *template_electron_Data = (TH1F*) f_histogram->Get( (std::string("template_electron_Data") + templateindex).c_str() );
    //TH1F *template_pion_Data     = (TH1F*) f_histogram->Get( (std::string("template_pion_Data")     + templateindex).c_str() );
    //TH1F *data_pass7_negative    = (TH1F*) f_histogram->Get( (std::string("data_pass7_negative")                   ).c_str() );

    if (ParametrilizedMode == "Yes"){
        TH1F *data_pass7_positive_original = (TH1F*) f_histogram->Get("data_pass7_positive");
        ProtonNumber = data_pass7_positive_original->GetEntries();
    }
    else if (ParametrilizedMode == "No"){
        ProtonNumber = data_pass7_positive->GetEntries();
    }

    ProtonNumber = ProtonNumber * RescaleFactor;

    return {data_pass7_positive, template_electron_Data, template_pion_Data, data_pass7_negative, ProtonNumber};

}


void PrintFitResult(std::vector<double> result, std::vector<double> ResultError, double CHI2dof, double ProtonNumber){
    std::cout << "Tempalte Fit Number Result: "      << result                      << std::endl;
    std::cout << "Tempalte Fit Number ResultError: " << ResultError                 << std::endl;
    std::cout << "chi2dof is: "                      << CHI2dof                     << std::endl;
    std::cout << "Proton Numbers:"                   << ProtonNumber                << std::endl;
    std::cout << "PbarOverP Ratio:"                  << result[0]/ProtonNumber      << std::endl;
    std::cout << "PbarOverP RatioError:"             << ResultError[0]/ProtonNumber << std::endl;
}


void Plot_Pattern_Percentage(std::string richcut, std::string issversion, std::string binmerge, std::vector<double> v_pattern0_percentage, std::vector<double> v_pattern1_percentage, std::vector<double> v_pattern2_percentage, std::vector<double> v_pattern4_percentage, std::vector<double> subbincenter, TChain *fpass7_negative, TCut RichBetaCut, TCut pattern0_percentage, TCut pattern1_percentage, TCut pattern2_percentage, TCut pattern4_percentage, TCut AllCut){
    // Calculate 
    TTree *tpass7_negative_all = (TTree*)fpass7_negative->GetTree();
    gROOT->cd();

    TTree* tpass7_negative_afterRICHcut = tpass7_negative_all->CopyTree(RichBetaCut);
    gROOT->cd();

    TTree* tpass7_negative_pattern0 = tpass7_negative_all->CopyTree(pattern0_percentage);
    v_pattern0_percentage.push_back((double)(tpass7_negative_pattern0->GetEntries())/(double)(tpass7_negative_afterRICHcut->GetEntries()));
    gROOT->cd();

    TTree* tpass7_negative_pattern1 = tpass7_negative_all->CopyTree(pattern1_percentage);
    v_pattern1_percentage.push_back((double)(tpass7_negative_pattern1->GetEntries())/(double)(tpass7_negative_afterRICHcut->GetEntries()));
    gROOT->cd();

    TTree* tpass7_negative_pattern2 = tpass7_negative_all->CopyTree(pattern2_percentage);
    v_pattern2_percentage.push_back((double)(tpass7_negative_pattern2->GetEntries())/(double)(tpass7_negative_afterRICHcut->GetEntries()));
    gROOT->cd();

    TTree* tpass7_negative_pattern4 = tpass7_negative_all->CopyTree(pattern4_percentage);
    v_pattern4_percentage.push_back((double)(tpass7_negative_pattern4->GetEntries())/(double)(tpass7_negative_afterRICHcut->GetEntries()));
    gROOT->cd();

    TTree* tpass7_negative = tpass7_negative_all->CopyTree(AllCut);
    
    TGraph *g_fullspan_percentage = new TGraph(v_pattern0_percentage.size(), subbincenter.data(), v_pattern0_percentage.data());
    TGraph *g_pattern1_percentage = new TGraph(v_pattern1_percentage.size(), subbincenter.data(), v_pattern1_percentage.data());
    TGraph *g_pattern2_percentage = new TGraph(v_pattern2_percentage.size(), subbincenter.data(), v_pattern2_percentage.data());
    TGraph *g_pattern4_percentage = new TGraph(v_pattern4_percentage.size(), subbincenter.data(), v_pattern4_percentage.data());    

    // Plot
    TCanvas allpattern_percentage("allpattern_percentage","allpattern_percentage",1000,500);

    g_fullspan_percentage->SetMarkerSize(2);
    g_fullspan_percentage->SetMarkerStyle(8);
    g_fullspan_percentage->SetMarkerColor(2);
    g_pattern1_percentage->SetMarkerSize(2);
    g_pattern1_percentage->SetMarkerStyle(8);
    g_pattern1_percentage->SetMarkerColor(4);
    g_pattern2_percentage->SetMarkerSize(2);
    g_pattern2_percentage->SetMarkerStyle(8);
    g_pattern2_percentage->SetMarkerColor(3);
    g_pattern4_percentage->SetMarkerSize(2);
    g_pattern4_percentage->SetMarkerStyle(8);
    g_pattern4_percentage->SetMarkerColor(1);
    g_fullspan_percentage->SetTitle("");

    g_fullspan_percentage->Draw("AP");
    g_pattern1_percentage->Draw("P same");
    g_pattern2_percentage->Draw("P same");
    g_pattern4_percentage->Draw("P same");

    g_fullspan_percentage->GetXaxis()->SetLabelSize(0.04);
    g_fullspan_percentage->GetXaxis()->SetTitleSize(0.04);
    g_fullspan_percentage->GetYaxis()->SetLabelSize(0.04);
    g_fullspan_percentage->GetYaxis()->SetTitleSize(0.04);
    g_fullspan_percentage->GetXaxis()->SetTitle("Rigidity (GV)");
    g_fullspan_percentage->GetYaxis()->SetTitle("Patterns ratio in Inner Tracker Only");
    g_fullspan_percentage->GetYaxis()->SetRangeUser(0,0.6);

    TLegend * patternleg = new TLegend(0.2,0.45,0.33,0.63); //(xmin, ymin, xmax, ymax)
    patternleg->SetFillColor(0);
    patternleg->AddEntry(g_fullspan_percentage,"Pattern0","lp");
    patternleg->AddEntry(g_pattern1_percentage,"Pattern1","lp");
    patternleg->AddEntry(g_pattern2_percentage,"Pattern2","lp");
    patternleg->AddEntry(g_pattern4_percentage,"Pattern4","lp");
    gStyle->SetLegendTextSize(0.03);
    patternleg->Draw();

    allpattern_percentage.Update();
    allpattern_percentage.SaveAs(( std::string("/Time_Averaged_ratio/binmerge") + binmerge + std::string("/pattern_percentage_") + richcut + "_" + issversion + std::string(".pdf") ).c_str());

}


void Plot_Error_Compare(TH1D h_error_relative, TH1D h_Published_Statistic_Error_Relative, TH1D h_Published_Systematic_Error_Relative, std::string trackerpattern, std::string richcut, std::string issversion, std::string binmerge, std::string trdeff){

    TCanvas ccerror("ccerror","ccerror",1000,500);

    h_error_relative.Draw("");
    h_error_relative.SetTitle("");
    h_error_relative.GetYaxis()->SetRangeUser(0,5.6);
    h_error_relative.GetYaxis()->SetLimits(0,5.6);
    h_error_relative.GetXaxis()->SetTitle("Rigidity (GV)");
    h_error_relative.GetYaxis()->SetTitle("Relative Error (%)");
    h_error_relative.GetXaxis()->SetLabelSize(0.04);
    h_error_relative.GetXaxis()->SetTitleSize(0.04);
    h_error_relative.GetYaxis()->SetLabelSize(0.04);
    h_error_relative.GetYaxis()->SetTitleSize(0.04);
    h_error_relative.SetFillColor(0);
    h_error_relative.SetLineColor(1);
    h_Published_Statistic_Error_Relative.Draw("same");
    h_Published_Statistic_Error_Relative.SetLineColor(2);
    h_Published_Systematic_Error_Relative.Draw("same");
    h_Published_Systematic_Error_Relative.SetLineColor(4);

    TLegend * errorleg = new TLegend(0.4,0.75,0.7,0.85); //(xmin, ymin, xmax, ymax)
    errorleg->SetFillColor(0);
    if (issversion == "pass7.8"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (Full Range)","lp");}
    else if (issversion == "2016paper"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (Published Range)","lp");}
    else if (issversion == "PhysicsReport"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (PhysicsReport)","lp");}
    errorleg->AddEntry(&h_Published_Statistic_Error_Relative,"Published_Statistic_Error_Relative","lp");
    errorleg->AddEntry(&h_Published_Systematic_Error_Relative,"Published_Systematic_Error_Relative","lp");
    gStyle->SetLegendTextSize(0.03);
    errorleg->Draw();

    ccerror.Update();
    ccerror.SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/error_compare_") + trackerpattern + "_" + richcut + "_" + issversion + "binmerge"+ binmerge + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}


void Plot_AntiprotonMCTemplate(TH1F *template_antiproton, std::string trackerpattern, std::string richcut, std::string issversion, std::vector<double> binning, int index, std::string binmerge){
    TCanvas * c1 = new TCanvas;
    template_antiproton->Draw("");
    c1->SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/template_antiproton_") + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string(".pdf")).c_str());
}



void Plot_ISSPositiveData(TH1F *data_pass7_positive, std::string trackerpattern, std::string richcut, std::string issversion, std::vector<double> binning, int index, std::string binmerge, std::string trdeff){
    TCanvas * c11 = new TCanvas;
    data_pass7_positive->Draw("");
    c11->SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/template_protondata_") + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}

void Plot_ElectronDataTemplate(TH1F *template_electron_Data, std::string trackerpattern, std::string richcut, std::string issversion, std::vector<double> binning, int index, std::string binmerge, std::string trdeff){
    TCanvas * c12 = new TCanvas;
    template_electron_Data->Draw("");
    c12->SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/template_electron_Data_") + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}

void Plot_PionDataTemplate(TH1F *template_pion_Data, std::string trackerpattern, std::string richcut, std::string issversion, std::vector<double> binning, int index, std::string binmerge, std::string trdeff){
    TCanvas * c13 = new TCanvas;
    template_pion_Data->Draw("");
    c13->SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/template_pion_Data_") + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}


void Plot_ElectronMCTemplate(TH1F *template_electron, std::string trackerpattern, std::string richcut, std::string issversion, std::vector<double> binning, int index, std::string binmerge, std::string trdeff){
    TCanvas * c2 = new TCanvas;
    template_electron->Draw("");
    c2->SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/template_electron_") + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}


void Plot_ISSNegativeData(TH1F *data_pass7_negative, std::string trackerpattern, std::string richcut, std::string issversion, std::vector<double> binning, int index, std::string binmerge, std::string trdeff){
    TCanvas * c3 = new TCanvas;
    data_pass7_negative->Draw("");
    c3->SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/data_") + trackerpattern + "_" + richcut + "_" + issversion + "_" + doubleToString(binning[index]) + "_" + doubleToString(binning[index+stoi(binmerge)]) + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}


void Plot_ratio(TGraphErrors gRatioError, TGraph *g_ratio, std::string trackerpattern, std::string richcut, std::string issversion, std::string binmerge, std::string trdeff){
    TCanvas cc1("cc1","cc1",1000,500);
    gRatioError.Draw("AP");
    g_ratio->SetMarkerStyle(8);
    g_ratio->SetMarkerColor(2);
    g_ratio->Draw("P same");
    gPad->SetLogx();
    cc1.Update();
    cc1.SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/intermediate_ratio_") + trackerpattern + "_" + richcut + "_" + issversion + "binmerge"+ binmerge + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}


void Plot_ratio_with_effective_acceptance(TGraphErrors gRatioError, TGraphErrors *g_ratio_with_effective_acceptance, std::string issversion, std::string trackerpattern, std::string richcut, std::string binmerge, std::string trdeff){
    TCanvas cc12("cc12","cc12",1000,500);
    gRatioError.SetTitle("");
    gRatioError.Draw("AP");
    gRatioError.GetXaxis()->SetTitle("Rigidity (GV)");
    gRatioError.GetYaxis()->SetTitle("Antiproton ratio");
    gRatioError.GetXaxis()->SetLabelSize(0.04);
    gRatioError.GetXaxis()->SetTitleSize(0.04);
    gRatioError.GetYaxis()->SetLabelSize(0.04);
    gRatioError.GetYaxis()->SetTitleSize(0.04);
    g_ratio_with_effective_acceptance->SetMarkerStyle(8);
    g_ratio_with_effective_acceptance->SetMarkerColor(2);
    g_ratio_with_effective_acceptance->Draw("P same");
    gPad->SetLogx();
    gRatioError.GetXaxis()->SetMoreLogLabels();

    TLegend * leg = new TLegend(0.15,0.7,0.6,0.9); //(xmin, ymin, xmax, ymax)
    leg->SetFillColor(0);
    if (issversion == "pass7.8"){
        leg->AddEntry(g_ratio_with_effective_acceptance,"This analysis (Full Range) (Statistic Error only)","lp");}
    else if (issversion == "2016paper"){
        leg->AddEntry(g_ratio_with_effective_acceptance,"This analysis (Publisehd Range) (Statistic Error only)","lp");}
    else if (issversion == "PhysicsReport"){
        leg->AddEntry(g_ratio_with_effective_acceptance,"This analysis (PhysicsReport Range) (Statistic Error only)","lp");}
    leg->AddEntry(&gRatioError,"AMS PRL paper in 2016","lp");
    gStyle->SetLegendTextSize(0.03);
    leg->Draw();

    cc12.Update();
    cc12.SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/intermediate_ratio_with_effective_acceptance_") + trackerpattern + "_" + richcut + "_" + issversion + "binmerge"+ binmerge + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}

void Plot_Chi2(TGraph *g_chi2dof, std::string trackerpattern, std::string richcut, std::string issversion, std::string binmerge, std::string trdeff){
    TCanvas ccchi2("ccchi2","ccchi2",1000,500);

    g_chi2dof->Draw("AP");

    g_chi2dof->SetTitle("");
    g_chi2dof->GetXaxis()->SetTitle("Rigidity (GV)");
    g_chi2dof->GetYaxis()->SetTitle("Chi2/dof");

    g_chi2dof->GetXaxis()->SetTitleFont(62);
    g_chi2dof->GetXaxis()->SetTitleSize(0.045);
    g_chi2dof->GetYaxis()->SetTitleFont(62);
    g_chi2dof->GetYaxis()->SetTitleSize(0.045);

    g_chi2dof->GetXaxis()->SetLabelFont(62);
    g_chi2dof->GetXaxis()->SetLabelSize(0.05);
    g_chi2dof->GetYaxis()->SetLabelFont(62);
    g_chi2dof->GetYaxis()->SetLabelSize(0.05);
    g_chi2dof->SetMarkerStyle(15);
    g_chi2dof->SetMarkerColor(1);

    g_chi2dof->GetYaxis()->SetRangeUser(0,2);
    ccchi2.Update();

    ccchi2.SaveAs((std::string("Time_Averaged_ratio/binmerge") + binmerge + std::string("/chi2dof_") + trackerpattern + "_" + richcut + "_" + issversion + "binmerge"+ binmerge + std::string("_TRDEff_") + trdeff + std::string(".pdf")).c_str());
}


//// For Time Dependent ONLY
std::tuple<std::string, std::string, int> LoadTimeStapeSetings(std::string timemode){
    std::string bartals = "";
    std::string timeindex = "";
    int maxindex;
    if (timemode == "3BartalRotation")
        {
        bartals = std::to_string(AntiprotonNewBinning::mergestep_3B); 
        timeindex = "time.txt";
        maxindex = AntiprotonNewBinning::SplitTotal_3B;  
        }
    else if (timemode == "6BartalRotation")
        {
        bartals = std::to_string(AntiprotonNewBinning::mergestep_6B);
        timeindex = "time.txt";
        maxindex = AntiprotonNewBinning::SplitTotal_6B;
        }
    else if (timemode == "6months")
        {
        bartals = std::to_string(AntiprotonNewBinning::mergestep_6M);
        timeindex = "time_6months.txt";
        maxindex = AntiprotonNewBinning::SplitTotal_6M;
        }

    return {bartals, timeindex, maxindex};
}


void CreateTxtFiles(std::string timemode, std::vector<double> binning, int rigidityindex, std::string binmerge, std::string bartals, std::string trackerpattern, std::string richcut){

    chdir("./results");
    if (timemode == "3BartalRotation" || timemode == "6BartalRotation"){
        std::ofstream fitresults(std::string("fit_results_" + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + "BartalRotation.txt"), std::ios::trunc);
        std::ofstream fitresults_error(std::string("fit_results_error_" + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + "BartalRotation.txt"), std::ios::trunc);
        std::ofstream fitresults_proton_antiproton_error(std::string("fit_results_proton_antiproton_error_" + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + "BartalRotation.txt"), std::ios::trunc);
        std::ofstream time_stamp(std::string("timestamp_" + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + "BartalRotation.txt"), std::ios::trunc);
        std::ofstream fit_chi2dof(std::string("fit_chi2dof_" + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + "BartalRotation.txt"), std::ios::trunc);
        std::ofstream antiprotonnumber(std::string("antiprotonnumber_" + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + "BartalRotation.txt"), std::ios::trunc);
        std::ofstream protonnumber(std::string("protonnumber_" + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + "BartalRotation.txt"), std::ios::trunc);
    }
    else if(timemode == "6months"){
        std::ofstream fitresults(std::string("fit_results_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream fitresults_error(std::string("fit_results_error_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream fitresults_proton_antiproton_error(std::string("fit_results_proton_antiproton_error_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream time_stamp(std::string("timestamp_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream fit_chi2dof(std::string("fit_chi2dof_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream antiprotonnumber(std::string("antiprotonnumber_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream protonnumber(std::string("protonnumber_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
    }
    chdir("..");

}


std::tuple<std::vector<std::string>> LoadTimeStampBins(std::string timeindex){
    std::vector<std::string> v_time;
    std::string st;
    std::ifstream in2(timeindex);
    while (getline(in2,st))
        v_time.push_back(st);
    in2.close();

    return {v_time};
}


std::tuple<TH1F *>LoadISSPositiveData_TimeDependent(std::string lowpath, std::string binmerge, std::vector<double> binning, int rigidityindex, TCut PositiveCut, std::string timemode, int i, std::string bartals, int maxindex, double TrdLOW, double TrdHIGH, double TrdBINNUMBER){

    TChain *fpass7_positive = new TChain("tree");
    for (int index=i; index < i+stoi(bartals); index++){
        if (index==maxindex){
            std::cout<< "index="    << index    << std::endl;
            std::cout<< "maxindex=" << maxindex << std::endl;
            std::cout<< "break loop positive!"  << std::endl;
            break;}
        if (timemode == "3BartalRotation" || timemode == "6BartalRotation"){
            for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                std::string left = doubleToString(binning[rigidityindex+binindex]);
                std::string right = doubleToString(binning[rigidityindex+binindex+1]);
                if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
                    left = to_string_with_precision(binning[rigidityindex+binindex], 1);}
                if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
                    right = to_string_with_precision(binning[rigidityindex+binindex+1], 1);}
                fpass7_positive->AddFile(( lowpath + std::string("/rootfiles/positive/positive_") + left + std::string("_") + right + std::string("_") + std::to_string(index) + std::string(".root")).c_str());
            }
        }
        else if(timemode == "6months"){
            for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                std::string left = doubleToString(binning[rigidityindex+binindex]);
                std::string right = doubleToString(binning[rigidityindex+binindex+1]);
                if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
                    left = to_string_with_precision(binning[rigidityindex+binindex], 1);}
                if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
                    right = to_string_with_precision(binning[rigidityindex+binindex+1], 1);}
                fpass7_positive->AddFile(( lowpath + std::string("/rootfiles_6months/positive/positive_") + left + std::string("_") + right + std::string("_") + std::to_string(index) + std::string(".root")).c_str());
            }
        }
    }

    fpass7_positive->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_positive(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")")).c_str(), PositiveCut);
    TH1F *data_pass7_positive =  (TH1F*)gDirectory->Get("th1f_positive");
    data_pass7_positive->SetTitle("Antiproton (Proton data)");

    return {data_pass7_positive};
}


std::tuple<TH1F *, TH1F *, TH1F *>LoadISSNegativeData_TimeDependent(std::string lowpath, std::string binmerge, std::vector<double> binning, int rigidityindex, TCut NegativeCut, TCut ElectronTemplateDataCut, TCut PionTemplateDataCut, std::string timemode, int i, std::string bartals, int maxindex, double TrdLOW, double TrdHIGH, double TrdBINNUMBER){

    TChain *fpass7_negative = new TChain("tree");
    for (int index=i; index < i+stoi(bartals); index++){
        if (index==maxindex){
            std::cout<< "index="    << index    << std::endl;
            std::cout<< "maxindex=" << maxindex << std::endl;
            std::cout<< "break loop negative!"  << std::endl;
            break;}

        if (timemode == "3BartalRotation" || timemode == "6BartalRotation"){
            for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                std::string left = doubleToString(binning[rigidityindex+binindex]);
                std::string right = doubleToString(binning[rigidityindex+binindex+1]);
                if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
                    left = to_string_with_precision(binning[rigidityindex+binindex], 1);}
                if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
                    right = to_string_with_precision(binning[rigidityindex+binindex+1], 1);}
                fpass7_negative->AddFile(( lowpath + std::string("/rootfiles/negative/negative_") + left + std::string("_") + right + std::string("_") + std::to_string(index) + std::string(".root")).c_str());
            }
        }
        else if(timemode == "6months"){
            for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                std::string left = doubleToString(binning[rigidityindex+binindex]);
                std::string right = doubleToString(binning[rigidityindex+binindex+1]);
                if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
                    left = to_string_with_precision(binning[rigidityindex+binindex], 1);}
                if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
                    right = to_string_with_precision(binning[rigidityindex+binindex+1], 1);}
                fpass7_negative->AddFile(( lowpath + std::string("/rootfiles_6months/negative/negative_") + left + std::string("_") + right + std::string("_") + std::to_string(index) + std::string(".root")).c_str());
            }
        }
    }

    fpass7_negative->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_negative(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")")).c_str(), NegativeCut);
    TH1F *data_pass7_negative =  (TH1F*)gDirectory->Get("th1f_negative");


    fpass7_negative->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_ElectronDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")")).c_str(), ElectronTemplateDataCut);
    TH1F *template_electron_Data =  (TH1F*)gDirectory->Get("th1f_ElectronDataTemplate");
    template_electron_Data->SetTitle("Electron");

    fpass7_negative->Draw( (std::string("(-1)*TrdLogLikelihoodRatioElectronProtonTracker>>th1f_PionDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(")")).c_str(), PionTemplateDataCut);
    TH1F *template_pion_Data =  (TH1F*)gDirectory->Get("th1f_PionDataTemplate");
    template_pion_Data->SetTitle("Pion");

    return {data_pass7_negative, template_electron_Data, template_pion_Data};
}


void SaveResult(std::string bartals, std::vector<double> binning, int rigidityindex, std::string binmerge, std::vector<double> result, double ProtonNumber, std::vector<double> ResultError, double CHI2dof, std::string trackerpattern, std::string richcut){

    if (stoi(bartals) == 3 || stoi(bartals) == 6){
        chdir("./results");
        std::ofstream fitresults(std::string("fit_results_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
        if (fitresults.is_open())
              {
              fitresults << result[0]/ProtonNumber << std::endl;
              fitresults.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream fitresultserror(std::string("fit_results_error_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") +  doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
        if (fitresultserror.is_open())
              {
              fitresultserror << ResultError[0]/ProtonNumber << std::endl;
              fitresultserror.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream fitresults_p_antip_error(std::string("fit_results_proton_antiproton_error_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") +  doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
        if (fitresults_p_antip_error.is_open())
              {
              fitresults_p_antip_error <<  ProtonNumber / pow(result[0],2) * ResultError[0] << std::endl; // only the second term.
              fitresults_p_antip_error.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream fitresultschi2dof(std::string("fit_chi2dof_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
        if (fitresultschi2dof.is_open())
              {
              fitresultschi2dof << CHI2dof << std::endl;
              fitresultschi2dof.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream antiprotonnumberstream(std::string("antiprotonnumber_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
        if (antiprotonnumberstream.is_open())
              {
              antiprotonnumberstream << result[0] << std::endl;
              antiprotonnumberstream.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream protonnumberstream(std::string("protonnumber_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
        if (protonnumberstream.is_open())
              {
              protonnumberstream << ProtonNumber << std::endl;
              protonnumberstream.close();
              }
        chdir("..");
    }
    else if (stoi(bartals) == 1)
    {
        chdir("./results");
        std::ofstream fitresults(std::string("fit_results_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") +  std::string("6months.txt"), std::ios::app);
        if (fitresults.is_open())
              {
              fitresults << result[0]/ProtonNumber << std::endl;
              fitresults.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream fitresultserror(std::string("fit_results_error_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
        if (fitresultserror.is_open())
              {
              fitresultserror << ResultError[0]/ProtonNumber << std::endl;
              fitresultserror.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream fitresults_p_antip_error(std::string("fit_results_proton_antiproton_error_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") +  doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("6months.txt"), std::ios::app);
        if (fitresults_p_antip_error.is_open())
              {
              fitresults_p_antip_error <<  ProtonNumber / pow(result[0],2) * ResultError[0] << std::endl; // only the second term.
              fitresults_p_antip_error.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream fitresultschi2dof(std::string("fit_chi2dof_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
        if (fitresultschi2dof.is_open())
              {
              fitresultschi2dof << CHI2dof << std::endl;
              fitresultschi2dof.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream antiprotonnumberstream(std::string("antiprotonnumber_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
        if (antiprotonnumberstream.is_open())
              {
              antiprotonnumberstream << result[0] << std::endl;
              antiprotonnumberstream.close();
              }
        chdir("..");

        chdir("./results");
        std::ofstream protonnumberstream(std::string("protonnumber_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
        if (protonnumberstream.is_open())
              {
              protonnumberstream << ProtonNumber << std::endl;
              protonnumberstream.close();
              }
        chdir("..");

        /*
        if (stoi(bartals) == 3)
            {
            chdir("./results");
            fstream timestamp(std::string("timestamp_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"),std::ios::app);
            if (timestamp.is_open())
                  {
                  timestamp << stoi(v_time.at(i+1))+1166400 << std::endl;
                  timestamp.close();
                  }
            chdir("..");
            }
        else if (stoi(bartals) == 6)
            {
            chdir("./results");
            fstream timestamp(std::string("timestamp_") + trackerpattern + "_" + richcut + "_" + doubleToString(binning[rigidityindex]) + std::string("_") + doubleToString(binning[rigidityindex+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"),std::ios::app);
            if (timestamp.is_open())
                  {
                  timestamp << stoi(v_time.at(i+3)) << std::endl;
                  timestamp.close();
                  }
            chdir("..");
            }
        else if (stoi(bartals) == 1)
            {
            chdir("./results");
            fstream timestamp(std::string("timestamp_6months" + trackerpattern + "_" + richcut + ".txt"),std::ios::app);
            if (timestamp.is_open())
                  {
                  timestamp << stoi(v_time.at(i-1))+7776000 << std::endl;
                  timestamp.close();
                  }
            chdir("..");
            */

    }

}


std::tuple<TH1F *, TH1F *, TH1F *, TH1F *> LoadTemplates_TimeDependent( std::string binmerge, std::vector<double> binning, int rigidityindex, std::string timemode, int TimeIndex, std::string ParametrilizedMode, int GeneraterNumberIndex, std::string SysErrEffStudyMode, std::string trdeff, std::string FullRange){

    std::string left  = doubleToString(binning[rigidityindex]  );
    std::string right = doubleToString(binning[rigidityindex+stoi(binmerge)]);
    if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
        left  = to_string_with_precision(binning[rigidityindex]  , 1);}
    if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
        right = to_string_with_precision(binning[rigidityindex+stoi(binmerge)], 1);}

    std::string nameindex;
    std::string templateindex;
    if (ParametrilizedMode == "Yes"){
        nameindex     = "_new";
        templateindex = std::string("_RandomIndex_") + std::to_string(GeneraterNumberIndex);
    }
    else if (ParametrilizedMode == "No"){
        nameindex     = "";
        templateindex = "";
    }

    std::string filename;
    if (timemode == "3BartalRotation"){
        filename = "rootfiles_3BartalRotation_template";}
    else if (timemode == "6BartalRotation"){
        filename = "rootfiles_6BartalRotation_template";}
    else if (timemode == "6months"){
        filename = "rootfiles_6months_template";}

    std::string TemplateFileName;  
    if ( FullRange == "Yes"){
        TemplateFileName = "/TimeDependentTemplatesAndData_fullRange_";
    }
    else{
        TemplateFileName = "/TimeDependentTemplatesAndData_";
    }

    TFile *f_histogram = new TFile( ( filename + TemplateFileName + left + std::string("_") + right + std::string("_") + std::to_string(TimeIndex) + nameindex + std::string(".root")).c_str(), "READ");
    TH1F *data_pass7_negative    = (TH1F*) f_histogram->Get( ((std::string("data_pass7_negative")    + std::string("_TRDeff_") + trdeff + templateindex).c_str() ));
    TH1F *template_electron_Data = (TH1F*) f_histogram->Get( ((std::string("template_electron_Data") + std::string("_TRDeff_") + trdeff + templateindex).c_str() ));
    TH1F *template_pion_Data     = (TH1F*) f_histogram->Get( ((std::string("template_pion_Data")     + std::string("_TRDeff_") + trdeff + templateindex).c_str() ));
    TH1F *data_pass7_positive    = (TH1F*) f_histogram->Get( ((std::string("data_pass7_positive")    + std::string("_TRDeff_") + trdeff                ).c_str() ));

    data_pass7_negative   ->SetDirectory(0);
    template_electron_Data->SetDirectory(0);
    template_pion_Data    ->SetDirectory(0);
    data_pass7_positive   ->SetDirectory(0);

    f_histogram->Close();


    return {data_pass7_negative, template_electron_Data, template_pion_Data, data_pass7_positive};

}




