#include "TemplateFitterforLowEnergy2D.hh"

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


std::tuple<TGraphErrors *> LoadEffectiveAcceptance(int number){
    TFile *f_eff_antiproton = new TFile("EffectiveAcceptance_B1042_antipr.pl1.1800_7.6_all.root");
    TFile *f_eff_proton     = new TFile("EffectiveAcceptance_B1042_pr.pl1.1800_7.6_all.root");
    TGraphAsymmErrors *Acceptance_antiproton = (TGraphAsymmErrors*)f_eff_antiproton->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");
    TGraphAsymmErrors *Acceptance_proton     = (TGraphAsymmErrors*)f_eff_proton->Get("QualityCuts/effectiveAcceptanceAfterAllCuts");

    Double_t xa[number],Aa[number],ea[number],xp[number],Ap[number],ep[number];

    TGraphErrors *Effective_Acceptance = new TGraphErrors(number);
    for (int p = 0; p < number; p++) {
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


std::tuple<std::vector<double>> LoadTRDLogLikelihoodCut(std::string TRDefficiency, std::string binmerge){
    std::vector<double> v_TRDcut;
    std::ifstream TRDfile;
    TRDfile.open( std::string("TRDLogLikelihood_CutValue_eff_") + TRDefficiency + std::string("_") + binmerge + std::string(".txt") );
    std::string trdcutvalue;
    assert(v_TRDcut.empty());
    while (getline(TRDfile, trdcutvalue)) {
    v_TRDcut.push_back(stod(trdcutvalue));
    }
    TRDfile.close();
    return {v_TRDcut};
}


std::tuple<std::vector<double>> LoadTOFBetaCut(std::string TOFefficiency, std::string binmerge){
    std::vector<double> v_TOFcut;
    std::ifstream TOFfile;
    TOFfile.open( std::string("TOFBeta_CutValue_eff_") + TOFefficiency + std::string("_") + binmerge + std::string(".txt") );
    std::string tofcutvalue;
    assert(v_TOFcut.empty());
    while (getline(TOFfile, tofcutvalue)) {
    v_TOFcut.push_back(stod(tofcutvalue));
    }
    TOFfile.close();
    return {v_TOFcut};
}


std::tuple<TGraphErrors *, TH1D, TH1D> LoadPublishedResult(std::vector<double> bincenter, std::string rigidity_start, std::string rigidity_end, std::vector<double> subbincenter, std::vector<double> subbinedge){
    int publishedresult_showpoint = 21;  //16points:1-17:1.08-5.635, 21points:1.08-8.87GV, 24points:1.08-11.5, 29points:1.08-17.3.
    std::vector<double> lowbincenter(bincenter.begin()+1, bincenter.begin()+publishedresult_showpoint+1); //first bin (0.8-1.0GV) is used for underflow.
    std::vector<double> published(AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioPRL.begin()+publishedresult_showpoint);
    std::vector<double> publishederror(AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin(), AntiprotonNewBinning::AntiprotonResults::PublishedRatioErrorPRL.begin()+publishedresult_showpoint);
    std::vector<double> v_StatisticError_relative(AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.begin()+stoi(rigidity_start), AntiprotonNewBinning::AntiprotonResults::PublishedRatioStatisticRelativeErrorPRL.begin()+stoi(rigidity_end));
    std::vector<double> v_SystematicError_relative(AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.begin()+stoi(rigidity_start), AntiprotonNewBinning::AntiprotonResults::PublishedRatioSystematicRelativeErrorPRL.begin()+stoi(rigidity_end));
    for (long unsigned int i=0; i<published.size(); i++){
      published.at(i) = published.at(i)*100000;
      publishederror.at(i) = publishederror.at(i)*100000;
    }
    TGraphErrors *g_ratio_published = new TGraphErrors(publishedresult_showpoint, lowbincenter.data(), published.data(), 0, publishederror.data());

    TGraph *g_Published_Statistic_Error_Relative = new TGraph(subbincenter.size(), subbincenter.data(), v_StatisticError_relative.data());
    TH1D h_Published_Statistic_Error_Relative = TH1D("", "", subbincenter.size(), subbinedge.data());
    Utilities::ConvertToHistogram(g_Published_Statistic_Error_Relative, h_Published_Statistic_Error_Relative);

    TGraph *g_Published_Systematic_Error_Relative = new TGraph(subbincenter.size(), subbincenter.data(), v_SystematicError_relative.data());
    TH1D h_Published_Systematic_Error_Relative = TH1D("", "", subbincenter.size(), subbinedge.data());
    Utilities::ConvertToHistogram(g_Published_Systematic_Error_Relative, h_Published_Systematic_Error_Relative);

    return {g_ratio_published, h_Published_Statistic_Error_Relative, h_Published_Systematic_Error_Relative};
}


std::tuple<TH2F *, TH2F *, TH2F *>LoadISSPositiveData(std::string binmerge, std::vector<double> binning, int index, TCut TOFBETACut, TCut TrdCut, TCut PositiveCut, TCut PositiveCut_test, std::string issversion, std::string TestMode, int TrdBINNUMBER, double TrdLOW, double TrdHIGH, int RECITOFBETANUMBER, double RECITOFBETALOW, double RECITOFBETAHIGH){

    TChain *fpass7_positive = new TChain("AntiprotonLowEnergyTree");
    if (issversion == "pass7.8"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            if (TestMode == "Yes"){
                fpass7_positive->AddFile((std::string("B1130_pass7_7.8_all_Tree_positive_") + left + std::string("_") + right + std::string("_test.root")).c_str());}
            else{
                fpass7_positive->AddFile((std::string("B1130_pass7_7.8_all_Tree_positive_") + left + std::string("_") + right + std::string(".root")).c_str());}
            }
   }
    else if (issversion == "2016paper"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            if (TestMode == "Yes"){
                fpass7_positive->AddFile((std::string("B1130_pass7_7.8_all_Tree_positive_May2015_") + left + std::string("_") + right + std::string("_test.root")).c_str());}
            else{
                fpass7_positive->AddFile((std::string("B1130_pass7_7.8_all_Tree_positive_May2015_") + left + std::string("_") + right + std::string(".root")).c_str());}
            }
    }
    else if (issversion == "PhysicsReport"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            if (TestMode == "Yes"){
                fpass7_positive->AddFile((std::string("B1130_pass7_7.8_all_Tree_positive_Nov2017_") + left + std::string("_") + right + std::string("_test.root")).c_str());}
            else{
                fpass7_positive->AddFile((std::string("B1130_pass7_7.8_all_Tree_positive_Nov2017_") + left + std::string("_") + right + std::string(".root")).c_str());}
            }
    }

    // Template Antiproton (from ISS positve data) (AntiprotonCut, PositiveCut)
    std::cout<< "Making Antiproton template histgram now." << std::endl;
    fpass7_positive->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_positive_template(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(",") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), TOFBETACut && TrdCut && PositiveCut);
    TH2F *TofTRD_data_pass7_positive_template = (TH2F*)gDirectory->Get("th2f_positive_template");
    TofTRD_data_pass7_positive_template->SetTitle("Antiproton");

    // ISS Positive data (Proton number)
    std::cout<< "Making Proton histgram now." << std::endl;
    fpass7_positive->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_positive_data(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(",") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), TOFBETACut && TrdCut && PositiveCut);
    fpass7_positive->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_positive_data_test(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(",") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), TOFBETACut && TrdCut && PositiveCut_test);
    TH2F *TofTRD_data_pass7_positive_data = (TH2F*)gDirectory->Get("th2f_positive_data");
    TH2F *TofTRD_data_pass7_positive_data_test = (TH2F*)gDirectory->Get("th2f_positive_data_test");
    TofTRD_data_pass7_positive_data->SetTitle("ISSpositive_data");    

    fpass7_positive->Reset();

    return {TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data, TofTRD_data_pass7_positive_data_test};
}


std::tuple<TH2F *, TH2F *, TH2F *>LoadISSNegativeData(std::string binmerge, std::vector<double> binning, int index, TCut TOFBETACut, TCut TrdCut, TCut NegativeCut, TCut ElectronTemplateDataCut, TCut PionTemplateDataCut, std::string issversion, std::string TestMode, int TrdBINNUMBER, double TrdLOW, double TrdHIGH, int RECITOFBETANUMBER, double RECITOFBETALOW, double RECITOFBETAHIGH){

    TChain *fpass7_negative = new TChain("AntiprotonLowEnergyTree");
    if (issversion == "pass7.8"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            fpass7_negative->AddFile((std::string("B1130_pass7_7.8_all_Tree_negative_") + left + std::string("_") + right + std::string(".root")).c_str());}
    }
    else if (issversion == "2016paper"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            fpass7_negative->AddFile((std::string("B1130_pass7_7.8_all_Tree_negative_May2015_") + left + std::string("_") + right + std::string(".root")).c_str());}
    }
    else if (issversion == "PhysicsReport"){
        for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
            std::string left = doubleToString(binning[index+binindex]);
            std::string right = doubleToString(binning[index+binindex+1]);
            if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                left = to_string_with_precision(binning[index+binindex], 1);}
            if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                right = to_string_with_precision(binning[index+binindex+1], 1);}
            fpass7_negative->AddFile((std::string("B1130_pass7_7.8_all_Tree_negative_Nov2017_") + left + std::string("_") + right + std::string(".root")).c_str());}
    }

    // ISS Negative data
    std::cout<< "Making ISS negaitve data histgram now." <<std::endl;
    std::cout<< "Total events before cuts: " << fpass7_negative->GetEntries() <<std::endl;
    fpass7_negative->Draw((std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_negative(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), TOFBETACut && TrdCut && NegativeCut);
    TH2F *TofTRD_data_pass7_negative = (TH2F*)gDirectory->Get("th2f_negative");
    TofTRD_data_pass7_negative->SetTitle("ISSnegative");
    std::cout<< "Total events after cuts (X):" << TofTRD_data_pass7_negative->ProjectionX()->GetEntries() <<std::endl;
    std::cout<< "Total events after cuts (Y):" << TofTRD_data_pass7_negative->ProjectionY()->GetEntries() <<std::endl;

    // Template Electron (from data)
    std::cout<< "Making Electron template histgram now." <<std::endl;
    fpass7_negative->Draw((std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_ElectronDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), TOFBETACut && TrdCut && ElectronTemplateDataCut);
    TH2F *TofTRD_template_ElectronData = (TH2F*)gDirectory->Get("th2f_ElectronDataTemplate");
    TofTRD_template_ElectronData->SetTitle("Electron");
    
    // Template Pion (from data)
    std::cout<< "Making Pion template histgram now." <<std::endl;
    fpass7_negative->Draw((std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_PionDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), TOFBETACut && TrdCut && PionTemplateDataCut);
    TH2F *TofTRD_template_PionData = (TH2F*)gDirectory->Get("th2f_PionDataTemplate");
    TofTRD_template_PionData->SetTitle("Pion");

    fpass7_negative->Reset();

    return {TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData};
}


std::tuple<TH2F *, TH2F *, TH2F *, TH2F *, TH2F *, TH2F *, double>LoadTemplates(std::string binmerge, std::vector<double> binning, int index, std::string issversion, std::string ParametrilizedMode, int numberindex, std::string trdeff, std::string tofeff, std::string FullRange, double RECITOFBETALOW = 0, double RECITOFBETAHIGH = 0, double RECITOFBETANUMBER = 0, double TrdLOW = 0, double TrdHIGH = 0, double TrdBINNUMBER = 0){

    std::string left;
    std::string right;
    double ProtonNumber;

    left  = doubleToString(binning[index]);
    right = doubleToString(binning[index+stoi(binmerge)]);
    if ( left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) ){
        left = to_string_with_precision(binning[index], 1);}
    if ( right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) ){
        right = to_string_with_precision(binning[index+stoi(binmerge)], 1);}

    //std::cout<< "left:" << left << ", right:" << right << std::endl;

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

    TFile *f_histogram = new TFile( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/") + TemplateRootName + left + std::string("_") + right + std::string("_") + issversion + nameindex + std::string(".root")).c_str(), "READ" );

    std::string TemplateName;
    if (FullRange == "No"){
        TemplateName = ( std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff);
    }
    else{
        TemplateName = ( std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + std::string("_RECITOFBETALOW_") + to_string_with_precision(RECITOFBETALOW, 2) + std::string("_RECITOFBETAHIGH_") + to_string_with_precision(RECITOFBETAHIGH, 2) + std::string("_RECITOFBETANUMBER_") + to_string_with_precision(RECITOFBETANUMBER, 0) + std::string("_TrdLOW_") + to_string_with_precision(TrdLOW, 2) + std::string("_TrdHIGH_") + to_string_with_precision(TrdHIGH, 2) + std::string("_TrdBINNUMBER_") + to_string_with_precision(TrdBINNUMBER, 0) );
    }

    TH2F *TofTRD_data_pass7_positive_template  = (TH2F*) f_histogram->Get( ((std::string("TofTRD_data_pass7_positive_template")  + TemplateName + templateindex).c_str() ));
    TH2F *TofTRD_template_ElectronData         = (TH2F*) f_histogram->Get( ((std::string("TofTRD_template_ElectronData")         + TemplateName + templateindex).c_str() ));
    TH2F *TofTRD_template_PionData             = (TH2F*) f_histogram->Get( ((std::string("TofTRD_template_PionData")             + TemplateName + templateindex).c_str() ));
    TH2F *TofTRD_data_pass7_positive_data      = (TH2F*) f_histogram->Get( ((std::string("TofTRD_data_pass7_positive_data")      + TemplateName ).c_str() ));
    TH2F *TofTRD_data_pass7_positive_data_test = (TH2F*) f_histogram->Get( ((std::string("TofTRD_data_pass7_positive_data_test") + TemplateName ).c_str() ));
    TH2F *TofTRD_data_pass7_negative           = (TH2F*) f_histogram->Get( ((std::string("TofTRD_data_pass7_negative")           + TemplateName ).c_str() ));
    //f_histogram->Close();
   

    ProtonNumber = TofTRD_data_pass7_positive_data->GetEntries();

    return {TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data, TofTRD_data_pass7_positive_data_test, TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData, ProtonNumber};

}


void PrintFitResult(std::vector<double> result_tof, TH2F *TofTRD_data_pass7_positive_data, double RescaleFactor, TH2F *TofTRD_data_pass7_positive_data_test, double CHI2dof_tof, double Chi2_tof, int NDF_tof){

    std::cout << "tof fit result:"                  << result_tof                                                         << std::endl;
    std::cout << "proton total number after cut: "  << TofTRD_data_pass7_positive_data->GetEntries() * RescaleFactor      << std::endl; 
    std::cout << "proton total number before cut: " << TofTRD_data_pass7_positive_data_test->GetEntries() * RescaleFactor << std::endl; 
    std::cout << "proton total cuts eff: "          << TofTRD_data_pass7_positive_data_test->GetEntries()/TofTRD_data_pass7_positive_data->GetEntries() * 100 << " %." << std::endl; 
    std::cout << "Chi2_tof:"                        << Chi2_tof                                                           << std::endl;
    std::cout << "NDF_tof:"                         << NDF_tof                                                            << std::endl;
    std::cout << "Chi2/dof_tof:"                    << CHI2dof_tof                                                        << std::endl;
    std::cout << "ratio from 1./tofbeta fit: "      << result_tof[0] / (TofTRD_data_pass7_positive_data->GetEntries() * RescaleFactor) * 100000  << std::endl;
}


void Plot_TemplatePion_Data(std::string binmerge, std::vector<double> binning, int index, std::string issversion, TH2F *TofTRD_template_PionData, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency, std::string AnalysisMode, int TimeIndex){
    TCanvas * c112 = new TCanvas;
    TofTRD_template_PionData->Draw("COLZ");
    if (AnalysisMode == "TimeAveraged"){
        if (ParametrilizedMode == "Yes"){
            c112->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_template_PionData_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str());
        }
        else if (ParametrilizedMode == "No"){
            c112->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_template_PionData_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
        }
    }
    else if (AnalysisMode == "TimeDependent"){
        c112->SaveAs( (std::string("results/fitplot/TofTRD_template_PionData_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::to_string(TimeIndex) + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
    }
}


void Plot_TemplateElectron_Data(std::string binmerge, std::vector<double> binning, int index, std::string issversion, TH2F *TofTRD_template_ElectronData, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency, std::string AnalysisMode, int TimeIndex){
    TCanvas * c113 = new TCanvas;
    TofTRD_template_ElectronData->Draw("COLZ");
    if (AnalysisMode == "TimeAveraged"){
        if (ParametrilizedMode == "Yes"){
            c113->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_template_ElectronData_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str());
        }
        else if (ParametrilizedMode == "No"){
            c113->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_template_ElectronData_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") +issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
        }
    }
    else if (AnalysisMode == "TimeDependent"){
        c113->SaveAs( (std::string("results/fitplot/TofTRD_template_ElectronData_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::to_string(TimeIndex) + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
    }
}


void Plot_ISSNegative(std::string binmerge, std::vector<double> binning, int index, std::string issversion, TH2F *TofTRD_data_pass7_negative, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency, std::string AnalysisMode, int TimeIndex){
    TCanvas * c12 = new TCanvas;
    TofTRD_data_pass7_negative->Draw("COLZ");
    if (AnalysisMode == "TimeAveraged"){
        if (ParametrilizedMode == "Yes"){
            c12->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_data_pass7_negative_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str());
        }
        else if (ParametrilizedMode == "No"){
            c12->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_data_pass7_negative_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
        }
    }
    else if (AnalysisMode == "TimeDependent"){
        c12->SaveAs( (std::string("results/fitplot/TofTRD_data_pass7_negative_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::to_string(TimeIndex) + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
    }

}


void Plot_TemplateAntiproton_Data(std::string binmerge, std::vector<double> binning, int index, std::string issversion, TH2F *TofTRD_data_pass7_positive_template, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency, std::string AnalysisMode, int TimeIndex){
    TCanvas * c14 = new TCanvas;
    TofTRD_data_pass7_positive_template->Draw("COLZ");
    if (AnalysisMode == "TimeAveraged"){
        if (ParametrilizedMode == "Yes"){
            c14->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_data_pass7_positive_template_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str());
        }
        else if (ParametrilizedMode == "No"){
            c14->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TofTRD_data_pass7_positive_template_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
        }
    }
    else if (AnalysisMode == "TimeDependent"){
        c14->SaveAs( (std::string("results/fitplot/TofTRD_data_pass7_positive_template_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::to_string(TimeIndex) + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
    }

}


/*
void Plot_FitResult2D(std::string binmerge, std::vector<double> binning, int index, std::string issversion, MYUtilities::TemplateFitter2D templateFitter_toftrd_2D){
    TCanvas * c8 = templateFitter_toftrd_2D.CreateResultDrawing("Fit_Result_tof",1000,500);
    c8->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/FitResult_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string(".pdf")).c_str());
}
*/

/*
void Plot_FitResult_TrdProjection(std::string binmerge, std::vector<double> binning, int index, std::string issversion, MYUtilities::TemplateFitter2D templateFitter_toftrd_2D, double RECITOFBETALOW, double RECITOFBETAHIGH){
    TCanvas * c81 = templateFitter_toftrd_2D.CreateResultDrawingXprojection("TrdLikelihood_X",1000,800, RECITOFBETALOW, RECITOFBETAHIGH);
    //TCanvas * c81 = templateFitter_toftrd_2D.CreateResultDrawingXprojection("TrdLikelihood_X",1000,800, -0.03, RECITOFBETAHIGH); // One bin to show for 2.2-2.4GV Generam meeting 06.2021
    c81->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TrdLikelihood_X_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string(".pdf")).c_str());

    TCanvas * c812 = templateFitter_toftrd_2D.CreateResultDrawingXprojection_LogY("TrdLikelihood_X",1000,800, RECITOFBETALOW, RECITOFBETAHIGH);
    //TCanvas * c812 = templateFitter_toftrd_2D.CreateResultDrawingXprojection_LogY("TrdLikelihood_X",1000,800, -0.03, RECITOFBETAHIGH); // One bin to show for 2.2-2.4GV Generam meeting 06.2021
    c812->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/TrdLikelihood_X_tof_LogY_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string(".pdf")).c_str());
}
*/

/*
void Plot_FitResult_TofBetaProjection(std::string binmerge, std::vector<double> binning, int index, std::string issversion, MYUtilities::TemplateFitter2D templateFitter_toftrd_2D, double TrdLOW, double TrdHIGH){
    TCanvas * c82 = templateFitter_toftrd_2D.CreateResultDrawingYprojection("1/TofBeta_Y",1000,800, TrdLOW, TrdHIGH);
    //TCanvas * c82 = templateFitter_toftrd_2D.CreateResultDrawingYprojection("1/TofBeta_Y",1000,800, 0.92, TrdHIGH); // One bin to show for 2.2-2.4GV Generam meeting 06.2021
    c82->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/1_TofBeta_Y_tof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string(".pdf")).c_str());

    TCanvas * c822 = templateFitter_toftrd_2D.CreateResultDrawingYprojection_LogY("1/TofBeta_Y",1000,800, TrdLOW, TrdHIGH);
    //TCanvas * c822 = templateFitter_toftrd_2D.CreateResultDrawingYprojection_LogY("1/TofBeta_Y",1000,800, 0.92, TrdHIGH);  // One bin to show for 2.2-2.4GV Generam meeting 06.2021    
    c822->SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/tof/1_TofBeta_Y_tof_LogY_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + issversion + std::string(".pdf")).c_str());
}
*/




void Plot_Low_Ratio_Tof(TGraphErrors *g_ratio_published, TGraph *g_ratio_tof, std::string binmerge, std::string issversion, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency){
    TCanvas cc3("cc3","cc3",1000,500);

    g_ratio_published->Draw("AP *");
    g_ratio_tof->Draw("same P");

    g_ratio_tof->SetMarkerStyle(15);
    g_ratio_tof->SetMarkerColor(2);

    g_ratio_published->SetTitle("");
    g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
    g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");

    cc3.Update();
    if (ParametrilizedMode == "No"){
        cc3.SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/low_ratio_tof_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str() );
    }
    else if (ParametrilizedMode == "Yes"){
        cc3.SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/low_ratio_tof_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str() );
    }
}


void Plot_Low_Ratio_Tof_WithAcceptanceCorrection(TGraphErrors *g_ratio_published, TGraphErrors *g_ratio_tof_with_effective, std::string binmerge, std::string issversion, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency){
    TCanvas cc33("cc33","cc33",1000,500);

    g_ratio_published->Draw("AP *");
    g_ratio_tof_with_effective->Draw("same P");
    //g_IntermediateResult->Draw("same P");

    g_ratio_tof_with_effective->SetMarkerStyle(15);
    g_ratio_tof_with_effective->SetMarkerColor(2);
    //g_IntermediateResult->SetMarkerStyle(15);
    //g_IntermediateResult->SetMarkerColor(4);

    //g_ratio_published->GetYaxis()->SetRangeUser(0,19);
    g_ratio_published->GetYaxis()->SetLimits(0,19);
    //g_ratio_published->GetYaxis()->SetRange(0,19);
    g_ratio_published->SetTitle("");
    g_ratio_published->GetXaxis()->SetTitle("Rigidity (GV)");
    g_ratio_published->GetYaxis()->SetTitle("ratio (*10^5)");

    cc33.Update();
    if (ParametrilizedMode == "No"){
        cc33.SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/low_ratio_tof_with_effective_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str() );
    }
    else if (ParametrilizedMode == "Yes"){
        cc33.SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/low_ratio_tof_with_effective_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str() );
    }
}


void Plot_Chi2(TGraph *g_chi2_tof, std::string binmerge, std::string issversion, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency){
    TCanvas c1("c1","c1",1000,500);

    g_chi2_tof->Draw("AP");

    g_chi2_tof->SetTitle("");
    g_chi2_tof->GetXaxis()->SetTitle("Rigidity (GV)");
    g_chi2_tof->GetYaxis()->SetTitle("Chi2/dof");

    g_chi2_tof->GetXaxis()->SetTitleFont(62);
    g_chi2_tof->GetXaxis()->SetTitleSize(0.045);
    g_chi2_tof->GetYaxis()->SetTitleFont(62);
    g_chi2_tof->GetYaxis()->SetTitleSize(0.045);

    g_chi2_tof->GetXaxis()->SetLabelFont(62);
    g_chi2_tof->GetXaxis()->SetLabelSize(0.05);
    g_chi2_tof->GetYaxis()->SetLabelFont(62);
    g_chi2_tof->GetYaxis()->SetLabelSize(0.05);
    g_chi2_tof->SetMarkerStyle(15);
    g_chi2_tof->SetMarkerColor(1);

    g_chi2_tof->GetYaxis()->SetRangeUser(0,2);
    c1.Update();

    if (ParametrilizedMode == "No"){
        c1.SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/chi2_tof_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str() );
    }
    else if (ParametrilizedMode == "Yes"){
        c1.SaveAs( (std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/chi2_tof_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str() );
    }

}


void Plot_error_compare(TH1D h_error_relative, TH1D h_Published_Statistic_Error_Relative, TH1D h_Published_Systematic_Error_Relative, std::string binmerge, std::string issversion, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency){
    TCanvas ccerror("ccerror","ccerror",1000,500);

    h_error_relative.Draw("");
    h_Published_Statistic_Error_Relative.Draw("same");
    h_Published_Systematic_Error_Relative.Draw("same");

    h_error_relative.SetFillColor(0);
    h_error_relative.SetLineColor(1);
    h_Published_Statistic_Error_Relative.SetLineColor(2);
    h_Published_Systematic_Error_Relative.SetLineColor(4);

    h_error_relative.SetTitle("");
    h_error_relative.GetYaxis()->SetRangeUser(0,12);
    h_error_relative.GetXaxis()->SetTitle("Rigidity (GV)");
    h_error_relative.GetYaxis()->SetTitle("Relative Error (%)");
    h_error_relative.GetXaxis()->SetLabelSize(0.04);
    h_error_relative.GetXaxis()->SetTitleSize(0.04);
    h_error_relative.GetYaxis()->SetLabelSize(0.04);
    h_error_relative.GetYaxis()->SetTitleSize(0.04);

    TLegend * errorleg = new TLegend(0.4,0.75,0.7,0.85); //(xmin, ymin, xmax, ymax)
    errorleg->SetFillColor(0);
    if (issversion == "pass7.8"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (Full Range)","lp");}
    else if (issversion == "2016paper"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (2016paper Range)","lp");}
    else if (issversion == "PhysicsReport"){
        errorleg->AddEntry(&h_error_relative,"Statistic_Error_Relative (PhysicsReport Range)","lp");}
    errorleg->AddEntry(&h_Published_Statistic_Error_Relative,"Published_Statistic_Error_Relative","lp");
    errorleg->AddEntry(&h_Published_Systematic_Error_Relative,"Published_Systematic_Error_Relative","lp");
    gStyle->SetLegendTextSize(0.03);
    errorleg->Draw();

    ccerror.Update();
    if (ParametrilizedMode == "No"){
        ccerror.SaveAs((std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/error_compare_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
    }
    else if (ParametrilizedMode == "Yes"){
        ccerror.SaveAs((std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/error_compare_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str());
    }
}


void Plot_PbarNumbers(TGraph *g_antiproton_number, std::string binmerge, std::string issversion, std::string ParametrilizedMode, std::string TRDefficiency, std::string TOFefficiency){

    TCanvas ccPbarNumber("ccPbarNumber","ccPbarNumber",1000,500);

    g_antiproton_number->Draw("");

    g_antiproton_number->SetMarkerStyle(15);
    g_antiproton_number->SetMarkerColor(2);

    gPad->SetLogy();

    g_antiproton_number->SetTitle("");
    g_antiproton_number->GetXaxis()->SetTitle("Rigidity (GV)");
    g_antiproton_number->GetYaxis()->SetTitle("Pbar Numbers");
    g_antiproton_number->GetXaxis()->SetLabelSize(0.04);
    g_antiproton_number->GetXaxis()->SetTitleSize(0.04);
    g_antiproton_number->GetYaxis()->SetLabelSize(0.04);
    g_antiproton_number->GetYaxis()->SetTitleSize(0.04);
    g_antiproton_number->GetYaxis()->SetMoreLogLabels();

    TLegend * PbarNumberleg = new TLegend(0.4,0.75,0.7,0.85); //(xmin, ymin, xmax, ymax)
    PbarNumberleg->SetFillColor(0);
    if (issversion == "pass7.8"){
        PbarNumberleg->AddEntry(g_antiproton_number,"Pbar Numbers (Full Range)","lp");}
    else if (issversion == "2016paper"){
        PbarNumberleg->AddEntry(g_antiproton_number,"Pbar Numbers (2016paper Range)","lp");}
    else if (issversion == "PhysicsReport"){
        PbarNumberleg->AddEntry(g_antiproton_number,"Pbar Numbers (PhysicsReport Range)","lp");}
    gStyle->SetLegendTextSize(0.03);
    PbarNumberleg->Draw();

    ccPbarNumber.Update();
    if (ParametrilizedMode == "No"){
        ccPbarNumber.SaveAs((std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/PbarNumbers_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string(".pdf")).c_str());
    }
    else if (ParametrilizedMode == "Yes"){
        ccPbarNumber.SaveAs((std::string("Time_Averaged_ratio_Low/binmerge") + binmerge + std::string("/plots/PbarNumbers_") + issversion + std::string("_TRDEff_") + TRDefficiency + std::string("_TOFEff_") + TOFefficiency + std::string("_uncertainty.pdf")).c_str());
    }
}





