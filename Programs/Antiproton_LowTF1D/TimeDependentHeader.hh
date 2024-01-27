#include "AntiprotonAnalysisTools.hh"

std::tuple<std::vector<std::string>> LoadTimeStampBins(std::string timeindex){
    std::vector<std::string> v_time;
    std::string st;
    std::ifstream in2(timeindex);
    while (getline(in2,st))
        v_time.push_back(st);
    in2.close();

    return {v_time};
}


void CreateTxtFiles(std::string timemode, std::vector<double> binning, int index, std::string binmerge, std::string bartals){
    chdir("./results");
    if (timemode == "3BartalRotation" || timemode == "6BartalRotation"){
        std::ofstream fitresults(std::string("fit_results_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::trunc);
        std::ofstream fitresults_error(std::string("fit_results_error_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::trunc);
        std::ofstream time_stamp(std::string("timestamp_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::trunc);
        std::ofstream fit_chi2dof(std::string("fit_chi2dof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::trunc);
        std::ofstream antiprotonnumber(std::string("antiprotonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::trunc);
        std::ofstream protonnumber(std::string("protonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::trunc);
    }
    else if(timemode == "6months"){
        std::ofstream fitresults(std::string("fit_results_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream fitresults_error(std::string("fit_results_error_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream time_stamp(std::string("timestamp_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream fit_chi2dof(std::string("fit_chi2dof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream antiprotonnumber(std::string("antiprotonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
        std::ofstream protonnumber(std::string("protonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::trunc);
    }
    chdir("..");
}


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


std::tuple<TH2F *> LoadAntiprotonMCTemplate_TimeDependent(std::string binmerge, std::vector<double> binning, int index, TCut AntiprotonCut, int TrdBINNUMBER, double TrdLOW, double TrdHIGH, int RECITOFBETANUMBER, double RECITOFBETALOW, double RECITOFBETAHIGH){
    TChain *fantiproton = new TChain("AntiprotonLowEnergyTree");
    for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
        fantiproton->AddFile((std::string("B1042_antipr.pl1.1800_7.6_all_Tree_") + doubleToString(binning[index+binindex]) + std::string("_") + doubleToString(binning[index+binindex+1]) + std::string(".root")).c_str());    }
    fantiproton->Draw((std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_antiproton(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(",") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), AntiprotonCut);
    TH2F *Tof_template_antiproton = (TH2F*)gDirectory->Get("th2f_antiproton");
    Tof_template_antiproton->SetTitle("Antiproton");
    delete fantiproton;

    return {Tof_template_antiproton};
}


std::tuple<TH2F *> LoadElectronMCTemplate_TimeDependent(std::string binmerge, std::vector<double> binning, int index, TCut ElectronCut, int TrdBINNUMBER, double TrdLOW, double TrdHIGH, int RECITOFBETANUMBER, double RECITOFBETALOW, double RECITOFBETAHIGH){

    TChain *felectron = new TChain("AntiprotonLowEnergyTree");
    for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
        felectron->AddFile((std::string("B1091_el.pl1.0_25200_7.6_all_Tree_negative_") + doubleToString(binning[index+binindex]) + std::string("_") + doubleToString(binning[index+binindex+1]) + std::string(".root")).c_str());
    }
    felectron->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_electron(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(),  ElectronCut);
    TH2F *TofTRD_template_electron = (TH2F*)gDirectory->Get("th2f_electron");
    TofTRD_template_electron->SetTitle("Electron");
    delete felectron;

    return {TofTRD_template_electron};
}


std::tuple<TH2F *> LoadPionMCTemplate_TimeDependent(std::string binmerge, std::vector<double> binning, int index, TCut PionCut, int TrdBINNUMBER, double TrdLOW, double TrdHIGH, int RECITOFBETANUMBER, double RECITOFBETALOW, double RECITOFBETAHIGH){

    TChain *fpion = new TChain("AntiprotonLowEnergyTree");
    for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
        fpion->AddFile((std::string("B1220_pr.pl1phpsa.0550.4_00_7.8_all_Tree_negative_") + doubleToString(binning[index+binindex]) + std::string("_") + doubleToString(binning[index+binindex+1]) + std::string(".root")).c_str());}
    fpion->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_pion(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), PionCut);
    TH2F *TofTRD_template_pion = (TH2F*)gDirectory->Get("th2f_pion");
    TofTRD_template_pion->SetTitle("Pion");
    delete fpion;

    return {TofTRD_template_pion};
}


std::tuple<TH2F *, TH2F *, TH2F *> LoadISSNegativeData_TimeDependent(std::string binmerge, std::vector<double> binning, int index, TCut NegativeCut, TCut ElectronTemplateDataCut, TCut PionTemplateDataCut,  int TrdBINNUMBER, double TrdLOW, double TrdHIGH, int RECITOFBETANUMBER, double RECITOFBETALOW, double RECITOFBETAHIGH, int i, std::string bartals, int maxindex, std::string timemode){

        TChain *fpass7_negative = new TChain("tree");
        for (int treeindex=i; treeindex < i+stoi(bartals); treeindex++){
            /*
            if (treeindex==maxindex){
                std::cout<< "treeindex="           << treeindex << std::endl;
                std::cout<< "maxindex="            << maxindex << std::endl;
                std::cout<< "break loop negative!" << std::endl;
                break;
            }
            */
            if (timemode == "3BartalRotation" || timemode == "6BartalRotation"){
                for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                    std::string left = doubleToString(binning[index+binindex]);
                    std::string right = doubleToString(binning[index+binindex+1]);
                    if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                        left = to_string_with_precision(binning[index+binindex], 1);}
                    if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                        right = to_string_with_precision(binning[index+binindex+1], 1);}
                    fpass7_negative->AddFile((std::string("rootfiles/negative/negative_") + left + std::string("_") + right + std::string("_") + std::to_string(treeindex) + std::string(".root")).c_str());
                }
            }
            else if(timemode == "6months"){
                for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                    std::string left = doubleToString(binning[index+binindex]);
                    std::string right = doubleToString(binning[index+binindex+1]);
                    if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                        left = to_string_with_precision(binning[index+binindex], 1);}
                    if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                        right = to_string_with_precision(binning[index+binindex+1], 1);}
                    fpass7_negative->AddFile((std::string("rootfiles_6months/negative/negative_") + left + std::string("_") + right + std::string("_") + std::to_string(treeindex) + std::string(".root")).c_str());
                }
            }
        }

        fpass7_negative->Draw((std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_negative(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), NegativeCut);
        TH2F *TofTRD_data_pass7_negative = (TH2F*)gDirectory->Get("th2f_negative");
        TofTRD_data_pass7_negative->SetTitle("ISSnegative");

        fpass7_negative->Draw((std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_ElectronDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), ElectronTemplateDataCut);
        TH2F *TofTRD_template_ElectronData = (TH2F*)gDirectory->Get("th2f_ElectronDataTemplate");
        TofTRD_template_ElectronData->SetTitle("Electron");

        fpass7_negative->Draw((std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_PionDataTemplate(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(", ") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), PionTemplateDataCut);
        TH2F *TofTRD_template_PionData = (TH2F*)gDirectory->Get("th2f_PionDataTemplate");
        TofTRD_template_PionData->SetTitle("Pion");

    return {TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData};
}


std::tuple<TH2F *, TH2F *> LoadISSPositiveData_TimeDependent(std::string binmerge, std::vector<double> binning, int index, TCut AntiprotonCut, TCut PositiveCut, int TrdBINNUMBER, double TrdLOW, double TrdHIGH, int RECITOFBETANUMBER, double RECITOFBETALOW, double RECITOFBETAHIGH, int i, std::string bartals, int maxindex, std::string timemode){

        TChain *fpass7_positive = new TChain("tree");
        for (int treeindex=i; treeindex < i+stoi(bartals); treeindex++){
            /*
            if (treeindex==maxindex){
                std::cout<< "treeindex=" << treeindex << std::endl;
                std::cout<< "maxindex=" << maxindex   << std::endl;
                std::cout<< "break loop negative!"    << std::endl;
                break;
            }
            */
            if (timemode == "3BartalRotation" || timemode == "6BartalRotation"){
                for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                    std::string left = doubleToString(binning[index+binindex]);
                    std::string right = doubleToString(binning[index+binindex+1]);
                    if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                        left = to_string_with_precision(binning[index+binindex], 1);}
                    if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                        right = to_string_with_precision(binning[index+binindex+1], 1);}
                    fpass7_positive->AddFile((std::string("rootfiles/positive/positive_") + left + std::string("_") + right + std::string("_") + std::to_string(treeindex) + std::string(".root")).c_str());
                }
            }
            else if(timemode == "6months"){
                for (int binindex=0; binindex<stoi(binmerge); binindex=binindex+1){
                    std::string left = doubleToString(binning[index+binindex]);
                    std::string right = doubleToString(binning[index+binindex+1]);
                    if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
                        left = to_string_with_precision(binning[index+binindex], 1);}
                    if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
                        right = to_string_with_precision(binning[index+binindex+1], 1);}
                    fpass7_positive->AddFile((std::string("rootfiles_6months/positive/positive_") + left + std::string("_") + right + std::string("_") + std::to_string(treeindex) + std::string(".root")).c_str());
                }
            }
        }

        //fpass7_positive->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_positive_template(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(",") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), AntiprotonCut);
        fpass7_positive->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_positive_template(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(",") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), PositiveCut);
        TH2F *TofTRD_data_pass7_positive_template = (TH2F*)gDirectory->Get("th2f_positive_template");
        TofTRD_data_pass7_positive_template->SetTitle("ISSpositive_template");

        fpass7_positive->Draw( (std::string("1./TofBeta - 1/ ( sqrt (pow(Rigidity,2) / ( pow(0.938,2) + pow(Rigidity,2) )) ):TrdLogLikelihoodRatioElectronProtonTracker>>th2f_positive_data(") + std::to_string(TrdBINNUMBER) + std::string(", ") + std::to_string(TrdLOW) + std::string(",") + std::to_string(TrdHIGH) + std::string(", ") + std::to_string(RECITOFBETANUMBER) + std::string(", ") + std::to_string(RECITOFBETALOW) + std::string(", ") + std::to_string(RECITOFBETAHIGH) + std::string(")")).c_str(), PositiveCut);
        TH2F *TofTRD_data_pass7_positive_data = (TH2F*)gDirectory->Get("th2f_positive_data");
        TofTRD_data_pass7_positive_data->SetTitle("ISSpositive_data");

    return {TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data};
}



std::tuple<TH2F *, TH2F *, TH2F *, TH2F *, TH2F *>LoadTemplates_TimeDependent(std::string binmerge, std::vector<double> binning, int index, int TimeIndex, std::string timemode, std::string ParametrilizedMode, int GeneraterNumberIndex, std::string SysErrEffStudyMode, std::string trdeff, std::string tofeff){

    std::string left = doubleToString(binning[index]);
    std::string right = doubleToString(binning[index+stoi(binmerge)]);
    if ( left == doubleToString(1) || left == doubleToString(11) || left == doubleToString(12) || left == doubleToString(13) || left == doubleToString(18)){
        left = to_string_with_precision(binning[index], 1);}
    if ( right == doubleToString(1) || right == doubleToString(11) || right == doubleToString(12) || right == doubleToString(13) || right == doubleToString(18) ){
        right = to_string_with_precision(binning[index+stoi(binmerge)], 1);}

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

    TFile *f_histogram = new TFile( ( filename + std::string("/TimeDependentTemplatesAndData_") + left + std::string("_") + right + std::string("_") + std::to_string(TimeIndex) + nameindex + std::string(".root")).c_str(), "READ");

    TH2F *TofTRD_template_ElectronData        = (TH2F*) f_histogram->Get( ((std::string("TofTRD_template_ElectronData")        + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + templateindex).c_str() )); 
    TH2F *TofTRD_data_pass7_positive_template = (TH2F*) f_histogram->Get( ((std::string("TofTRD_data_pass7_positive_template") + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + templateindex).c_str() ));
    TH2F *TofTRD_template_PionData            = (TH2F*) f_histogram->Get( ((std::string("TofTRD_template_PionData")            + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff + templateindex).c_str() ));
    TH2F *TofTRD_data_pass7_positive_data     = (TH2F*) f_histogram->Get( ((std::string("TofTRD_data_pass7_positive_data")     + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff                ).c_str() ));
    TH2F *TofTRD_data_pass7_negative          = (TH2F*) f_histogram->Get( ((std::string("TofTRD_data_pass7_negative")          + std::string("_TRDeff_") + trdeff + std::string("_TOFeff_") + tofeff                ).c_str() ));

    TofTRD_template_ElectronData       ->SetDirectory(0);
    TofTRD_data_pass7_positive_template->SetDirectory(0);
    TofTRD_template_PionData           ->SetDirectory(0);
    TofTRD_data_pass7_positive_data    ->SetDirectory(0);
    TofTRD_data_pass7_negative         ->SetDirectory(0); 

    f_histogram->Close();

    return {TofTRD_data_pass7_negative, TofTRD_template_ElectronData, TofTRD_template_PionData, TofTRD_data_pass7_positive_template, TofTRD_data_pass7_positive_data};
}





void SaveResult(std::string bartals, std::vector<double> binning, int index, std::string binmerge, std::vector<double> result_tof, double ProtonNumber, std::vector<double> ResultError_tof, double CHI2dof_tof){
        if (stoi(bartals) == 3 || stoi(bartals) == 6){
            chdir("./results");
            std::ofstream fitresults(std::string("fit_results_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
            if (fitresults.is_open()){
                fitresults << result_tof[0]/ProtonNumber << std::endl;
                fitresults.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream fitresults_error(std::string("fit_results_error_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
            if (fitresults_error.is_open()){
                fitresults_error << ResultError_tof[0]/ProtonNumber << std::endl;
                fitresults_error.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream fit_chi2dof(std::string("fit_chi2dof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
            if (fit_chi2dof.is_open()){
                fit_chi2dof << CHI2dof_tof <<std::endl;
                fit_chi2dof.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream antiprotonnumber(std::string("antiprotonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
            if (antiprotonnumber.is_open()){
                antiprotonnumber<< result_tof[0] << std::endl;
                antiprotonnumber.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream protonnumber(std::string("protonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + bartals + std::string("BartalRotation.txt"), std::ios::app);
            if (protonnumber.is_open()){
                protonnumber<< ProtonNumber <<std::endl;
                protonnumber.close();
            }
            chdir("..");
        }
        else if (stoi(bartals) == 1){
            chdir("./results");
            std::ofstream fitresults(std::string("fit_results_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
            if (fitresults.is_open()){
                fitresults << result_tof[0]/ProtonNumber << std::endl;
                fitresults.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream fitresults_error(std::string("fit_results_error_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
            if (fitresults_error.is_open()){
                fitresults_error << ResultError_tof[0]/ProtonNumber << std::endl;
                fitresults_error.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream fit_chi2dof(std::string("fit_chi2dof_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
            if (fit_chi2dof.is_open()){
                fit_chi2dof << CHI2dof_tof <<std::endl;
                fit_chi2dof.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream antiprotonnumber(std::string("antiprotonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
            if (antiprotonnumber.is_open()){
                antiprotonnumber<< result_tof[0] << std::endl;
                antiprotonnumber.close();
            }
            chdir("..");

            chdir("./results");
            std::ofstream protonnumber(std::string("protonnumber_") + doubleToString(binning[index]) + std::string("_") + doubleToString(binning[index+stoi(binmerge)]) + std::string("_") + std::string("6months.txt"), std::ios::app);
            if (protonnumber.is_open()){
                protonnumber<< ProtonNumber <<std::endl;
                protonnumber.close();
            }
            chdir("..");
        }
}


//// Cuts definations

TCut ExtrapolatedPhotoElectronsFromTracker ="NPhotoElectrons-ExtrapolatedRichExpectedPhotoElectronsProton<20 && NPhotoElectrons>-999 && ExtrapolatedRichTileIndex == TileIndex || NPhotoElectrons==-999 || ExtrapolatedRichTileIndex != TileIndex";
TCut ProtonCCMVABDCut = "ProtonCCMVABDT > 0.9";

//Trigger
TCut TriggerPhysics = "TriggerFlags != 1 && TriggerFlags != 64 ";
//TRD
//TCut TrdLikelihoodCut = "TrdLogLikelihoodRatioElectronProtonTracker > 0.7 || TrdLogLikelihoodRatioElectronProtonTracker==-1.5";    // Antiproton:1to1.2, Electron:0.5
TCut TrdLikelihoodHeProtonCut = "TrdLogLikelihoodRatioProtonHeliumTracker < 0.1";
//ECAL
TCut EcalBDT_EnergyDCut = "EcalBDT_EnergyD < -0.9";
//RICH
TCut RichBetaCut = "RichBeta==0 || RichIsNaF==1";
//TCut RichBetaCut = "RichIsNaF==1";
TCut RichNaFCut = "RichIsNaF==1";
TCut RichAglCut = "RichIsNaF==0";
TCut RichPhotoElectron = "NPhotoElectrons-NExpectedPhotoElectrons<10  || NPhotoElectrons==-999";
TCut RichChargeCut = "RichCharge<2 || RichCharge==0";
//TOF
//TCut BetaConverted = "BetaConverted<0.9";
TCut TofMassonecharge = "TofMassonecharge>0.5";
//TCut TOFBETALikelihood = "TOFBETALikelihood<0.7";
// Tracker

// New cuts added:
TCut TrdSegmentsXZNumberCut = "TrdSegmentsXZNumber==1";
TCut TrdSegmentsYZNumberCut = "TrdSegmentsYZNumber==1";
TCut TRDVTracksSizeCut = "TRDVTracksSize==1";
TCut TrdNumberOfHitsCut = "TrdNumberOfHits<50";
TCut ACCHitsCut = "ACCHits==0";













