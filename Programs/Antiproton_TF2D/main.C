#include "ExampleAnalysisTree.hh"
// ACsoft includes
#include "AnalysisEvent.hh"
#include "ConfigHandler.hh"
#include "EventFactory.hh"
#include "FileManager.hh"
#include "Selector.hh"
#include "SelectionParser.hh"
#include "TreeWriter.hh"
#include "Environment.hh"
#include "ObjectManager.hh"
#include "McSpectrumScaler.hh"
#include <iostream>
#include <string>
#include <cassert>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>
#include <vector>
#include <TCanvas.h>
#include <sstream>
#include <string>
#include <unistd.h>
// Antiproton lib
#include "TemplateFitterforProtonCC.hh"
#include "TemplateFitter2DforProtonCC.hh"
#include "AntiprotonAnalysisTools.hh"
#include "AntiprotonBinning.hh"
// ROOT
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
// Others
#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>
#include <stdlib.h>


using namespace std;

#define INFO_OUT_TAG "Antiproton_TF2D"
#include "debugging.hh"


template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


// Derive ISS/MC uncertainties funcitons.
double LoadRelativeResult(string cccut, string ccnumber, string trdnumber, string issname, int RigidityIndex, string binningversion, string pattern, string issversion, string NNsuffix){
    TFile *f_CCfreeResult = TFile::Open((string("FitResult_Pattern_") + pattern + NNsuffix + string("_") + issversion + string("_") + binningversion + string("_CCFree.root")).c_str(), "READ");
    std::vector<Double_t> *relativevalue_all;
    double relativevalue;
    f_CCfreeResult->GetObject((string("CCprotonNumber_Relative_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion") + issname ).c_str(), relativevalue_all);
    relativevalue = relativevalue_all->at(RigidityIndex);
    return relativevalue;
}

double LoadUncertaintyBand(string cccut, int RigidityIndex, string pattern, string issversion, string binningversion, string NNsuffix){
    //TFile *f_CCfreeResult = TFile::Open((string("../FixedCC/UncertaintyBand_Pattern_") + pattern + NNsuffix + string("_") + issversion + string("_") + binningversion + string(".root")).c_str(), "READ");
    TFile *f_CCfreeResult = TFile::Open((string("../FixedCC/UncertaintyBand_Pattern_") + string("0") + string("_VGG16NN") + string("_") + string("pass7.8") + string("_") + binningversion + string(".root")).c_str(), "READ");
    std::vector<Double_t> *UncertaintyBand_RigidityAll;
    double UncertaintyBand;
    f_CCfreeResult->GetObject("ConfidenceInterval", UncertaintyBand_RigidityAll);  //Linear Fit
    //f_CCfreeResult->GetObject("KMPConfidenceInterval_67_P4", UncertaintyBand_RigidityAll);    //KMP FitP4
    UncertaintyBand = UncertaintyBand_RigidityAll->at(RigidityIndex);
    cout<< "UncertaintyBand:" << UncertaintyBand <<endl;
    return UncertaintyBand;
}


// Global Variables
double ccFixedRelative;
double double_FixedCC_ISSMC_Ratio;


//// Main.
int main(int argc, char* argv[]) {

    Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
    config.ReadCommandLine(argc, argv);
    config.SetProgramHelpText("Antiproton_TF2D",
                              "Template fit for antiproton signal determination in High Rigidity Range.");
    config.AddHelpExample("Antiproton_TF2D", "--binningversion 450version --issversion pass7.8 --pattern 0 --fitmethod ccfree --Rigidityofdataset data_negative");

    string binningversion = "";
    config.GetValue("OPTIONS", "binningversion", binningversion,
                    "The binningversion is");

    string issversion = "";
    config.GetValue("OPTIONS", "issversion", issversion,
                    "The issversion is");

    string pattern = "";
    config.GetValue("OPTIONS", "pattern", pattern,
                    "The choosen tracker pattern is:");

    string fitmethod = "";
    config.GetValue("OPTIONS", "FitMethod", fitmethod,
                    "The FitMethod: ccfree or ccfixed");

    string Rigidityofdataset = "";
    config.GetValue("OPTIONS", "Rigidityofdataset", Rigidityofdataset,
                    "The Sign of ISS data to perform a template fit");

    string ifVGGNN = "";
    config.GetValue("OPTIONS", "ifVGGNN", ifVGGNN,
                    "If you use VGGNN estimator, turn ifVGGNN to yes");

    if (binningversion == "" || issversion == "" || pattern == "" || fitmethod == "" || Rigidityofdataset=="") {
      WARN_OUT << "Some arguments are not given! Please check the help example! " << std::endl;
      return EXIT_FAIL_CONFIG;
    }


    //// Parameters
    std::vector<string> CCcut_array{};
    std::vector<string> ccnumber_array{};
    std::vector<string> trdnumber_array{};
    std::vector<string> RigidityBin_array{};
    std::vector<double> bincenter;
    std::vector<double> subbincenter;

    if (binningversion == "450version"){
        string bin_array[31] = {"14.1_15.3","15.3_16.6","16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175","175_211","211_259","259_450"};
        RigidityBin_array.insert(RigidityBin_array.end(), bin_array, bin_array+31);
        bincenter = AntiprotonNewBinning::NewBinning::AntiprotonBinCenter450().Bins();
        subbincenter.assign(bincenter.begin()+27, bincenter.end()); // 14.7-354.5?
    }
    if (binningversion == "525version"){
        string bin_array[32] = {"14.1_15.3","15.3_16.6","16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175","175_211","211_250","250_330","330_525"};
        RigidityBin_array.insert(RigidityBin_array.end(), bin_array, bin_array+32);
        bincenter = AntiprotonNewBinning::NewBinning::AntiprotonBinCenter525_zhili().Bins();
        subbincenter.assign(bincenter.begin()+27, bincenter.end()); // 14.7-427.5?
    }

    string NNsuffix;
    if ( ifVGGNN == "No" ){
        NNsuffix = "";}
    else if ( ifVGGNN == "Yes" ){
        NNsuffix = "_VGG16NN";}    


    //// Old
    //CCcut_array = {"0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95"};
    //CCcut_array = {"0.00", "0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65", "0.70"};
    //CCcut_array = {"0.00", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45"};
    //ccnumber_array = {"20","20", "20", "9"};
    //trdnumber_array = {"20","16", "12", "11"};  

    //// New
    //CCcut_array = {"0.00", "0.20", "0.35", "0.40", "0.65", "0.70", "0.80", "0.90"};
    //ccnumber_array = {"9", "11", "13", "19", "20", "21", "25", "27", "29"};
    //trdnumber_array = {"9", "11", "13", "19", "20", "21", "25", "27", "29"};


    //CCcut_array     = {"0.00", "0.20", "0.35", "0.40", "0.65", "0.70", "0.80", "0.90"};
    //ccnumber_array  = {"30","20", "20", "20", "9"};
    //trdnumber_array = {"20","20", "16", "12", "11"};

    // Official to Thesis
    //CCcut_array     = {"0.65", "0.20"};
    //ccnumber_array  = {"9"   , "20"};
    //trdnumber_array = {"11"  , "20"};

    // CC Plot test
    // Possible option:
    // CCCut: np.array([0, 0.2, 0.35, 0.4, 0.65, 0.7, 0.8, 0.9])
    // Bin number: [(30,20), (20,20), (20,16), (20,12), (9,11)]
    //CCcut_array     = {"0.00", "0.20", "0.40", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95"};
    //ccnumber_array  = {   "9",    "9",    "9",    "9",    "9",    "9",    "9",    "9",    "9",    "9",    "9"};
    //trdnumber_array = {  "11",   "11",   "11",   "11",   "11",   "11",   "11",   "11",   "11",   "11",   "11"};    

    CCcut_array     = {"0.65", "0.20"};
    ccnumber_array  = {"9"  , "20"};
    trdnumber_array = {"11" , "20"};


    //// Set up, Load root file and Create fit result root file
    if (Rigidityofdataset == "data_negative"){
        chdir(( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/data/templatefit/negative/rootfiles/")).c_str() );}
    else if (Rigidityofdataset == "data_positive"){
        chdir(( string(getenv("HPCHIGHENERGYDATADIR")) + string("/ISS_anylsis/data/templatefit/positive/rootfiles/")).c_str() );}
    string issname = "";
    string ccnumber;
    string trdnumber;

    TFile *file_ProtonTemplate_Data; 
    TFile *file_ElectronTemplate_Data;
    //TFile *file_ElectronTemplate_DataP0;
    TFile *file_ISS_positive;
    TFile *file_ISS_negative;
    TFile *file_ChargeConfusedProtomTemplate_MC = new TFile( ( string("Histo_ChargeConfusedProtomTemplate_MC_Pattern_") + pattern + NNsuffix + string(".root") ).c_str());
    TFile *file_ElectronTemplate_MC             = new TFile( ( string("Histo_ElectronTemplate_MC_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
    TFile *file_FitResult;

    if (issversion == "pass7.8"){
        issname = "";
        file_ProtonTemplate_Data     = new TFile(( string("Histo_ProtonTemplate_Datapass7.8_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        file_ElectronTemplate_Data   = new TFile(( string("Histo_ElectronTemplate_Datapass7.8_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        //file_ElectronTemplate_DataP0 = new TFile(( string("Histo_ElectronTemplate_Datapass7.8_Pattern_0") + NNsuffix + string(".root")).c_str());
        file_ISS_positive            = new TFile(( string("Histo_ISS_positivepass7.8_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        file_ISS_negative            = new TFile(( string("Histo_ISS_negativepass7.8_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
    }
    else if (issversion == "published2016"){
        issname = "_May2015";
        file_ProtonTemplate_Data     = new TFile(( string("Histo_ProtonTemplate_Datapublished2016_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        file_ElectronTemplate_Data   = new TFile(( string("Histo_ElectronTemplate_Datapublished2016_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        //file_ElectronTemplate_DataP0 = new TFile(( string("Histo_ElectronTemplate_Datapublished2016_Pattern_0") + NNsuffix + string(".root")).c_str());
        file_ISS_positive            = new TFile(( string("Histo_ISS_positivepublished2016_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        file_ISS_negative            = new TFile(( string("Histo_ISS_negativepublished2016_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
    }
    else if (issversion == "PhyRep2021"){
        issname = "_Nov2017";
        file_ProtonTemplate_Data     = new TFile(( string("Histo_ProtonTemplate_DataPhyRep2021_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        file_ElectronTemplate_Data   = new TFile(( string("Histo_ElectronTemplate_DataPhyRep2021_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        //file_ElectronTemplate_DataP0 = new TFile(( string("Histo_ElectronTemplate_DataPhyRep2021_Pattern_0") + NNsuffix + string(".root")).c_str());
        file_ISS_positive            = new TFile(( string("Histo_ISS_positivePhyRep2021_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
        file_ISS_negative            = new TFile(( string("Histo_ISS_negativePhyRep2021_Pattern_") + pattern + NNsuffix + string(".root")).c_str());
    }

    if (fitmethod == "ccfree"){
         file_FitResult  = new TFile( (string("../FitResult/FitResult_Pattern_") + pattern + NNsuffix + string("_") + issversion + ("_") + binningversion + string("_CCFree.root")).c_str(), "RECREATE");
    }
    else if (fitmethod == "ccfixed"){
         file_FitResult = new TFile( (string("../FitResult/FitResult_Pattern_") + pattern + NNsuffix + string("_") + issversion + ("_") + binningversion + string("_CCFixed.root")).c_str(), "RECREATE");
    }
    chdir("..");


    //// Loop
    for (int p=0; p<ccnumber_array.size(); p++){
        ccnumber  = ccnumber_array.at(p);
        trdnumber = trdnumber_array.at(p);

    /*
    for (int p=0; p<ccnumber_array.size(); p++){
        ccnumber = ccnumber_array.at(p);

        for (int k=0; k<trdnumber_array.size(); k++){
           trdnumber = trdnumber_array.at(k); 
    */

            for (auto cccut : CCcut_array){
                int RigidityIndex = 0;
                std::vector<double> chi2dof_array;
                std::vector<double> AntiprotonNumber_array;
                std::vector<double> CCprotonNumber_array;
                std::vector<double> ElectronNumber_array; 
                std::vector<double> AntiprotonNumberError_array;
                std::vector<double> CCprotonNumberError_array;
                std::vector<double> ElectronNumberError_array;

                std::vector<double> AntiprotonNumber_Relative_array;
                std::vector<double> CCprotonNumber_Relative_array;
                std::vector<double> ElectronNumber_Relative_array;
                std::vector<double> AntiprotonNumber_RelativeError_array;
                std::vector<double> CCprotonNumber_RelativeError_array;
                std::vector<double> ElectronNumber_RelativeError_array;

                std::vector<double> AntiprotonGlobalCorrelationFactors_array;
                std::vector<double> CCprotonGlobalCorrelationFactors_array;
                std::vector<double> ElectronGlobalCorrelationFactors_array;

                std::vector<double> ErrorMatrix_array;


                for (auto RigidityBin : RigidityBin_array){

                    cout<< "\n" <<endl;
                    cout<< "********************************************************************" <<endl;
                    cout<< "Rigidity range for this fit is :" << RigidityBin <<endl;
                    cout<< "cccut is:" + cccut + '\n' + "ccnumber is:" + ccnumber + '\n' + "trdnumber is:" + trdnumber <<endl;
                   
                    // Load Fixed Relative Result (To calculate CC error)
                    if (fitmethod == "ccfixed"){
                        chdir("./FitResult");
                        ccFixedRelative = LoadRelativeResult(cccut, ccnumber, trdnumber, issname, RigidityIndex, binningversion, pattern, issversion, NNsuffix);
                        double_FixedCC_ISSMC_Ratio = LoadUncertaintyBand(cccut, RigidityIndex, pattern, issversion, binningversion, NNsuffix);
                        chdir("..");
                    }

                    // Load Template Histograms
                    TH2D *template_correct = (TH2D*)file_ProtonTemplate_Data->Get( (string("ProtonTemplate_Data_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber).c_str() );
                    template_correct->SetTitle("Antiproton");

                    TH2D *template_confused = (TH2D*)file_ChargeConfusedProtomTemplate_MC->Get( (string("ChargeConfusedProtomTemplate_MC_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber).c_str() );
                    template_confused->SetTitle("CCProton");
                    template_confused->SetName("CCProton"); 

                    TH2D *template_electron_MC = (TH2D*)file_ElectronTemplate_MC->Get( (string("ElectronTemplate_MC_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber).c_str() );
                    template_electron_MC->SetTitle("Electron");

                    TH2D *template_electron_Data = (TH2D*)file_ElectronTemplate_Data->Get( (string("ElectronTemplate_Data_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber).c_str() );
                    template_electron_Data->SetTitle("Electron");

                    //TH2D *template_electron_Data_P0 = (TH2D*)file_ElectronTemplate_DataP0->Get( (string("ElectronTemplate_Data_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber).c_str() );
                    //template_electron_Data_P0->SetTitle("Electron");

                    // Nomarlizing template histograms
                    Double_t factor = 1;
                    template_correct       ->Scale(factor/template_correct->Integral());
                    template_confused      ->Scale(factor/template_confused->Integral());
                    template_electron_MC   ->Scale(factor/template_electron_MC->Integral());
                    template_electron_Data ->Scale(factor/template_electron_Data->Integral());

                    // Load Data Histograms
                    TH2D *ISSNegativedata = (TH2D*)file_ISS_negative->Get((string("ISS_negative_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_" + trdnumber)).c_str());

                    TH2D *ISSPositivedata = (TH2D*)file_ISS_positive->Get((string("ISS_positive_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_" + trdnumber)).c_str());

                    MYUtilities::TemplateFitter2D templateFitter(-1); // printlevel: output: 0 standard, -1 quiet
                    templateFitter.AddTemplateHistogram(template_correct);
                    templateFitter.AddTemplateHistogram(template_confused);
                    //templateFitter.AddTemplateHistogram(template_electron_Data);
                    //templateFitter.AddTemplateHistogram(template_electron_Data_P0);
                    templateFitter.AddTemplateHistogram(template_electron_MC);
                   
                    if (Rigidityofdataset == "data_negative"){
                        templateFitter.SetDataHistogram(ISSNegativedata);}
                    else if (Rigidityofdataset == "data_positive"){
                        templateFitter.SetDataHistogram(ISSPositivedata);}

                    // Check numbers
                    cout<< "template_correct: "     << template_correct->GetEntries()     << endl;
                    cout<< "template_confused: "    << template_confused->GetEntries()    << endl;
                    cout<< "template_electron_MC: " << template_electron_MC->GetEntries() << endl;
                    cout<< "ISSNegativedata: "      << ISSNegativedata->GetEntries()      << endl;


                    // Perform Fit
                    if (fitmethod == "ccfixed"){
                        cout << "ccFixedRelative: "     << ccFixedRelative            << endl;
                        cout << "FixedCC_ISSMC_Ratio: " << double_FixedCC_ISSMC_Ratio << endl;
                        cout << "ccFixedRelative with uncertainty:" +  to_string_with_precision( ccFixedRelative * (1 - double_FixedCC_ISSMC_Ratio) ) << endl;
                        templateFitter.FixParameter(1, ccFixedRelative * (1 - double_FixedCC_ISSMC_Ratio) ); // Fixed CCProton numbers for fit. (templateindex, fixed value): fixed value is relative and must be between 0 and 1.   +/- for fit?
                    }

                    templateFitter.Fit(1);  // 0:Ch2 fit, 1:Likelihood Fit


                    // Template fit result
                    vector<double> result, ResultError, RelativeResult, RelativeResultError, GlobalCorrelationFactors, ErrorMatrix;

                    assert(result             .empty());
                    assert(ResultError        .empty());
                    assert(RelativeResult     .empty());
                    assert(RelativeResultError.empty());
                    assert(RelativeResult     .empty());

                    result.assign(templateFitter.GetAbsoluteResult().begin()          , templateFitter.GetAbsoluteResult().end());       //MIGRAD minimization for result, Hesse:compute asymptotic errors. 
                    ResultError.assign(templateFitter.GetAbsoluteResultError().begin(), templateFitter.GetAbsoluteResultError().end());  // symmetric MIGRAD uncertainty
                    RelativeResult.assign(templateFitter.GetRelativeResult().begin()  , templateFitter.GetRelativeResult().end());
                    RelativeResultError.assign(templateFitter.GetRelativeResultError().begin()          , templateFitter.GetRelativeResultError().end()); // symmetric MIGRAD uncertainty
                    GlobalCorrelationFactors.assign(templateFitter.GetGlobalCorrelationFactors().begin(), templateFitter.GetGlobalCorrelationFactors().end()); // global correlation factors as determined by TMinuit from the error matrix
                    //ErrorMatrix.assign(templateFitter.GetErrorMatrix().begin()                          , templateFitter.GetErrorMatrix().end());

                    double Chi2    = templateFitter.Chi2(); 
                    int NDF        = templateFitter.NDF();
                    double CHI2dof = Chi2/NDF;

                    chi2dof_array.insert(chi2dof_array.end(), CHI2dof);
                    AntiprotonNumber_array.insert(AntiprotonNumber_array.end()          , result.at(0));
                    CCprotonNumber_array.insert(CCprotonNumber_array.end()              , result.at(1));
                    ElectronNumber_array.insert(ElectronNumber_array.end()              , result.at(2));
                    AntiprotonNumberError_array.insert(AntiprotonNumberError_array.end(), ResultError.at(0));
                    CCprotonNumberError_array.insert(CCprotonNumberError_array.end()    , ResultError.at(1));
                    ElectronNumberError_array.insert(ElectronNumberError_array.end()    , ResultError.at(2));

                    AntiprotonNumber_Relative_array.insert(AntiprotonNumber_Relative_array.end(), RelativeResult.at(0));
                    CCprotonNumber_Relative_array.insert(CCprotonNumber_Relative_array.end()    , RelativeResult.at(1));
                    ElectronNumber_Relative_array.insert(ElectronNumber_Relative_array.end()    , RelativeResult.at(2));
                    AntiprotonNumber_RelativeError_array.insert(AntiprotonNumber_RelativeError_array.end(), RelativeResultError.at(0));
                    CCprotonNumber_RelativeError_array.insert(CCprotonNumber_RelativeError_array.end()    , RelativeResultError.at(1));
                    ElectronNumber_RelativeError_array.insert(ElectronNumber_RelativeError_array.end()    , RelativeResultError.at(2));

                    AntiprotonGlobalCorrelationFactors_array.insert(AntiprotonGlobalCorrelationFactors_array.end(), GlobalCorrelationFactors.at(0));
                    CCprotonGlobalCorrelationFactors_array.insert(CCprotonGlobalCorrelationFactors_array.end()    , GlobalCorrelationFactors.at(1));
                    ElectronGlobalCorrelationFactors_array.insert(ElectronGlobalCorrelationFactors_array.end()    , GlobalCorrelationFactors.at(2));

                    ErrorMatrix_array.insert(ErrorMatrix_array.end(), templateFitter.GetErrorMatrix().begin()                          , templateFitter.GetErrorMatrix().end());
                    //ErrorMatrix_array.push_back("\n"); // needed ??? ErrorMatrix_array is vector<double> 


                    // create Fit Result and Projections plots (Only for CC free fit)
                    if (fitmethod == "ccfree"){
                        if (templateFitter.HasBeenFit()){ 
                            chdir("./Fit_plots_CCfree");
                            //char path_buf[160];
                            //printf("current working directory: %s\n", getcwd(path_buf, sizeof(path_buf)));

                            string namefit          = string("FitResult_Pattern_")     + pattern + NNsuffix + string("_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + issname + string(".pdf");
                            string nameccprojection = string("CCprojection_Pattern_")  + pattern + NNsuffix + string("_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + issname + string(".pdf");
                            string nameTRD          = string("TRDprojection_Pattern_") + pattern + NNsuffix + string("_") + RigidityBin + string("_cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + issname + string(".pdf");

                            int BinRemovedNumber_FromLeft = 0;
                            int RebinNumber_X = 1;
                            int BinRemovedNumber_FromRight = 0;
                            int RebinNumber_Y = 1;
                            TCanvas * c1 = templateFitter.CreateResultDrawing           ("Fit_Result",800,600);
                            c1->SaveAs((char*)namefit.c_str());
                            
                            TCanvas * c2 = templateFitter.CreateResultDrawingXprojection("CCprojection" , 800, 600, 0.7, 1.6, BinRemovedNumber_FromLeft , RebinNumber_X); //0.8, 1.2
                            c2->SaveAs((char*)nameccprojection.c_str());
                            TCanvas * c3 = templateFitter.CreateResultDrawingYprojection("TRDprojection", 800, 600, 0.75, 1.0, BinRemovedNumber_FromRight, RebinNumber_Y); //0.4, 1.0
                            c3->SaveAs((char*)nameTRD.c_str());
                            
                            //TCanvas * c2 = templateFitter.CreateResultDrawingXprojection("CCprojection" , 800, 600, 0.6, 1.6); //0.8, 1.2
                            //c2->SaveAs((char*)nameccprojection.c_str());
                            //TCanvas * c3 = templateFitter.CreateResultDrawingYprojection("TRDprojection", 800, 600, 0.6, 1.0); //0.4, 1.0
                            //c3->SaveAs((char*)nameTRD.c_str());

                            chdir("..");
                        }
                    }

                    templateFitter.Clear();

                    RigidityIndex = RigidityIndex +1;
                } // End of rigidity loop


                // Save Template fit Result
                chdir("./FitResult");
                TGraph *g_CHI2dof = new TGraph(RigidityBin_array.size(), subbincenter.data(), chi2dof_array.data());

                file_FitResult->cd();
                g_CHI2dof->Write();

                file_FitResult->WriteObject(&chi2dof_array, (string("chi2dof_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );

                file_FitResult->WriteObject(&AntiprotonNumber_array, (string("AntiprotonNumber_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&CCprotonNumber_array  , (string("CCprotonNumber_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&ElectronNumber_array  , (string("ElectronNumber_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );

                file_FitResult->WriteObject(&AntiprotonNumberError_array, (string("AntiprotonNumberError_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&CCprotonNumberError_array  , (string("CCprotonNumberError_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&ElectronNumberError_array  , (string("ElectronNumberError_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );

                file_FitResult->WriteObject(&AntiprotonNumber_Relative_array, (string("AntiprotonNumber_Relative_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&CCprotonNumber_Relative_array  , (string("CCprotonNumber_Relative_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&ElectronNumber_Relative_array  , (string("ElectronNumber_Relative_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );

                file_FitResult->WriteObject(&AntiprotonNumber_RelativeError_array, (string("AntiprotonNumber_RelativeError_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&CCprotonNumber_RelativeError_array  , (string("CCprotonNumber_RelativeError_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&ElectronNumber_RelativeError_array  , (string("ElectronNumber_RelativeError_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );

                file_FitResult->WriteObject(&AntiprotonGlobalCorrelationFactors_array, (string("AntiprotonGlobalCorrelationFactors_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&CCprotonGlobalCorrelationFactors_array  , (string("CCprotonGlobalCorrelationFactors_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );
                file_FitResult->WriteObject(&ElectronGlobalCorrelationFactors_array  , (string("ElectronGlobalCorrelationFactors_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );

                file_FitResult->WriteObject(&ErrorMatrix_array, (string("ErrorMatrix_") + string("cccut_") + cccut + string("_CCN_") + ccnumber + string("_TRDN_") + trdnumber + string("_ISSVersion")+ issname).c_str() );

                chdir("..");
        /*
            } // End of cccut loop
        }  //End of trd number loop
        */
        }  //End of cccut and trd number loop
    } // End of cc number loop

    return EXIT_SUCCESS;
}
