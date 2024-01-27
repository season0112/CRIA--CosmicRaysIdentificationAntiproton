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
#include "Utilities.hh"
#include "Quantity.hh"
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

#include "TemplateFitterforProtonCC.hh"
#include "TemplateFitter2DforProtonCC.hh"
#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>
#include <stdlib.h>

#include "TemplateFitter2DforProtonCC_CCfixed_dev.hh"
#include "AntiprotonAnalysisTools.hh"

using namespace std;

#define INFO_OUT_TAG "Antiproton_HighTF2D_Unbinned"
#include "debugging.hh"

// Others:
struct TemplateFit2DParameters {
  TemplateFit2DParameters()
    : ccProb(0.0) {
  }

  Utilities::Quantity correctFraction;
  Utilities::Quantity confusedFraction;
  Utilities::Quantity electronFraction;
  double ccProb;
};


struct TemplateFit2DResults {
  TemplateFit2DResults()
    : nEvents(0.0)
    , freeParameters(0)
    , chiSquare(0.0)
    , hFitTotal(nullptr)
    , hFitCorrect(nullptr)
    , hFitConfused(nullptr)
    , hFitElectron(nullptr) {

  }

  ~TemplateFit2DResults() {

    Reset();
  }


  void Reset() {
    nEvents = 0.0;
    freeParameters = 0;
    chiSquare = 0.0;

    delete hFitTotal;
    hFitTotal = nullptr;

    delete hFitCorrect;
    hFitCorrect = nullptr;

    delete hFitConfused;
    hFitConfused = nullptr;

    delete hFitElectron;
    hFitElectron = nullptr;

    fitParameters = TemplateFit2DParameters();
    correct = { 0.0, 0.0 };
    confused = { 0.0, 0.0 };
    electron = { 0.0, 0.0 };
  }

  TemplateFit2DParameters fitParameters;
  Utilities::Quantity correct;
  Utilities::Quantity confused;
  Utilities::Quantity electron;

  double nEvents;
  unsigned int freeParameters;
  double chiSquare;

  TH2* hFitTotal;
  TH2* hFitCorrect;
  TH2* hFitConfused;
  TH2* hFitElectron;

};



// class TemplateFit2D
class TemplateFit2D {
public:
  TemplateFit2D(const TH2* data, const TH2* template_correct, const TH2* template_confused, const TH2* template_electron_Data)
    : fdata(data)
    , ftemplate_correct(template_correct)
    , ftemplate_confused(template_confused)
    , ftemplate_electron_Data(template_electron_Data)
    , nEventsNeg1D(0.0)
    , nEventsPos1D(0.0) {

  }

  bool Minimize(){
    return true;
  }

  bool Perform2DTemplateFit(){
    return Minimize();
  }

private:
  TemplateFit2DParameters GenerateParameters(){
    TemplateFit2DParameters parameters;
    return parameters;
}

  ROOT::Math::Minimizer* SetupMinimizer(){
    static ROOT::Math::Minimizer* sMinimizer = nullptr;
    return sMinimizer;
  }

private:
  const TH2* fdata;
  const TH2* ftemplate_correct;
  const TH2* ftemplate_confused;
  const TH2* ftemplate_electron_Data;
  double nEventsNeg1D;
  double nEventsPos1D;
};



/*
// PerformTemplateFits function
int PerformTemplateFits(const char* resultDirectory, const char* dataAndTemplateInputPath, const char* mcChargeConfusionInputPath, const char* inputFileSuffix, int useEnergyBin) {
int PerformTemplateFits(const TH2* data,
                        const TH2* template_correct,
                        const TH2* template_confused,
                        const TH2* template_electron_Data){
}
*/


// Main function
int main(int argc, char* argv[]) {

  Utilities::ConfigHandler& config = Utilities::ConfigHandler::GetGlobalInstance();
  config.ReadCommandLine(argc, argv);

  config.SetProgramHelpText("Antiproton_HighTF2D_Unbinned",
                            "Template fit for antiproton signal determination in High Rigidity Range.");

  config.AddHelpExample("Antiproton_HighTF2D_Unbinned", "--binningversion 450version --Rigidityofdataset data_negative");

  std::string binningversion = "";
  config.GetValue("OPTIONS", "binningversion", binningversion,
                  "The binningversion is");

  std::string issversion = "";
  config.GetValue("OPTIONS", "issversion", issversion,
                  "The issversion is");

  std::string sign_of_data = "";
  config.GetValue("OPTIONS", "Rigidityofdataset", sign_of_data,
                  "The sign_of_data is:(data_negative or data_positive)");

  float trdlow;
  config.GetValue("OPTIONS", "trdlow value", trdlow,
                  "trdlow value");

  float trdhigh;
  config.GetValue("OPTIONS", "trdhigh value", trdhigh,
                  "trdhigh value");

  if (binningversion == "") {
    WARN_OUT << "No binningversion is given! Please give a binningversion." << std::endl;
    return EXIT_FAIL_CONFIG;
  }

  if (sign_of_data == "") {
    WARN_OUT << "No Rigidityofdataset is given! Please give a Rigidityofdataset." << std::endl;
    return EXIT_FAIL_CONFIG;
  }

  std::vector<std::string> se_array{};
  if (binningversion == "450version"){
    string bin_array[31] = {"14.1_15.3","15.3_16.6","16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175","175_211","211_259","259_450"};
     se_array.insert(se_array.end(),bin_array,bin_array+31);
  }
  if (binningversion == "525version"){
    //string bin_array[32] = {"14.1_15.3","15.3_16.6","16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175","175_211","211_259","259_330","330_525"};
    string bin_array[32] = {"14.1_15.3","15.3_16.6","16.6_18","18_19.5","19.5_21.1","21.1_22.8","22.8_24.7","24.7_26.7","26.7_28.8","28.8_31.1","31.1_33.5","33.5_36.1","36.1_38.9","38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93","93_108","108_125","125_147","147_175","175_211","211_250","250_330","330_525"};
     se_array.insert(se_array.end(),bin_array,bin_array+32);
  }

//string cccut_array[16] = {"0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65","0.70","0.75","0.80","0.85","0.90","0.95"};
string cccut_array[11] = {"0.20","0.25","0.30","0.35","0.40","0.45","0.50","0.55","0.60","0.65", "0.70"};
//string cccut_array[1] = {"0.45"};

//string ccnumber_array[2] = {"20","9"};
//string trdnumber_array[2] = {"12","11"};

string ccnumber_array[4] = {"20","20", "20", "9"};
string trdnumber_array[4] = {"20","16", "12", "11"};

//wrapper<string> w(ccnumber_array,trdnumber_array);
//for (string &ccnumber:ccnumber_array && string &trdnumber:trdnumber_array)

//// Clean the previous result 
if (binningversion == "450version"){
if(NULL==opendir("/results_450version"))
   mkdir("./results_450version",0777);
chdir("./results_450version");
char* currentpath = getcwd(NULL, 0); 
DeleteFile(currentpath);
printf("All files deleted, continue...");
chdir("..");
}
if (binningversion == "525version"){
if(NULL==opendir("/results_525version"))
   mkdir("./results_525version",0777);
chdir("./results_525version");
char* currentpath = getcwd(NULL, 0); 
DeleteFile(currentpath);
printf("All files deleted, continue...");
chdir("..");
}

std::string ccnumber;
std::string trdnumber;
for (int p=0; p<4; p++)
  {
  ccnumber = ccnumber_array[p];
  trdnumber = trdnumber_array[p];
for (auto se : se_array)
    {
        std::cout<< "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!the rigidity range for this fit is :" << se <<std::endl;
    for (string &cccut:cccut_array)
        {
        chdir("./rootfiles");
        //double TRDlow = atof(argv[2]);
        //double TRDhigh = atof(argv[3]);
        //double CClow = atof(argv[4]);

        TFile *file = new TFile();
        if (issversion == "pass7.8")
            file = new TFile((string("histo_") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string(".root")).c_str());
        else if (issversion == "published2016")
            file = new TFile((string("histo_") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_May2015") + string(".root")).c_str());
        else if (issversion == "PhyRep2021")
            file = new TFile((string("histo_") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_Nov2017") + string(".root")).c_str());
        TH2D *template_correct = (TH2D*)file->Get("template_correct");
        template_correct->SetTitle("Antiproton");
        TH2D *template_confused = (TH2D*)file->Get("template_confused");
        template_confused->SetTitle("ChargeConfusedProton");
        TH2D *template_electron = (TH2D*)file->Get("template_electron");
        template_electron->SetTitle("Electron");
        TH2D *template_electron_Data = (TH2D*)file->Get("template_electron_Data");
        template_electron_Data->SetTitle("Electron");
        TH2D *data = (TH2D*)file->Get(sign_of_data.c_str());



        TemplateFit2D templateFitter(data, template_correct, template_confused, template_electron_Data);
 
        //int result = PerformTemplateFits();
        //TemplateFit2DFitFunction fitFunction(data, template_correct, template_confused, template_electron_Data);
        //TemplateFit2DParameters fitParameters = GenerateParameters(identifier, chargeConfusionGraph, startValues);

        /*
        MYUtilities::TemplateFitter2D templateFitter(0);
        templateFitter.AddTemplateHistogram(template_correct);
        templateFitter.AddTemplateHistogram(template_confused);
        templateFitter.AddTemplateHistogram(template_electron_Data);
        templateFitter.SetDataHistogram(data);

        templateFitter.SetStartValue(0, 243);
        templateFitter.SetStartValue(1, 2700);
        templateFitter.SetStartValue(2, 1300);
        templateFitter.Fit(1);

         // store those
        vector<double> result, ResultError;
        assert(result.empty());
        assert(ResultError.empty());
        result.assign(templateFitter.GetAbsoluteResult().begin(), templateFitter.GetAbsoluteResult().end());
        ResultError.assign(templateFitter.GetAbsoluteResultError().begin(), templateFitter.GetAbsoluteResultError().end());
        double Chi2 = templateFitter.Chi2(); 
        int NDF = templateFitter.NDF();
        double CHI2dof = Chi2/NDF;
        */ 

 
        chdir("..");

        /*
        if (binningversion == "450version"){ 
        if(NULL==opendir("/results_450version"))
           mkdir("./results_450version",0777);
        chdir("./results_450version");
        }
        if (binningversion == "525version"){
        if(NULL==opendir("/results_525version"))
           mkdir("./results_525version",0777);
        chdir("./results_525version");
        }

        if (issversion == "pass7.8"){
            ofstream fitresults(string("fit_results") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string(".txt"), ios::app);
            if (fitresults.is_open()){
                fitresults << result << endl;
                fitresults.close();}
            ofstream fiterror(string("fit_results_error") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string(".txt"), ios::app);
            if (fiterror.is_open()){
                fiterror << ResultError << endl;
                fiterror.close();}
            ofstream fitchi2(string("chi2") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string(".txt"), ios::app);
            if (fitchi2.is_open()){
                fitchi2 << CHI2dof << endl;
                fitchi2.close();}
        }
        else if (issversion == "published2016"){
            ofstream fitresults(string("fit_results") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_May2015.txt"), ios::app);
            if (fitresults.is_open()){
                fitresults << result << endl;
                fitresults.close();}
            ofstream fiterror(string("fit_results_error") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_May2015.txt"), ios::app);
            if (fiterror.is_open()){
                fiterror << ResultError << endl;
                fiterror.close();}
            ofstream fitchi2(string("chi2") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_May2015.txt"), ios::app);
            if (fitchi2.is_open()){
                fitchi2 << CHI2dof << endl;
                fitchi2.close();}
        }
        else if (issversion == "PhyRep2021"){
            ofstream fitresults(string("fit_results") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_Nov2017.txt"), ios::app);
            if (fitresults.is_open()){
                fitresults << result << endl;
                fitresults.close();}
            ofstream fiterror(string("fit_results_error") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_Nov2017.txt"), ios::app);
            if (fiterror.is_open()){
                fiterror << ResultError << endl;
                fiterror.close();}
            ofstream fitchi2(string("chi2") + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_Nov2017.txt"), ios::app);
            if (fitchi2.is_open()){
                fitchi2 << CHI2dof << endl;
                fitchi2.close();}
        }

        chdir("..");

        //////////////////////// create Fit Result and Projections plots //////////////////////////////////
        if (templateFitter.HasBeenFit()){ 
            chdir("./fit_plots");
            char path_buf[160];
            printf("current working directory: %s\n", getcwd(path_buf, sizeof(path_buf)));
            if (issversion == "pass7.8"){
                std::string namefit = string("Fit_Result") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string(".pdf");
                std::string nameccprojection = string("CCprojection") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_trdlow_pro:") + argv[2] + string(".pdf");
                std::string nameTRD = string("TRDprojection") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_cclow_pro:") + argv[4]  + string(".pdf");
                TCanvas * c1 = templateFitter.CreateResultDrawing("Fit_Result",800,600);
                c1->SaveAs((char*)namefit.c_str());
                TCanvas * c2 = templateFitter.CreateResultDrawingXprojection("CCprojection",800,600, 0.8, 1.2);
                c2->SaveAs((char*)nameccprojection.c_str());
                TCanvas * c3 = templateFitter.CreateResultDrawingYprojection("TRDprojection",800,600, 0.4, 1.0);
                c3->SaveAs((char*)nameTRD.c_str());
                }
            else if (issversion == "published2016"){
                std::string namefit = string("Fit_Result") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_May2015.pdf");
                std::string nameccprojection = string("CCprojection") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_trdlow_pro:") + argv[2] + string("_May2015.pdf");
                std::string nameTRD = string("TRDprojection") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_cclow_pro:") + argv[4]  + string("_May2015.pdf");
                TCanvas * c1 = templateFitter.CreateResultDrawing("Fit_Result",800,600);
                c1->SaveAs((char*)namefit.c_str());
                TCanvas * c2 = templateFitter.CreateResultDrawingXprojection("CCprojection",800,600, 0.8, 1.2);
                c2->SaveAs((char*)nameccprojection.c_str());
                TCanvas * c3 = templateFitter.CreateResultDrawingYprojection("TRDprojection",800,600, 0.6, 1.0);
                c3->SaveAs((char*)nameTRD.c_str());
                }
            else if (issversion == "PhyRep2021"){
                std::string namefit = string("Fit_Result") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_Nov2017.pdf");
                std::string nameccprojection = string("CCprojection") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_trdlow_pro:") + argv[2] + string("_Nov2017.pdf");
                std::string nameTRD = string("TRDprojection") + se + string("cccut") + cccut + string("CCN") + ccnumber + string("TRDN") + trdnumber + string("_cclow_pro:") + argv[4]  + string("_Nov2017.pdf");
                TCanvas * c1 = templateFitter.CreateResultDrawing("Fit_Result",800,600);
                c1->SaveAs((char*)namefit.c_str());
                TCanvas * c2 = templateFitter.CreateResultDrawingXprojection("CCprojection",800,600, 0.8, 1.2);
                c2->SaveAs((char*)nameccprojection.c_str());
                TCanvas * c3 = templateFitter.CreateResultDrawingYprojection("TRDprojection",800,600, 0.6, 1.0);
                c3->SaveAs((char*)nameTRD.c_str());
                }
            chdir("..");
            }
        */
        }
  } 
  }
  return EXIT_SUCCESS;
}
