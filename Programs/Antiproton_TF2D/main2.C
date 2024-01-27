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

#include "TemplateFitterforProtonCC.hh"
#include "TemplateFitter2DforProtonCC.hh"
#include "dirent.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/io.h>

using namespace std;

#define INFO_OUT_TAG "Antiproton_TF2D_numberfree"
#include "debugging.hh"

int main(int argc, char* argv[]) {
(void) argc;

string se_array[19] = {"38.9_41.9","41.9_45.1","45.1_48.5","48.5_52.2","52.2_56.1","56.1_60.3","60.3_64.8","64.8_69.7","69.7_74.9","74.9_80.5","80.5_93.0","93.0_108.0","108.0_125.0","125.0_147.0","147_175","175_211","211_259","259_330","330_525"};
string ccnumber_array[5] = {"9","11","13","15","17"};
string trdnumber_array[4] = {"10","11","12","13"};


for (string &se:se_array)
  {
for (string &ccnumber:ccnumber_array)
    {
for (string &trdnumber:trdnumber_array)
      {
printf("current working directory: %s\n", getcwd(nullptr, 0));
chdir("./rootfiles");
printf("current working directory: %s\n", getcwd(nullptr, 0));
std::string s1 = "histo_";
std::string s2 = ".root";
std::string fitresultname = "fit_results";
std::string fiterrorname = "fit_results_error";
std::string chi2name = "chi2";
std::string cccut = argv[5];////////////////////////////////////////////////// index check
//std::string ccnumber = argv[5];
//std::string trdnumber = argv[6];
std::string suffix = ".txt";
std::string cccuttext = "cccut";
std::string ccnumbertext = "CCN";
std::string trdnumbertext = "TRDN";
std::string fit_results = fitresultname + cccuttext + cccut + ccnumbertext + ccnumber + trdnumbertext + trdnumber + suffix;
std::string fit_results_error = fiterrorname + cccuttext + cccut + ccnumbertext + ccnumber + trdnumbertext + trdnumber + suffix;
std::string chi2 = chi2name + cccuttext + cccut + ccnumbertext + ccnumber + trdnumbertext + trdnumber + suffix;
std::string s3 = s1 + se + cccuttext + cccut + ccnumbertext + ccnumber + trdnumbertext +  trdnumber + s2;
const char* filename = (char*)s3.c_str();

//double TRDlow = atof(argv[2]);
//double TRDhigh = atof(argv[3]);
//double CClow = atof(argv[4]);
std::string TRDlow_pro_text = "_trdlow_pro:";
std::string CClow_pro_text = "_cclow_pro:";

cout << filename << endl;

TFile *f = new TFile(filename);
TH2D *template_correct = (TH2D*)f->Get("template_correct");
TH2D *template_confused = (TH2D*)f->Get("template_confused");
TH2D *template_electron = (TH2D*)f->Get("template_electron");
TH2D *data = (TH2D*)f->Get(argv[1]); 

MYUtilities::TemplateFitter2D templateFitter(0);
templateFitter.AddTemplateHistogram(template_correct);
templateFitter.AddTemplateHistogram(template_confused);
templateFitter.AddTemplateHistogram(template_electron);
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
 
cout << result << endl;
cout << ResultError << endl;
cout << Chi2 << endl;
cout << NDF << endl;
cout << CHI2dof << endl;

chdir("..");
if(NULL==opendir("/results"))
   mkdir("./results",0777);
chdir("./results");

ofstream fitresults(fit_results, ios::app);
if (fitresults.is_open())
        {
      fitresults << result << endl;
      fitresults.close();
        }
ofstream fiterror(fit_results_error, ios::app);
if (fiterror.is_open())
        {
      fiterror << ResultError << endl;
      fiterror.close();
        }
ofstream fitchi2(chi2, ios::app);
if (fitchi2.is_open())
        {
      fitchi2 << CHI2dof << endl;
      fitchi2.close();
        }
chdir("..");
/*
//////////////////////// create Fit Result and Projections plots //////////////////////////////////
chdir("./fit_plots");
printf("current working directory: %s\n", getcwd(NULL, NULL));

std::string name1 = "Fit_Result";
std::string name2 = "CCprojection";
std::string name3 = "TRDprojection";
std::string name4 = se;
std::string name5 = ".pdf";
std::string namefit = name1 + name4 + cccuttext + cccut + ccnumbertext + ccnumber + trdnumbertext + trdnumber + name5;
std::string nameccprojection = name2 + name4 + cccuttext + cccut + ccnumbertext + ccnumber + trdnumbertext + trdnumber + TRDlow_pro_text + argv[2] + name5;
std::string nameTRD = name3 + name4 + cccuttext + cccut + ccnumbertext + ccnumber + trdnumbertext + trdnumber + CClow_pro_text + argv[4]  +name5;

TCanvas * c1 = templateFitter.CreateResultDrawing("Fit_Result",800,600);
c1->SaveAs((char*)namefit.c_str());
TCanvas * c2 = templateFitter.CreateResultDrawingXprojection("CCprojection",800,600,TRDlow,TRDhigh);
c2->SaveAs((char*)nameccprojection.c_str());
TCanvas * c3 = templateFitter.CreateResultDrawingYprojection("TRDprojection",800,600,CClow,1.0);
c3->SaveAs((char*)nameTRD.c_str());

chdir("..");
*/
      }
    }
  } 
  return EXIT_SUCCESS;
}
