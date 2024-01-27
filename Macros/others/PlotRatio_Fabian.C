#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include "BinningDefinition.hh"
#include "BinningTools.hh"
#include "Utilities.hh"

int PlotRatio_Fabian() {

 TFile* fFile = TFile::Open("/home/bo791269/Software/AntiprotonAnalysis/ReferenceFiles/TimeDependentAnalysis/LeptonFluxes_Zimmermann_04_09_2018_2D.root", "OPEN");
  assert(fFile);
  TH2D* histoBinning = static_cast<TH2D*>(fFile->Get("hAmsSecondsAboveGeoCutoffVsTime"));
  assert(histoBinning);
  auto energyBinning = Binning::FromTAxis(histoBinning->GetYaxis());
  auto timeBinning = Binning::FromTAxis(histoBinning->GetXaxis());
  TGraphAsymmErrors** grRatio = new TGraphAsymmErrors*[timeBinning.NumberOfBins()];
  for (unsigned int i = 0; i < timeBinning.NumberOfBins(); ++i) {
    std::stringstream ss;
    ss << (int)timeBinning.LowEdge(i + 1) << "_" << (int)timeBinning.UpEdge(i + 1);
    grRatio[i] = static_cast<TGraphAsymmErrors*>(fFile->Get(("SxRatioTools_PosiOverElec_" + ss.str() + "/grPosiOverElec_" + ss.str() + "_PosiOverElec").c_str()));
    assert(grRatio[i]);
    TGraph* grRelErrStat = static_cast<TGraphAsymmErrors*>(fFile->Get(("SxRatioTools_PosiOverElec_" + ss.str() + "/grPosiOverElec_" + ss.str() + "_RelErrorStat").c_str()));
    assert(grRelErrStat);
    TGraph* grRelErrSys = static_cast<TGraphAsymmErrors*>(fFile->Get(("SxRatioTools_PosiOverElec_" + ss.str() + "/grPosiOverElec_" + ss.str() + "_RelErrorUncorrSyst").c_str()));
    assert(grRelErrSys);
    for (int j = 0; j < grRatio[i]->GetN(); ++j) {
      const double error = Utilities::QuadraticSum(grRelErrStat->GetY()[j], grRelErrSys->GetY()[j]) * grRatio[i]->GetY()[j];
      grRatio[i]->SetPointError(j, grRatio[i]->GetEXlow()[j], grRatio[i]->GetEXhigh()[j], error, error);
    }
  }

  const int xbin = energyBinning.FindBin(5.2) - 1; // example energy bin to plot
 
  cout<< xbin <<endl;
 
  TGraphAsymmErrors* gr = new TGraphAsymmErrors(timeBinning.NumberOfBins());
  for (unsigned int i = 0; i < timeBinning.NumberOfBins(); ++i) {
    assert(grRatio[i]->GetX()[xbin] > 5.0 && grRatio[i]->GetX()[xbin] < 5.49);
    gr->SetPoint(i, timeBinning.Value(i + 1), grRatio[i]->GetY()[xbin]);
    gr->SetPointError(i, 0.0, 0.0, grRatio[i]->GetEYlow()[xbin], grRatio[i]->GetEYhigh()[xbin]);
    
  }

  TCanvas* c = new TCanvas();
  gr->SetMarkerStyle(20);
  gr->SetTitle("Positron to Electron Ratio from 5.0 to 5.49 GeV");
  gr->GetYaxis()->SetTitle("R_{e}");
  gr->GetYaxis()->SetRangeUser(0.04, 0.08);
  gr->Draw("AP");
  c->Print("Ratio_vs_Time.pdf");

  return 0;
}


