#include "TemplateFitter2DforProtonCC_CCfixed_dev.hh"

#include <cassert>
#include <TH2.h>
#include <TMath.h>

TemplateFit2DFitFunction::TemplateFit2DFitFunction(const TH2* data,
                                                   const TH2* template_correct, 
                                                   const TH2* template_confused,
                                                   const TH2* template_electron_Data)
  : fdata(data)
  , ftemplate_correct(template_correct)
  , ftemplate_confused(template_confused)
  , ftemplate_electron_Data(template_electron_Data) {

  fNumberOfEvents = data->Integral();
}

double TemplateFit2DFitFunction::PdfValue(int xBin, int yBin) const {

  double nCorrect   = fCorrectFraction  * fNumberOfEvents;
  double nConfused  = fConfusedFraction * fNumberOfEvents;
  double nElectron  = fElectronFraction * fNumberOfEvents;

  double pdf = nCorrect           *  ftemplate_correct->GetBinContent(xBin, yBin) 
             + nConfused          *  ftemplate_confused->GetBinContent(xBin, yBin)
             + nElectron          *  ftemplate_electron_Data->GetBinContent(xBin, yBin);

  assert(fNumberOfEvents >= 0);
  assert(!std::isinf(pdf));
  assert(!std::isnan(pdf));

  return pdf / fNumberOfEvents;
}

double TemplateFit2DFitFunction::DoEval(const double* parameters) const {

  const_cast<TemplateFit2DFitFunction*>(this)->UpdateParameters(parameters);

  double minusLogL = 0.0;
  for (int i = 1 ; i <= fdata->GetNbinsX(); ++i) {
    for (int j = 1; j <= fdata->GetNbinsY(); ++j) {
      int binContent = TMath::Nint(fdata->GetBinContent(i, j));
      if (!binContent)
        continue;

      double pdfValue = PdfValue(i, j);
      if (pdfValue <= 0.0)
        continue;

      minusLogL += (binContent * (-std::log(pdfValue)));
    }
  }

  // term from extended ML, to take normalization into account
  // http://hepunx.rl.ac.uk/~adye/thesis/html/node49.html
  // or Cowan pp 84f
  minusLogL += (fCorrectFraction + fConfusedFraction + fElectronFraction) * fNumberOfEvents;

  return minusLogL;
}

TH2* TemplateFit2DFitFunction::CalculateHistogram(int contributingCategories) const {

  double nCorrect   = fCorrectFraction  * fNumberOfEvents;
  double nConfused  = fConfusedFraction * fNumberOfEvents;
  double nElectron  = fElectronFraction * fNumberOfEvents;

  TH2* histogram = dynamic_cast<TH2*>(fdata->Clone(Form("hTemplateFitResult2D_%d", contributingCategories)));
  histogram->Reset();
  histogram->SetTitle("2D template fit result");
  histogram->SetStats(0);

  for (int xBin = 1; xBin <= histogram->GetNbinsX(); ++xBin) {
    for (int yBin = 1; yBin <= histogram->GetNbinsY(); ++yBin) {
      double prediction = 0.0;
      if (contributingCategories & int(ContributingCategories::Correct))
        prediction += nCorrect * ftemplate_correct->GetBinContent(xBin, yBin);

      if (contributingCategories & int(ContributingCategories::Confused))
        prediction += nConfused * ftemplate_confused->GetBinContent(xBin, yBin);

      if (contributingCategories & int(ContributingCategories::Electron))
        prediction += nElectron * ftemplate_electron_Data->GetBinContent(xBin, yBin);

      histogram->SetBinContent(xBin, yBin, prediction);
    }
  }

  return histogram;
}

TH2* TemplateFit2DFitFunction::CalculateTotalHistogram() const {

  int total = int(ContributingCategories::Correct)
            | int(ContributingCategories::Confused)
            | int(ContributingCategories::Electron);

  return CalculateHistogram(total);
}

TH2* TemplateFit2DFitFunction::DrawZeroProbabilityEvents() const {

  TH2* histogram = dynamic_cast<TH2*>(fdata->Clone("hZeroProbability"));
  histogram->Reset();
  histogram->SetStats(0);

  for (int xBin = 1; xBin <= histogram->GetNbinsX(); ++xBin) {
    for (int yBin = 1; yBin <= histogram->GetNbinsY(); ++yBin) {
      double pdfCorrect = ftemplate_correct->GetBinContent(xBin, yBin);
      double pdfConfused = ftemplate_confused->GetBinContent(xBin, yBin);
      double pdfElectron = ftemplate_electron_Data->GetBinContent(xBin, yBin);

      if (pdfCorrect <= 0.0 && pdfConfused <= 0.0 && pdfElectron <= 0.0)
        histogram->SetBinContent(xBin, yBin, fdata->GetBinContent(xBin, yBin));
    }
  }

  histogram->SetTitle(Form("%g events with zero probability", histogram->Integral()));
  return histogram;
}

void TemplateFit2DFitFunction::UpdateParameters(const double* parameters) {

  fCorrectFraction = parameters[0];
  fConfusedFraction = parameters[1];
  fElectronFraction = parameters[2];
  fChargeConfusionProbability = parameters[3];
}
