#include "UnfoldingSystematicErrors.hh"

#include <TGraph.h>

#include <Utilities.hh>

#define INFO_OUT_TAG "UnfoldingSystematicErrors"
#include <debugging.hh>

TGraph* UnfoldingSystematicErrors::MakeUncertaintyFromToymcBias() {

  static const int npe = 8;
  Double_t ex1[npe] = {
    0.5,
    0.75,
    1.0,
    1.5,
    2.0,
    100.0,
    500.0,
    1000.0
  };

  Double_t ey1[npe] = {
    0.35,
    0.04,
    0.01,
    0.004,
    0.002,
    0.002,
    0.004,
    0.01
  };

  TGraph* graphError = new TGraph(npe, ex1, ey1);
  graphError->SetName("grSystErrorFromToyMcBias");

  return graphError;
}

TGraph* UnfoldingSystematicErrors::MakeUncertaintyFromTestbeam() {

  static const int npe = 8;
  Double_t ex1[npe] = {
    0.5,
    1.0,
    2.0,
    3.0,
    20.0,
    50.0,
    400.0,
    1000.0
  };

  Double_t ey1[npe] = {
    0.03,
    0.01,
    0.004,
    0.002,
    0.002,
    0.002,
    0.002,
    0.002
  };

  TGraph* graphError = new TGraph(npe, ex1, ey1);
  graphError->SetName("grSystErrorFromTestbeam");

  return graphError;
}

TGraph* UnfoldingSystematicErrors::MakeTotalSystErrorFromUnfolding(const std::vector<double>& xPoints) {

  TGraph* grTotalErr = new TGraph;

  TGraph* grToyMcBias = MakeUncertaintyFromToymcBias();
  TGraph* grMigr = MakeUncertaintyFromTestbeam();

  for (double x : xPoints) {

    double totalErr = Utilities::QuadraticSum(Utilities::LogarithmicEval(grToyMcBias, x), Utilities::LogarithmicEval(grMigr, x));
    grTotalErr->SetPoint(grTotalErr->GetN(), x, totalErr);
  }

  delete grToyMcBias;
  delete grMigr;

  return grTotalErr;
}
