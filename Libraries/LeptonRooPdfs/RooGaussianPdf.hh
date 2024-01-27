#ifndef RooGaussianPdfPdf_hh
#define RooGaussianPdfPdf_hh

#if defined(__GNUC__) &&  __GNUC__ >= 6
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Waddress"
#pragma GCC diagnostic ignored "-Wnonnull-compare"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooRealVar.h>
#pragma GCC diagnostic pop
#else
#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooRealVar.h>
#endif

class RooGaussianPdf : public RooAbsPdf {
public:
  RooGaussianPdf() { }
  RooGaussianPdf(const char* name, const char* title, RooAbsReal& _x, RooAbsReal& _mean, RooAbsReal& _sigma);
  RooGaussianPdf(const RooGaussianPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { return new RooGaussianPdf(*this, newname); }
  virtual ~RooGaussianPdf() {  }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

protected:
  RooRealProxy x;
  RooRealProxy mean;
  RooRealProxy sigma;

  virtual Double_t evaluate() const;

private:
  ClassDef(RooGaussianPdf, 1)
};

#endif
