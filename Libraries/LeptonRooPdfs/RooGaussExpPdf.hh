#ifndef RooGausExpPdf_hh
#define RooGausExpPdf_hh

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

class RooGaussExpPdf : public RooAbsPdf {
public:
  RooGaussExpPdf() { }
  RooGaussExpPdf(const char *name, const char *title,
		 RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _k);

  RooGaussExpPdf(const RooGaussExpPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { return new RooGaussExpPdf(*this, newname); }
  virtual ~RooGaussExpPdf() { }

protected:
  RooRealProxy x;
  RooRealProxy mu;
  RooRealProxy sigma;
  RooRealProxy k;
  virtual Double_t evaluate() const;

private:
  ClassDef(RooGaussExpPdf, 1)
};

#endif
