#ifndef RooExpGausExpPdf_hh
#define RooExpGausExpPdf_hh

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

class RooExpGaussExpPdf : public RooAbsPdf {
public:
  RooExpGaussExpPdf() { }
  RooExpGaussExpPdf(const char *name, const char *title,
                    RooAbsReal& _x,
                    RooAbsReal& _mu, RooAbsReal& _sigma,
                    RooAbsReal& _kLow, RooAbsReal& _kHigh);

  RooExpGaussExpPdf(const RooExpGaussExpPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { return new RooExpGaussExpPdf(*this, newname); }
  virtual ~RooExpGaussExpPdf() { }

protected:
  RooRealProxy x;
  RooRealProxy mu;
  RooRealProxy sigma;
  RooRealProxy kLow;
  RooRealProxy kHigh;
  virtual Double_t evaluate() const;

private:
  ClassDef(RooExpGaussExpPdf, 1)
};

#endif
