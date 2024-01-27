#ifndef RooNovosibirskPdf_hh
#define RooNovosibirskPdf_hh

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

class RooNovosibirskPdf : public RooAbsPdf {
public:
  RooNovosibirskPdf() { }
  RooNovosibirskPdf(const char* name, const char* title, RooAbsReal& _x, RooAbsReal& _mu, RooAbsReal& _sigma, RooAbsReal& _tau);
  RooNovosibirskPdf(const RooNovosibirskPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { return new RooNovosibirskPdf(*this, newname); }
  virtual ~RooNovosibirskPdf() { }

protected:
  RooRealProxy x;
  RooRealProxy mu;
  RooRealProxy sigma;
  RooRealProxy tau;
  virtual Double_t evaluate() const;

private:
  ClassDef(RooNovosibirskPdf, 3)
};

#endif
