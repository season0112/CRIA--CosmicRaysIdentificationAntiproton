#ifndef RooLangausPdf_hh
#define RooLangausPdf_hh

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

class RooLangausPdf : public RooAbsPdf {
public:
  RooLangausPdf() { }
  RooLangausPdf(const char *name, const char *title,
		 RooAbsReal& _x, RooAbsReal& _widthLandau, RooAbsReal& _mpvLandau, RooAbsReal& _totalArea, RooAbsReal& _sigmaGauss);

  RooLangausPdf(const RooLangausPdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { return new RooLangausPdf(*this, newname); }
  virtual ~RooLangausPdf() { }

protected:
  RooRealProxy x;
  RooRealProxy widthLandau;
  RooRealProxy mpvLandau;
  RooRealProxy totalArea;
  RooRealProxy sigmaGauss;
  virtual Double_t evaluate() const;

private:
  ClassDef(RooLangausPdf, 1)
};

#endif
