#ifndef RooAsymmLandauLikePdf_hh
#define RooAsymmLandauLikePdf_hh

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

class RooAsymmLandauLikePdf : public RooAbsPdf {
public:
  RooAsymmLandauLikePdf() { }
  RooAsymmLandauLikePdf(const char *name, const char *title,
		        RooAbsReal& _x, RooAbsReal& _peak, RooAbsReal& _width, RooAbsReal& _asymm);

  RooAsymmLandauLikePdf(const RooAsymmLandauLikePdf& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const { return new RooAsymmLandauLikePdf(*this, newname); }
  virtual ~RooAsymmLandauLikePdf() { }

protected:
  RooRealProxy x;
  RooRealProxy peak;
  RooRealProxy width;
  RooRealProxy asymm;
  virtual Double_t evaluate() const;

private:
  ClassDef(RooAsymmLandauLikePdf, 1)
};

#endif
