#ifndef CompositeParameterizationFunction_hh
#define CompositeParameterizationFunction_hh

class TF1;

class CompositeParameterizationFunction {

public:
  CompositeParameterizationFunction(double Emin, double Emax);
  ~CompositeParameterizationFunction();

  double Evaluate(double* x, double* p);

  void SetUnderOverflow(double underflowProb, double overflowProb);

private:
  void SetParameters(double landauMpv, double landauSigma, double cbAlpha, double cbN, double cbMu, double cbSigma, double landauFraction);

  double fLandauMpv = 0.0;
  double fLandauSigma = 0.0;
  double fCBAlpha = 0.0;
  double fCBN = 0.0;
  double fCBMu = 0.0;
  double fCBSigma = 0.0;
  double fLandauFraction = 0.0; // integral over Landau part of pdf

  double fUnderflowProb = 0.0;
  double fOverflowProb = 0.0;

  double fEmin = 0.0;
  double fEmax = 0.0;

  TF1* fLandauPart = nullptr; // helper function for Landau integration
};

#endif
