#ifndef TemplateFitter2DforProtonCC_CCfixed_dev_hh
#define TemplateFitter2DforProtonCC_CCfixed_dev_hh

#include <TMinuitMinimizer.h>

class TH2;

class TemplateFit2DFitFunction : public ROOT::Math::IMultiGenFunction {
public:
  enum class ContributingCategories {
    Correct     = 1 << 1,
    Confused    = 1 << 2,
    Electron    = 1 << 3,
  };

  TemplateFit2DFitFunction(const TH2* data,
                           const TH2* template_correct,
                           const TH2* template_confused,
                           const TH2* template_electron_Data);

  /** Number of parameters */
  virtual unsigned int NDim() const { return 5; }  

  /** Defines the actual function call. */
  virtual double DoEval(const double* parameters) const;

  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const { return new TemplateFit2DFitFunction(*this); }

  double PdfValue(int xbin, int ybin) const;
  TH2* CalculateHistogram(int contributingCategories) const;
  TH2* CalculateTotalHistogram() const;
  TH2* DrawZeroProbabilityEvents() const;

  void UpdateParameters(const double*);

private:
  const TH2* fdata { nullptr };

  const TH2* ftemplate_correct { nullptr };
  const TH2* ftemplate_confused { nullptr };
  const TH2* ftemplate_electron_Data { nullptr }; 

  double fNumberOfEvents { 0.0 };
  double fCorrectFraction { 0.0 }; 
  double fConfusedFraction { 0.0 };
  double fElectronFraction { 0.0 };
  double fChargeConfusionProbability { 0.0 };

};

#endif
