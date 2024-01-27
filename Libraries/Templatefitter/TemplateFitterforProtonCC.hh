#ifndef TemplateFitterforProtonCC_hh
#define TemplateFitterforProtonCC_hh

#include <TObject.h>
#include <string>
#include <vector>

class TH1F;
class TH1;
class TCanvas;
class TMinuit;
class TH2D;

namespace MYUtilities {

/** The TemplateFitter class
  *
  * Class to fit an arbitrary number of template distributions to a data distribution.
  *
  * Variable bin widths are fine.
  *
  * You can use a Chisquared or Likelihood fit. The Chisquared fit also works fine if you have filled your histograms with SetBinContent and SetBinError.
  * The likelihood fit will always assume Poisson errors on the BinContents, no matter how the data and template histos are filled.
  */
class TemplateFitter {

public:
  /** Default constructor
    *
    * \param printlevel TMinuit output: 0 standard, -1 quiet
    */
  TemplateFitter(int printlevel);

  /** Default destructor
    *
    * calls Clear()
    */
  virtual ~TemplateFitter();

  /** Find out if fit has been performed successfully */
  bool HasBeenFit() const { return fHasBeenFit; }

  /** Re-sets the printlevel of Minuit after initialization */
  void SetPrintLevel(int pl) { fPrintLevel = pl; }

  /** Set Data Histogram */
  void SetDataHistogram(const TH1* hist) { fInputData = hist; }

  /** Data Histogram. */
  const TH1* DataHistogram() const { return fInputData; }

  /** Index of histogram bin at the lower fit range. */
  int FitRangeLow() const { return fFitRangeLow; }

  /** Index of histogram bin at the upper fit range. */
  int FitRangeHigh() const { return fFitRangeHigh; }

  /** Number of degrees of freedom. */
  int NDF() const { return fNDF; }

  /** Chi square of the fit. */
  double Chi2() const { return fChi2; }

  /** Fit result histogram. */
  const TH1* FitResult() const { return fFitResult; }

  /** Histograms of fit results of individual template distributions. */
  const std::vector<TH1*>& ResultHistograms() const { return fResultHistos; }

  /** Add Template Histogram
    *
    * Add as many templates as you wish. The parameter order is the same as the histogram order, i.e. the first fit parameter belongs to the first template shape.
    * You can also provide a starting value for the fit to be used with this template.
    */
  virtual void AddTemplateHistogram(const TH1* hist, double start = -1.0);

  /**
    * Set Start Value for the fit for template i.
    * Attention: Start Values are relative! Must be between 0 and 1.
    */
  void SetStartValue(unsigned int i, double start);

  /** Set start values automagically.
    *
    * Algorithm:
    *  - assume templates have (nearly) symmetrical shape around the respective mean values
    *  - and overlap between templates beyond mean of the others is small
    *  - then we can place two cuts, one at the mean of each template, count the number of events beyond each cut and multiply by 2
    *  - this allows us to calculate an approximate fraction for each contribution
    *  - in the end, rescale fraction such that their sum is unity
    *
    * \attention Call this function right before the call to Fit().
    */
  bool AutoStartValues();

  /** Fix Parameter i to the value start
    *
    * \attention start is relative and must be between 0 and 1.
    */
  void FixParameter(unsigned int i, double start);

  /** Free parameter i.
    *
    * The start value set before is kept. To reset the start value, use SetStartValue or ResetStartValue
    */
  void FreeParameter(unsigned int i);

  /**
    * Reset Start Value for parameter i. This "frees" the start value and it is set to the default (1/Ntemplates).
    */
  void ResetStartValue(unsigned int i);

  void UseUnboundedFit(bool val = true) { fUnboundedFit = val; }

  /** Get Relative Result of fit
    *
    * This represents the relative amount of a template form found in the data.
    * Example: If the first entry is 0.2, then 20% of the data distribution is due to first template.
    * The sum of all sould add up to unity within uncertainties.
    */
  const std::vector<double>& GetRelativeResult() const;

  /** Get Absolute Result of fit
    *
    * This represents the absolute amount of a template form found in the data.
    * Example: If the first entry is 700, then 700 events in the data are from first template.
    *
    * This getter respects the fit ranges for 1D fits: If the first entry is 700, then 700 events in the data within the fit range are from first template.
    */
  const std::vector<double>& GetAbsoluteResult() const;

  /** Get Absolute Result of fit considering the whole histogram range.
    *
    * This represents the absolute amount of a template form found in the data.
    * Example: If the first entry is 700, then 700 events in the data are from first template.
    *
    * This getter always considers the whole template range: Even though only a sub-range is fit, if the first entry is 700, then 700 events in the whole data are from first template.
    */
  const std::vector<double>& GetAbsoluteFullRangeResult() const;


  /** Get MIGRAD uncertainty of relative result of fit.
    *
    * This is the symmetric MIGRAD uncertainty of the values of GetRelativeResult
    */
  const std::vector<double>& GetRelativeResultError() const;

  /** Get MIGRAD uncertainty of absolute result of fit.
    *
    * This is the symmetric MIGRAD uncertainty of the values of GetAbsoluteResult
    */
  const std::vector<double>& GetAbsoluteResultError() const;

  /** Get MIGRAD uncertainty of absolute result of fit over the whole range.
    *
    * This is the symmetric MIGRAD uncertainty of the values of GetAbsoluteFullRangeResult
    */
  const std::vector<double>& GetAbsoluteFullRangeResultError() const;

  /** Get upper uncertainty of relative result of fit. This is the MINOS error.
    *
    * This is the MINOS uncertainty of the values of GetRelativeResult
    */
  const std::vector<double>& GetRelativeResultErrorUpper() const;

  /** Get upper uncertainty of absolute result of fit. This is the MINOS error.
    *
    * This is the MINOS uncertainty of the values of GetAbsoluteResult
    */
  const std::vector<double>& GetAbsoluteResultErrorUpper() const;

  /** Get upper uncertainty of absolute result of fit over the whole range. This is the MINOS error.
    *
    * This is the MINOS uncertainty of the values of GetAbsoluteFullRangeResult
    */
  const std::vector<double>& GetAbsoluteFullRangeResultErrorUpper() const;


  /** Get lower uncertainty of relative result of fit. This is the MINOS error.
    *
    * This is the MINOS uncertainty of the values of GetRelativeResult
    */
  const std::vector<double>& GetRelativeResultErrorLower() const;

  /** Get lower uncertainty of absolute result of fit. This is the MINOS error.
    *
    * This is the MINOS uncertainty of the values of GetAbsoluteResult
    */
  const std::vector<double>& GetAbsoluteResultErrorLower() const;

  /** Get lower uncertainty of absolute result of fit over the whole range. This is the MINOS error.
    *
    * This is the MINOS uncertainty of the values of GetAbsoluteFullRangeResult
    */
  const std::vector<double>& GetAbsoluteFullRangeResultErrorLower() const;


  /** Get global correlation factors as determined by TMinuit from the error matrix. */
  const std::vector<double>& GetGlobalCorrelationFactors() const;

  /** Get MIGRAD error matrix.
    *
    * The returned vector has N*N elements for N fitted templates. The index of
    * the matrix entry (row,col) in the vector is [row*N + col].
    */
  const std::vector<double>& GetErrorMatrix() const;

  /** Clone fit result
    *
    * Returns the resulting histogram of the fit as a clone. This is the weighted sum of the templates that describe the data best.
    * The user has to delete this, the TemplateFitter doesn't
    */
  TH1F* CloneCombinedResultHistogram(const char*) const;


  /** Create Canvas with result drawn
    *
    * Parameters are canvasname, width, and height.
    *
    * Returns a canvas with the fit result drawn into it. The user has to delete this, the class doesn't.
    * Supply a canvas name to prevent canvases with one default name to be overwritten.
    */
  TCanvas* CreateResultDrawing(std::string canvasname, int width, int height, bool defaultStyle = true) const;

  /** Generates a contour plot of the two specified parameters of the fit.
    *
    * It is drawn for the parameters par1 and par2.
    * The function returns a TCanvas object.
    * The contours of the 1-, 2-, and 3-sigma uncertainties are drawn.
    */
  TCanvas* CreateContourPlot(int par1, int par2, std::string canvasname, int width, int height) const;

  /** Clear function
    *
    * Resets all input histograms. Use this function before you add new data or template histograms.
    */
  void Clear();

  /** Perform fit
    *
    * Calls the fit appropriate fit function. You can choose with the parameter i between chi-square (i==0) or likelihood fit (i==1)(default).
    *
    * The two doubles give the fit range to be used. Fit range is only available in 1D case; the 2D fitter ignores ranges.
    * If you specify ranges, please make sure to select the correct AbsoluteResult for your application:
    * GetAbsoluteResult() returns the numbers of events within the range.
    * GetAbsoluteFullRangeResult() returns the numbers of events over the whole range, despite the fit considering only a range.
    *
    */
  virtual double Fit(int i = 1, double lowBound = 0, double highBound = 0);


  /** Internal fit function for Chisquared fit. */
  static void chisqr_fit_func(int& npar, double* gin, double& f, double* par, int iflag);

  /** Internal fit function for Likelihood fit. */
  static void llh_fit_func(int& npar, double* gin, double& f, double* par, int iflag);

protected:
  struct TemplateFitObject : public TObject {

    std::vector<std::vector<double>> vec;
  };

  int TemplateColor(int i) const;
  int TemplateFillStyle(int i) const;

  /** Tests the compatibility of histograms
    *
    * Checks if the data and template histograms are compatible. They need the same number of bins, binwidth, and bin edges.
    */
  virtual bool TestHistogramCompatibility() const;

  /** Internal routine to calculate Chisquare values.
    *
    * This is used to determine goodness-of-fit for the likelihood fit and to fill the contour plot.
    */
  double ManualChiSqr(int& NDF, std::vector<double> par) const;

  bool fHasBeenFit; ///< Flag if fit has been done

  const TH1* fInputData;                   ///< Container for data histogram that is given to Fitter.
  std::vector<const TH1*> fInputTemplates; ///< Container for templates that are given to Fitter.
  std::vector<double> fStartValues;        ///< Container for start values for the fit.
  std::vector<bool> fFixParameter;         ///< Container of flags to fix parameters

  TH1* fFitResult;                 ///< Histogram of fit result combined as sum from all templates
  std::vector<TH1*> fResultHistos; ///< Histograms of fit results of individual template distributions
  std::vector<double> fRelative;   ///< Fit Result in relative amounts

private:
  /** Internal, private Routine to execute the fit either by Chisquared or -Likelihood minimization.
    * Before the fit a check on consistency of the histograms is done (Nbins, BinWidth,BinEdges)
    *
    * Option 0 in Fit():
    *
    * Applies a Chisquared fit to match the data histogram with the template histograms and free amplitudes for the templates.
    * If you used SetBinError to create your histograms, use this function!
    *
    * Option 1 in Fit():
    *
    * Fit routine using a likelihood function to determine the parameters. The fit function
    * is taken from Cowan, "Statistical Data Analysis", p. 89, Eq. 6.46. It always assumes
    * Poisson statistics in Data and Template Histograms.
    *
    * The fit results can afterwards be accessed by getter functions.
    * ATTENTION: If you specify a range, please choose the appropriate getter function for your application!
    * Relative weights as well as absolute values are provided.
    *
    * The function returns the reduced Chisquared.
    */
  double ExecuteFit(void (*function)(int&, double*, double&, double*, int));
  bool CheckFit() const;

  /** Internal routine to translate your input histograms in data that the fit uses.
    *
    * Here, fDataSet and fRelError are filled. The relative error is the poisson error, or, if SetBinError was used, the set bin error is filled.
    * \attention No matter what you do, the likelihood fit will use poisson errors for data and templates! In case you use the likelihood fit with templates prepared with SetBinError, the returned uncertainty is not the one you want. And the contour plot will be not the one you want.
    * IF YOU USED SetBinError TO ASSEMBLE THE TEMPLATES USE THE chisquare Fit (Option 0 in the Fit() routine)
    *
    */
  void PrepareFitData();

  TMinuit* fMinimizer; ///< The ROOT Minuit Minimizer Object. It does all the hard work.
  int fPrintLevel;     ///< Print Level for Minuit output
  bool fUnboundedFit;  ///< Do not constrain relative amplitudes to [0,1]

  int fFitRangeLow;  ///< Index of histogram bin at the lower fit range.
  int fFitRangeHigh; ///< Index of histogram bin at the upper fit range.

  int fNDF;     ///< Number of degrees of freedom
  double fChi2; ///< Chi square of the fit

  std::vector<double> fAbsolute;              ///< Fit result in absolute numbers
  std::vector<double> fAbsoluteAllRange;      ///< Fit result in absolute numbers as if no fit range was specified.
  std::vector<double> fRelativeError;         ///< Uncertainty of fit result in relative amounts (MIGRAD)
  std::vector<double> fAbsoluteError;         ///< Uncertainty of fit result of absolute numbers (MIGRAD)
  std::vector<double> fAbsoluteAllRangeError; ///< Uncertainty of fit result of absolute numbers (MIGRAD) as if no fit range was specified.

  std::vector<double> fRelativeErrorUpper;         ///< Uncertainty of fit result in relative amounts (upper MINOS uncertainty)
  std::vector<double> fAbsoluteErrorUpper;         ///< Uncertainty of fit result of absolute numbers (upper MINOS uncertainty)
  std::vector<double> fAbsoluteAllRangeErrorUpper; ///< Uncertainty of fit result of absolute numbers (upper MINOS uncertainty) as if no fit range was specified.
  std::vector<double> fRelativeErrorLower;         ///< Uncertainty of fit result in relative amounts (lower MINOS uncertainty)
  std::vector<double> fAbsoluteErrorLower;         ///< Uncertainty of fit result of absolute numbers (lower MINOS uncertainty)
  std::vector<double> fAbsoluteAllRangeErrorLower; ///< Uncertainty of fit result of absolute numbers (lower MINOS uncertainty) as if no fit range was specified.

  std::vector<double> fGlobalCorrelationFactor; ///< Global Correlation Factor as calculated from MIGRAD error matrix
  std::vector<double> fErrorMatrix;             ///< MIGRAD error matrix

  std::vector<double> fBinWidth;              ///< Container used in fit function: Bin width (1D) or bin area (2D)
  std::vector<double> fNentries;              ///< Container used in fit function: Number of entries in data and template fits
  std::vector<std::vector<double>> fDataSet;  ///< Container used in fit function: Data of data and template histograms
  std::vector<std::vector<double>> fRelError; ///< Container used in fit function: Relative Error on entries in data and template histograms. Poisson error, except if SetBinError was used
};

} // namespace MYUtilities

#endif
