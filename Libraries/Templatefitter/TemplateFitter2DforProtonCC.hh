#ifndef TemplateFitter2DforProtonCC_hh
#define TemplateFitter2DforProtonCC_hh

#include "TemplateFitterforProtonCC.hh"

#include <string>
#include <vector>

class TH2F;
class TH1F;
class TH1;
class TH2;
class TCanvas;

namespace MYUtilities {

/** The TemplateFitter2D class
  *
  * Class to fit an arbitrary number of template distributions to a data distribution in two dimensions.
  *
  * You have to use TH2F. Variable bin widths are fine.
  *
  * NO FIT RANGES ARE SUPPORTED
  *
  * You can use a Chisquared or Likelihood fit.
  *
  */
class TemplateFitter2D : public TemplateFitter {

public:
  /** Default constructor
    *
    * \param printlevel for TMinuit output: 0 standard, -1 quiet
    */
  TemplateFitter2D(int printlevel);
  ~TemplateFitter2D();

  /** Add Template Histogram
    *
    * Add as many templates as you wish. The parameter order is the same as the histogram order, i.e. the first fit parameter belongs to the first template shape.
    * You can also provide a starting value for the fit to be used with this template.
    */
  virtual void AddTemplateHistogram(const TH1* hist, double start = -1.0);


  /** Clone fit result
    *
    * Returns the resulting histogram of the fit as a clone. This is the weighted sum of the templates that describe the data best.
    * The user has to delete this, the TemplateFitter doesn't
    */
  TH2F* CloneCombinedResultHistogram(const char*) const;

  /** Function to get the X projection of the resulting fit
    *
    * Returned is a TH1F* histogram. By default the whole Y range is used to generate the projections.
    * Ranges can be given for creating the projections.
    * Use the parameter hname to give the resulting histogram a unique identifier (Name)
    */
  TH1F* GetResultHistogramXprojection(const char*, double Ylow = 0, double Yhigh = 0) const;

  /** Function to get the Y projection of the resulting fit
    *
    * Returned is a TH1F* histogram. By default the whole X range is used to generate the projections.
    * Ranges can be given for creating the projections.
    * Use the parameter hname to give the resulting histogram a unique identifier (Name)
    */
  TH1F* GetResultHistogramYprojection(const char*, double Xlow = 0, double Xhigh = 0) const;


  /** Create Canvas with result drawn
    *
    * \param canvasname
    *
    * Returns a canvas with the fit result drawn into it. The user has  to delete this, the class doesn't.
    * Supply a canvas name to prevent canvases with one default name to be overwritten.
    *
    * The resulting canvas is filled with the data histogram, the fit result, the residual histogram, and fit information
    * The residual histogram is the quadratic difference between fit and data divided by the square of the data uncertainty. It assumes that the template has no uncertainty!
    */
  TCanvas* CreateResultDrawing(std::string canvasname, int width, int height) const;

  /** This function generates a TCanvas with the X projection of the fit result drawn into it.
    *
    * The canvas needs to be given a name as identifier as well as height and width.
    * By default the whole Y range is used to produce the projections, however you can also specify a range.
    *
    * The layout of the canvas follows the default drawing style of the 1D fitter.
    */
  TCanvas* CreateResultDrawingXprojection(std::string canvasname, int width, int height, double Ylow = 0, double Yhigh = 0, int BinRemovedNumber_FromLeft = 0, int RebinNumber_X = 1) const;

  /** This function generates a TCanvas with the Y projection of the fit result drawn into it.
    *
    * The canvas needs to be given a name as identifier as well as height and width.
    * By default the whole X range is used to produce the projections, however you can also specify a range.
    *
    * The layout of the canvas follows the default drawing style of the 1D fitter.
    */
  TCanvas* CreateResultDrawingYprojection(std::string canvasname, int width, int height, double Xlow = 0, double Xhigh = 0, int BinRemovedNumber_FromRight = 0, int RebinNumber_Y = 1) const;

  /** Perform Fit
    *
    * Calls the fit appropriate fit function. You can choose with the parameter i between chi-square (i==0) or likelihood fit (i==1)(default).
    *
    */
  virtual double Fit(int i = 1, double lowBound = 0, double highBound = 0);

protected:
  /** Tests the compatibility of histograms
    *
    * Checks if the data and template histograms are compatible. They need the same number of bins, binwidth, and bin edges.
    */
  virtual bool TestHistogramCompatibility() const;

private:
  const TH2* Cast(const TH1*) const;
};

} // namespace MYUtilities

#endif
