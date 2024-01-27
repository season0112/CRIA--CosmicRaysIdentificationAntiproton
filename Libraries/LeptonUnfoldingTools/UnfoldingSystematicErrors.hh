#ifndef UnfoldingSystematicErrors_hh_
#define UnfoldingSystematicErrors_hh_

#include <vector>

class TGraph;

class UnfoldingSystematicErrors {

public:
  UnfoldingSystematicErrors() { }
  ~UnfoldingSystematicErrors() { }

  /** Relative systematic uncertainty on unfolded flux: contribution from unfolding bias seen in toy MC. */
  static TGraph* MakeUncertaintyFromToymcBias();

  /** Relative systematic uncertainty on unfolded flux: contribution from uncertainty in migration matrix. */
  static TGraph* MakeUncertaintyFromTestbeam();

  /** Total relative systematic uncertainty on unfolded flux.
    *
    * The uncertainty is calculated as the quadratic sum of the individual contributions:
    *   - unfolding bias evaluated in toy MC
    *   - impact of uncertainty in migration matrix, estimated from testbeam, on the fluxes
    *
    * \param[in] xPoints Vector of points where the uncertainty is to be evaluated.
    */
  static TGraph* MakeTotalSystErrorFromUnfolding(const std::vector<double>& xPoints);
};

#endif
