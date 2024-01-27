#ifndef Matrix_Tools_hh_
#define Matrix_Tools_hh_

namespace Binning {
class Definition;
}

class TH1;
class TH2;

namespace MatrixTools {

enum class UnderOverflowHandling { Zero, Sum, WeightedAverage, ClosestBin };

TH2* MakeMatrixStrippedOfOutermostBins(const TH2* h);

TH1* MakeStrippedHistogram1D(const TH1* h, const Binning::Definition& subbinning, UnderOverflowHandling option);
TH2* MakeStrippedHistogram2D(const TH2* h, const Binning::Definition& subbinning, bool setUnderOverflows);

/** Transpose a quadratic matrix in place. */
TH2* Transpose(TH2* h);

}

#endif
