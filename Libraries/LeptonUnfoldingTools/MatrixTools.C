#include "MatrixTools.hh"

#include <TH1.h>
#include <TH2.h>
#include <TClass.h>
#include <TMath.h>

#include <BinningDefinition.hh>
#include <BinningTools.hh>
#include <Quantity.hh>
#include <Statistics.hh>

#define INFO_OUT_TAG "MatrixTools"
#include <debugging.hh>

namespace MatrixTools {

TH2* MakeMatrixStrippedOfOutermostBins(const TH2* h) {

  auto newBinsX = Binning::Tools::FromTAxis(h->GetXaxis()).Bins();
  newBinsX.erase(newBinsX.begin());
  newBinsX.pop_back();
  auto newBinsY = Binning::Tools::FromTAxis(h->GetYaxis()).Bins();
  newBinsY.erase(newBinsY.begin());
  newBinsY.pop_back();

  TH2* hStripped = nullptr;
  if (h->IsA()->InheritsFrom("TH2F"))
    hStripped = Make<TH2F>(Form("%s_stripped", h->GetName()), h->GetTitle(), Binning::Tools::FromVector(newBinsX), Binning::Tools::FromVector(newBinsY));
  else
    hStripped = Make<TH2D>(Form("%s_stripped", h->GetName()), h->GetTitle(), Binning::Tools::FromVector(newBinsX), Binning::Tools::FromVector(newBinsY));

  hStripped->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  hStripped->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

  for (int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
    for (int jbin = 1; jbin <= h->GetNbinsY(); ++jbin) {
      hStripped->SetBinContent(ibin-1, jbin-1, h->GetBinContent(ibin, jbin));
        hStripped->SetBinError(ibin-1, jbin-1, h->GetBinError(ibin, jbin));
    }
  }

  return hStripped;
}

TH2* Transpose(TH2* h) {

  if (!h)
    return h;

  if (Binning::Tools::FromTAxis(h->GetXaxis()) != Binning::Tools::FromTAxis(h->GetYaxis()))
    FATAL_OUT << "Cannot transpose histogram " << h->GetName() << " because x-axis and y-axis are different." << std::endl;

  std::string xtitle = h->GetXaxis()->GetTitle();
  h->GetXaxis()->SetTitle(h->GetYaxis()->GetTitle());
  h->GetYaxis()->SetTitle(xtitle.c_str());

  for (int xbin = 0; xbin <= h->GetNbinsX(); ++xbin) {
    for (int ybin = xbin + 1; ybin <= h->GetNbinsY() + 1; ++ybin){

      double upperBinContent = h->GetBinContent(xbin, ybin);
      double upperBinError = h->GetBinError(xbin, ybin);
      double lowerBinContent = h->GetBinContent(ybin, xbin);
      double lowerBinError = h->GetBinError(ybin, xbin);

      h->SetBinContent(xbin, ybin, lowerBinContent);
      h->SetBinError(xbin, ybin, lowerBinError);
      h->SetBinContent(ybin, xbin, upperBinContent);
      h->SetBinError(ybin, xbin, upperBinError);
    }
  }

  return h;
}

TH1* MakeStrippedHistogram1D(const TH1* h, const Binning::Definition& subbinning, UnderOverflowHandling option) {

  TH1* hStripped = nullptr;

  if (Binning::Tools::SubRange(Binning::Tools::FromTAxis(h->GetXaxis()), subbinning.Min(), subbinning.Max(), true) != subbinning)
    FATAL_OUT << "New binning is not a sub-range of original histogram axis!" << std::endl;

  if (h->IsA()->InheritsFrom("TH1F"))
    hStripped = Make<TH1F>(Form("%s_stripped", h->GetName()), h->GetTitle(), subbinning);
  else
    hStripped = Make<TH1D>(Form("%s_stripped", h->GetName()), h->GetTitle(), subbinning);

  hStripped->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  hStripped->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

  std::vector<double> underflows;
  std::vector<double> overflows;
  std::vector<double> ufe2s;
  std::vector<double> ofe2s;

  for (int ibin = 0; ibin <= h->GetNbinsX() + 1; ++ibin) {

    double x = h->GetXaxis()->GetBinCenter(ibin);
    double c = h->GetBinContent(ibin);
    double e = h->GetBinError(ibin);

    if (x < subbinning.Min()) {
      underflows.push_back(c);
      ufe2s.push_back(e*e);
    }
    else if (subbinning.IsInRange(x)) {
      auto newbin = subbinning.FindBin(x);
      hStripped->SetBinContent(newbin, c);
      hStripped->SetBinError(newbin, e);
    }
    else if (x > subbinning.Max()) {
      overflows.push_back(c);
      ofe2s.push_back(e*e);
    }
    else
      assert(false);
  }

  double underflow = 0.0;
  double overflow = 0.0;
  double ufe2 = 0.0;
  double ofe2 = 0.0;

  if (option == UnderOverflowHandling::Sum) {
    for (double x : underflows)
      underflow += x;
    for (double x : overflows)
      overflow += x;
    for (double x : ufe2s)
      ufe2 += x;
    for (double x : ofe2s)
      ofe2 += x;
  }
  else if (option == UnderOverflowHandling::WeightedAverage) {
    std::vector<double> ufes;
    std::vector<double> ofes;
    for (double x : ufe2s)
      ufes.push_back(std::sqrt(x));
    for (double x : ofe2s)
      ofes.push_back(std::sqrt(x));
    Utilities::Quantity uflow = Utilities::WeightedMean(underflows, ufes);
    Utilities::Quantity oflow = Utilities::WeightedMean(overflows, ofes);
    underflow = uflow.value;
    ufe2 = std::pow(uflow.uncertainty, 2);
    overflow = oflow.value;
    ofe2 = std::pow(oflow.uncertainty, 2);
  }
  else if (option == UnderOverflowHandling::ClosestBin) {
    if (!underflows.empty()) {
      underflow = underflows.back();
      ufe2 = ufe2s.back();
    }
    if (!overflows.empty()) {
      overflow = overflows.front();
      ofe2 = ofe2s.front();
    }
  }

  hStripped->SetBinContent(0, underflow);
  hStripped->SetBinError(0, TMath::Sqrt(ufe2));
  hStripped->SetBinContent(hStripped->GetNbinsX()+1, overflow);
  hStripped->SetBinContent(hStripped->GetNbinsX()+1, TMath::Sqrt(ofe2));

  return hStripped;
}

TH2* MakeStrippedHistogram2D(const TH2* h, const Binning::Definition& subbinning, bool setUnderOverflows) {

  TH2* hStripped = nullptr;

  if (Binning::Tools::SubRange(Binning::Tools::FromTAxis(h->GetXaxis()), subbinning.Min(), subbinning.Max(), true) != subbinning)
    FATAL_OUT << "New binning is not a sub-range of original histogram x-axis!" << std::endl;
  if (Binning::Tools::SubRange(Binning::Tools::FromTAxis(h->GetYaxis()), subbinning.Min(), subbinning.Max(), true) != subbinning)
    FATAL_OUT << "New binning is not a sub-range of original histogram y-axis!" << std::endl;

  if (h->IsA()->InheritsFrom("TH2F"))
    hStripped = Make<TH2F>(Form("%s_stripped", h->GetName()), h->GetTitle(), subbinning, subbinning);
  else
    hStripped = Make<TH2D>(Form("%s_stripped", h->GetName()), h->GetTitle(), subbinning, subbinning);

  hStripped->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  hStripped->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());

  int firstnewbinX = h->GetXaxis()->FindFixBin(subbinning.Value(1));
  int lastnewbinX = h->GetXaxis()->FindFixBin(subbinning.Value(subbinning.NumberOfBins()));
  int firstnewbinY = h->GetYaxis()->FindFixBin(subbinning.Value(1));
  int lastnewbinY = h->GetYaxis()->FindFixBin(subbinning.Value(subbinning.NumberOfBins()));

  for (int xbin = firstnewbinX; xbin <= lastnewbinX; ++xbin) {
    for (int ybin = firstnewbinY; ybin <= lastnewbinY; ++ybin) {

      double x = h->GetXaxis()->GetBinCenter(xbin);
      double y = h->GetYaxis()->GetBinCenter(ybin);
      double c = h->GetBinContent(xbin, ybin);
      double e = h->GetBinError(xbin, ybin);

      auto newbinX = subbinning.FindBin(x);
      auto newbinY = subbinning.FindBin(y);
      hStripped->SetBinContent(newbinX, newbinY, c);
      hStripped->SetBinError(newbinX, newbinY, e);
    }
  }

  if (setUnderOverflows) {

    for (int ybin = firstnewbinY; ybin <= lastnewbinY; ++ybin) {

      double y = h->GetYaxis()->GetBinCenter(ybin);
      auto newbinY = subbinning.FindBin(y);

      // A.1) fill underflows along x
      double underflow = 0.0;
      double ufe2 = 0.0;
      for (int xbin = 0; xbin < firstnewbinX; ++xbin) {
        underflow += h->GetBinContent(xbin, ybin);
        ufe2 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
      hStripped->SetBinContent(0, newbinY, underflow);
      hStripped->SetBinError(0, newbinY, TMath::Sqrt(ufe2));

      // A.2) fill overflows along x
      double overflow = 0.0;
      double ofe2 = 0.0;
      for (int xbin = lastnewbinX+1; xbin <= h->GetNbinsX()+1; ++xbin) {
        overflow += h->GetBinContent(xbin, ybin);
        ofe2 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
      hStripped->SetBinContent(hStripped->GetNbinsX()+1, newbinY, overflow);
      hStripped->SetBinError(hStripped->GetNbinsX()+1, newbinY, TMath::Sqrt(ofe2));
    }

    for (int xbin = firstnewbinX; xbin <= lastnewbinX; ++xbin) {

      double x = h->GetXaxis()->GetBinCenter(xbin);
      auto newbinX = subbinning.FindBin(x);

      // B.1) fill underflows along y
      double underflow = 0.0;
      double ufe2 = 0.0;
      for (int ybin = 0; ybin < firstnewbinY; ++ybin) {
        underflow += h->GetBinContent(xbin, ybin);
        ufe2 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
      hStripped->SetBinContent(newbinX, 0, underflow);
      hStripped->SetBinError(newbinX, 0, TMath::Sqrt(ufe2));

      // B.2) fill overflows along y
      double overflow = 0.0;
      double ofe2 = 0.0;
      for (int ybin = lastnewbinY+1; ybin <= h->GetNbinsY()+1; ++ybin) {
        overflow += h->GetBinContent(xbin, ybin);
        ofe2 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
      hStripped->SetBinContent(newbinX, hStripped->GetNbinsY()+1, overflow);
      hStripped->SetBinError(newbinX, hStripped->GetNbinsY()+1, TMath::Sqrt(ofe2));
    }

    // C) fill edge cases

    // C.1) bottom-left
    double c1 = 0.0;
    double e12 = 0.0;
    for (int xbin = 0; xbin < firstnewbinX; ++xbin) {
      for (int ybin = 0; ybin < firstnewbinY; ++ybin) {
        c1 += h->GetBinContent(xbin, ybin);
        e12 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
    }
    hStripped->SetBinContent(0, 0, c1);
    hStripped->SetBinError(0, 0, TMath::Sqrt(e12));

    // C.2) top-left
    double c2 = 0.0;
    double e22 = 0.0;
    for (int xbin = 0; xbin < firstnewbinX; ++xbin) {
      for (int ybin = lastnewbinY+1; ybin <= h->GetNbinsY()+1; ++ybin) {
        c2 += h->GetBinContent(xbin, ybin);
        e22 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
    }
    hStripped->SetBinContent(0, hStripped->GetNbinsY()+1, c2);
    hStripped->SetBinError(0, hStripped->GetNbinsY()+1, TMath::Sqrt(e22));

    // C.3) bottom-right
    double c3 = 0.0;
    double e32 = 0.0;
    for (int xbin = lastnewbinX+1; xbin <= h->GetNbinsX()+1; ++xbin) {
      for (int ybin = 0; ybin < firstnewbinY; ++ybin) {
        c3 += h->GetBinContent(xbin, ybin);
        e32 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
    }
    hStripped->SetBinContent(hStripped->GetNbinsX()+1, 0, c3);
    hStripped->SetBinError(hStripped->GetNbinsX()+1, 0, TMath::Sqrt(e32));

    // C.4) top-right
    double c4 = 0.0;
    double e42 = 0.0;
    for (int xbin = lastnewbinX+1; xbin <= h->GetNbinsX()+1; ++xbin) {
      for (int ybin = lastnewbinY+1; ybin <= h->GetNbinsY()+1; ++ybin) {
        c4 += h->GetBinContent(xbin, ybin);
        e42 += std::pow(h->GetBinError(xbin, ybin), 2);
      }
    }
    hStripped->SetBinContent(hStripped->GetNbinsX()+1, hStripped->GetNbinsY()+1, c4);
    hStripped->SetBinError(hStripped->GetNbinsX()+1, hStripped->GetNbinsY()+1, TMath::Sqrt(e42));
  }

  return hStripped;
}

}




