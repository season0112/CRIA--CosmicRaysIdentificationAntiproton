#ifndef Quantity_hh
#define Quantity_hh

#include <ostream>

namespace MYUtilities {

/** A simple structure for a physical value with an uncertainty associated to it. */
struct Quantity {
  Quantity(double Value = 0.0, double Uncertainty = 0.0)
    : value(Value)
    , uncertainty(Uncertainty) {
  }

  /// value
  double value;
  /// associated uncertainty
  double uncertainty;
};

std::ostream& operator<<(std::ostream& out, const Quantity& q);

} // namespace MYUtilities

#endif
