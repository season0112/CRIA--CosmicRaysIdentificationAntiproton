#include "Quantity.hh"

namespace MYUtilities {

std::ostream& operator<<(std::ostream& out, const MYUtilities::Quantity& q) {

  out << q.value << "+-" << q.uncertainty;
  return out;
}

} // namespace MYUtilities
