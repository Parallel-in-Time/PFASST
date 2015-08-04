#ifndef _PFASST__UTIL_HPP_
#define _PFASST__UTIL_HPP_

#include <algorithm>
#include <limits>
using namespace std;

#include "pfasst/logging.hpp"


namespace pfasst
{
  template<typename precision>
  static bool almost_equal(const precision& a, const precision& b,
                           const int digits = numeric_limits<precision>::digits)
  {
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return abs(a - b) < numeric_limits<precision>::epsilon() * abs(a + b) * digits
    // unless the result is subnormal
           || abs(a - b) < numeric_limits<precision>::min();
  }

  template<typename precision>
  static bool almost_zero(const precision& a)
  {
    return abs(a) < numeric_limits<precision>::epsilon();
  }
}  // ::pfasst

#endif  // _PFASST__UTIL_HPP_
