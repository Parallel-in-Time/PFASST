#ifndef _PFASST__UTIL_HPP_
#define _PFASST__UTIL_HPP_

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
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

  template<typename T>
  static string join(const vector<T>& vec, const string& sep)
  {
    stringstream out;
    out << "[";
    for (size_t i = 0; i < vec.size() - 1; ++i) {
      out << vec[i] << sep;
    }
    out << vec.back();
    out << "]";
    return out.str();
  }
}  // ::pfasst

#endif  // _PFASST__UTIL_HPP_
