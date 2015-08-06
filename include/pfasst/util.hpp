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
  /**
   * Compares absolute difference against zero.
   *
   * Takes the absolute value of the difference between @p a and @p b and compares it to the machine precision scaled
   * to the magnitude of @p a and @p b up to the desired precision in units in the last place.
   *
   * @tparam    precision value type of @p a and @p b
   * @param[in] a
   * @param[in] b
   * @param[in] digits    number of digits to compare
   * @ingroup Utilities
   */
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

  /**
   * Compares against smalles machine epsilon.
   *
   * Takes the absolute value of @p a and compares it against `std::numeric_limits<precision>::epsilon()`.
   *
   * @tparam    precision value type of @p a, e.g. `double`
   * @param[in] a         value to compare
   * @returns   `true` if @p a is closer to zero than the machine precision, `false` otherwise.
   *
   * @ingroup Utilities
   */
  template<typename precision>
  static bool almost_zero(const precision& a)
  {
    return abs(a) < numeric_limits<precision>::epsilon();
  }

  /**
   * A string representation of `std::vector`.
   *
   * Uses `operator<<` on each element of @p vec and seperates each element with @p sep .
   *
   * @tparam    T   type of the vector's elements
   * @param[in] vec vector to print
   * @param[in] sep seperator placed between each element of @p vec
   * @returns string representation of @p vec
   *
   * **Example**:
   * @code{.cpp}
   * std::vector<double> vec{0.0, 0.5, 0.1};
   * join(vec, ", ");  // ==> "[0.0, 0.5, 0.1]"
   * @endcode
   *
   * @ingroup Utilities
   */
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
