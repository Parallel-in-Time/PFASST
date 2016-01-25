#include "particle_util.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <type_traits>
using namespace std;

#include <pfasst/globals.hpp>
#include <pfasst/logging.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      template<typename precision>
      inline static vector<precision> cloud_component_factory(const size_t num_particles, const size_t dim)
      {
        vector<precision> out(num_particles * dim, precision(0.0));
        return out;
      }


      template<typename precision>
      inline static void zero(vector<precision>& data)
      {
        std::fill(data.begin(), data.end(), precision(0.0));
      }

      template<typename precision>
      inline static void zero(shared_ptr<vector<precision>>& data)
      {
        zero(*data.get());
      }


      template<typename precision>
      inline static vector<precision> cross_prod(const vector<precision>& first, const vector<precision>& second)
      {
        if (first.size() == 3 && first.size() == second.size()) {
          vector<precision> result(3, precision(0.0));
          cross_prod_1part<precision>(first.cbegin(), second.cbegin(), result.begin());
          return result;
        } else {
          return cross_prod_npart(first, second);
        }
      }

      template<typename precision>
      inline static void cross_prod_1part(typename vector<precision>::const_iterator __first,
                                          typename vector<precision>::const_iterator __second,
                                          typename vector<precision>::iterator __result)
      {
        __result[0] = __first[1] * __second[2] - __first[2] * __second[1];
        __result[1] = __first[2] * __second[0] - __first[0] * __second[2];
        __result[2] = __first[0] * __second[1] - __first[1] * __second[0];
      }

      template<typename precision>
      inline static vector<precision> cross_prod_npart(const vector<precision>& first, const vector<precision>& second)
      {
        assert(first.size() % 3 == 0 && second.size() % 3 == 0);  // make sure the particles have 3 spacial dimensions
        assert(first.size() == second.size() || second.size() == 3);
        const size_t npart = first.size() / 3;
        vector<precision> dest(first.size(), precision(0.0));
        if (first.size() == second.size()) {
          for (size_t p = 0; p < npart; ++p) {
            cross_prod_1part<precision>(first.cbegin() + (p * 3), second.cbegin() + (p * 3), dest.begin() + (p * 3));
          }
        } else if (second.size() == 3) {
          for (size_t p = 0; p < npart; ++p) {
            cross_prod_1part<precision>(first.cbegin() + (p * 3), second.cbegin(), dest.begin() + (p * 3));
          }
        }
        return dest;
      }


      template<typename precision>
      inline static vector<precision> kronecker(const vector<precision>& first,
                                                const vector<precision>& second)
      {
        const size_t npart = first.size();
        const size_t ndim = second.size();
        vector<precision> result(npart * ndim, precision(0.0));
        for (size_t p = 0; p < npart; ++p) {
          for (size_t d = 0; d < ndim; ++d) {
            result[p * ndim + d] = first[p] * second[d];
          }
        }
        return result;
      }


      template<typename precision>
      inline static vector<precision> cmp_wise_mul(const vector<precision>& first, const vector<precision>& second)
      {
        assert(first.size() == second.size());
        vector<precision> out(first.size(), precision(0.0));
        transform(first.cbegin(), first.cend(), second.cbegin(), out.begin(), std::multiplies<precision>());
        return out;
      }

      template<typename precision>
      inline static vector<precision> cmp_wise_div(const vector<precision>& first, const vector<precision>& second)
      {
        assert(first.size() == second.size());
        vector<precision> out(first.size(), precision(0.0));
        transform(first.cbegin(), first.cend(), second.cbegin(), out.begin(), std::divides<precision>());
        return out;
      }


      template<typename precision>
      inline static precision max(const vector<precision>& data)
      {
        return *max_element(data.cbegin(), data.cend());
      }

      template<typename precision>
      inline static precision max_abs(const vector<precision>& data)
      {
        return fabs(*max_element(data.cbegin(), data.cend(), [](const precision& a, const precision& b) { return fabs(a) < fabs(b); }));
      }


      template<typename precision>
      inline static precision norm_sq(const vector<precision>& data)
      {
        auto norm = norm_sq<precision>(data.cbegin(), data.cend());
        return norm;
      }

      template<typename precision>
      inline static precision norm_sq(typename vector<precision>::const_iterator __first,
                                      typename vector<precision>::const_iterator __second)
      {
        return inner_product(__first, __second, __first, precision(0.0));
      }

      template<typename precision>
      inline static vector<precision> norm_sq_npart(const vector<precision>& data, const size_t npart)
      {
        assert(data.size() % npart == 0);
        const size_t dim = data.size() / npart;
        vector<precision> norm(npart, precision(0.0));
        for (size_t p = 0; p < npart; ++p) {
          norm[p] = norm_sq<precision>(data.cbegin() + (p * dim), data.cbegin() + ((p+1) * dim));
        }
        return norm;
      }


      template<typename precision>
      inline static precision norm0(const vector<precision>& data)
      {
        return norm0<precision>(data.cbegin(), data.cend());
      }

      template<typename precision>
      inline static precision norm0(typename vector<precision>::const_iterator __first,
                                    typename vector<precision>::const_iterator __second)
      {
        return sqrt(inner_product(__first, __second, __first, precision(0.0)));
      }

      template<typename precision>
      inline static vector<precision> norm0_npart(const vector<precision>& data, const size_t npart)
      {
        assert(data.size() % npart == 0);
        const size_t dim = data.size() / npart;
        vector<precision> norm(npart, precision(0.0));
        for (size_t p = 0; p < npart; ++p) {
          norm[p] = norm0<precision>(data.cbegin() + (p * dim), data.cend() + ((p+1) * dim));
        }
        return norm;
      }


      ////
      // OPERATORS: ADD
      template<typename precision>
      inline vector<precision> operator+(const vector<precision>& first, const vector<precision>& second)
      {
        vector<precision> dest;
        if (first.size() == second.size()) {
          dest.resize(first.size());
          std::transform(first.cbegin(), first.cend(), second.cbegin(), dest.begin(), std::plus<precision>());
        } else if (first.size() % 3 == 0 && second.size() == 3) {
          dest.resize(first.size());
          for (size_t p = 0; p < first.size() / 3; ++p) {
            std::transform(first.cbegin() + (p * 3), first.cbegin() + ((p+1) * 3),
                           second.cbegin(), dest.begin() + (p * 3), std::plus<precision>());
          }
        } else {
          // try other way round
          // ATTENTION! recursion
          ML_LOG(WARNING, "Commutativity of addition primaly implemented other way round."
                          << " You should switch to avoid unneccessary function calls.");
          dest = second + first;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      inline vector<precision> operator+(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        dest += value;
        return dest;
      }

      template<typename precision, typename ValueT>
      inline vector<precision> operator+(const ValueT& value, const vector<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ML_LOG(WARNING, "Commutativity of addition primaly implemented other way round."
                        << " You should switch to avoid unneccessary function calls.");
        return vec + value;
      }
      // OPERATORS: ADD
      ////

      ////
      // OPERATORS: INPLACE-ADD
      template<typename precision>
      inline vector<precision>& operator+=(vector<precision>& first, const vector<precision>& second)
      {
        if (first.size() == second.size()) {
          transform(first.cbegin(), first.cend(), second.cbegin(), first.begin(), std::plus<precision>());
        } else if (first.size() % 3 == 0 && second.size() == 3) {
          for (size_t p = 0; p < first.size() / 3; ++p) {
            std::transform(first.cbegin() + (p * 3), first.cbegin() + ((p+1) * 3),
                           second.cbegin(), first.begin() + (p * 3), std::plus<precision>());
          }
        } else if (first.size() == 3 && second.size() % 3 == 0) {
          for (size_t p = 0; p < first.size() / 3; ++p) {
            std::transform(first.cbegin(), first.cend(), second.cbegin() + (p * 3),
                           first.begin(), std::plus<precision>());
          }
        }
        return first;
      }

      template<typename precision, typename ValueT>
      inline vector<precision>& operator+=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (auto&& elem : vec) {
          elem += value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-ADD
      ////

      ////
      // OPERATORS: MINUS
      template<typename precision>
      inline vector<precision> operator-(const vector<precision>& first, const vector<precision>& second)
      {
        vector<precision> dest;
        if (first.size() == second.size()) {
          dest.resize(first.size());
          std::transform(first.cbegin(), first.cend(), second.cbegin(), dest.begin(), std::minus<precision>());
        } else if (first.size() % 3 == 0 && second.size() == 3) {
          dest.resize(first.size());
          for (size_t p = 0; p < first.size() / 3; ++p) {
            std::transform(first.cbegin() + (p * 3), first.cbegin() + ((p+1) * 3), second.cbegin(),
                           dest.begin() + (p * 3), std::minus<precision>());
          }
        } else if (first.size() == 3 && second.size() % 3 == 0) {
          dest.resize(first.size());
          for (size_t p = 0; p < first.size() / 3; ++p) {
            std::transform(first.cbegin(), first.cend(), second.cbegin() + (p * 3),
                           dest.begin() + (p * 3), std::minus<precision>());
          }
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      inline vector<precision> operator-(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        dest -= value;
        return dest;
      }
      // OPERATORS: MINUS
      ////

      ////
      // OPERATORS: INPLACE-MINUS
      template<typename precision>
      inline vector<precision>& operator-=(vector<precision>& first, const vector<precision>& second)
      {
        if (first.size() == second.size()) {
          transform(first.cbegin(), first.cend(), second.cbegin(), first.begin(), std::minus<precision>());
        } else if (first.size() % 3 == 0 && second.size() == 3) {
          for (size_t p = 0; p < first.size() / 3; ++p) {
            std::transform(first.cbegin() + (p * 3), first.cbegin() + ((p+1) * 3),
                           second.cbegin(), first.begin() + (p * 3), std::minus<precision>());
          }
        } else if (first.size() == 3 && second.size() % 3 == 0) {
          for (size_t p = 0; p < first.size() / 3; ++p) {
            std::transform(first.cbegin(), first.cend(), second.cbegin() + (p * 3),
                           first.begin(), std::minus<precision>());
          }
        }
        return first;
      }

      template<typename precision, typename ValueT>
      inline vector<precision>& operator-=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (auto&& elem : vec) {
          elem -= value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-MINUS
      ////

      ////
      // OPERATORS: MUL
      template<typename precision, typename ValueT>
      inline vector<precision> operator*(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        dest *= value;
        return dest;
      }

      template<typename precision, typename ValueT>
      inline vector<precision> operator*(const ValueT& value, const vector<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ML_LOG(WARNING, "Commutativity of multiplication primaly implemented other way round."
                        << " You should switch to avoid unneccessary function calls.");
        return vec * value;
      }

      template<typename precision>
      inline vector<precision> operator* (const vector<precision>& vec, const vector<precision>& values)
      {
        assert(vec.size() % 3 == 0 && vec.size() / 3 == values.size());
        vector<precision> dest(vec.size(), precision(0.0));
        for (size_t p = 0; p < values.size(); ++p) {
          std::transform(vec.cbegin() + (p * 3), vec.cbegin() + ((p+1) * 3), dest.begin() + (p * 3),
                         [&](const precision& v) { return v * values[p]; });
        }
        return dest;
      }
      // OPERATORS: MUL
      ////

      ////
      // OPERATORS: INPLACE-MUL
      template<typename precision, typename ValueT>
      inline vector<precision>& operator*=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (auto&& elem : vec) {
          elem *= value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-MUL
      ////

      ////
      // OPERATORS: DIV
      template<typename precision, typename ValueT>
      inline vector<precision> operator/(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        dest /= value;
        return dest;
      }

      template<typename precision>
      inline vector<precision> operator/ (const vector<precision>& vec, const vector<precision>& values)
      {
        assert(vec.size() % 3 == 0 && vec.size() / 3 == values.size());
        vector<precision> dest(vec.size(), precision(0.0));
        for (size_t p = 0; p < values.size(); ++p) {
          std::transform(vec.cbegin() + (p * 3), vec.cbegin() + ((p+1) * 3), dest.begin() + (p * 3),
                         [&](const precision& v) { return v / values[p]; });
        }
        return dest;
      }
      // OPERATORS: DIV
      ////

      ////
      // OPERATORS: INPLACE-DIV
      template<typename precision, typename ValueT>
      inline vector<precision>& operator/=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (auto&& elem : vec) {
          elem /= value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-DIV
      ////
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst
