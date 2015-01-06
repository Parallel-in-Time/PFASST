#include "particle_util.hpp"

#include <algorithm>
#include <cassert>
#include <type_traits>
using namespace std;


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      template<typename precision>
      static vector<vector<precision>> cloud_component_factory(const size_t num_particles, const size_t dim)
      {
        vector<vector<precision>> out(num_particles);
        for (size_t p = 0; p < out.size(); ++p) {
          out[p].resize(dim);
        }
        return out;
      }


      template<typename precision>
      static vector<precision> cross_prod(const vector<precision>& first,
                                                     const vector<precision>& second)
      {
        assert(first.size() == 3 && first.size() == second.size());
        vector<precision> result(3);
        result[0] = first[1] * second[2] - first[2] * second[1];
        result[1] = first[2] * second[0] - first[0] * second[2];
        result[2] = first[0] * second[1] - first[1] * second[0];
        return result;
      }

      template<typename precision>
      static vector<vector<precision>> cross_prod(const vector<vector<precision>>& first,
                                                          const vector<vector<precision>>& second)
      {
        assert(first.size() == second.size());
        vector<vector<precision>> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] = cross_prod(first[i], second[i]);
        }
        return dest;
      }

      template<typename precision>
      static vector<vector<precision>> cross_prod(const vector<vector<precision>>& first,
                                                          const vector<precision>& second)
      {
        vector<vector<precision>> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] = cross_prod(first[i], second);
        }
        return dest;
      }


      template<typename precision>
      static vector<precision> cmp_wise_mul(const vector<precision>& first,
                                                       const vector<precision>& second)
      {
        vector<precision> out(first);
        assert(out.size() == first.size());
        for (size_t i = 0; i < out.size(); ++i) {
          out[i] = first[i] * second[i];
        }
        return out;
      }

      template<typename precision>
      static vector<precision> cmp_wise_div(const vector<precision>& first,
                                                       const vector<precision>& second)
      {
        vector<precision> out(first);
        assert(out.size() == first.size());
        for (size_t i = 0; i < out.size(); ++i) {
          out[i] = first[i] / second[i];
        }
        return out;
      }


      template<typename precision>
      static precision max(const vector<precision>& data)
      {
        precision _max = precision(0.0);
        for (auto elem : data)
        {
          _max = std::max(_max, elem);
        }
        return _max;
      }


      template<typename precision>
      static precision norm0(const vector<precision>& data)
      {
        precision norm = precision(0.0);
        for (auto elem : data) {
          norm += elem * elem;
        }
        return sqrt(norm);
      }

      template<typename precision>
      static vector<precision> norm0(const vector<vector<precision>>& data)
      {
        vector<precision> norm(data.size(), precision(0.0));
        for (size_t p = 0; p < data.size(); ++p) {
          norm[p] = norm0(data[p]);
        }
        return norm;
      }


      template<typename precision>
      static void zero(vector<precision>& data)
      {
        for (auto iter = data.begin(); iter != data.end(); ++iter) {
          *iter = precision(0.0);
        }
      }

      template<typename precision>
      static void zero(vector<vector<precision>>& data)
      {
        for (auto iter = data.begin(); iter != data.end(); ++iter) {
          zero(*iter);
        }
      }

      template<typename precision>
      static void zero(shared_ptr<vector<vector<precision>>>& data)
      {
        zero(*data.get());
      }


      ////
      // OPERATORS: ADD
      template<typename precision>
      vector<precision> operator+(const vector<precision>& first,
                                             const vector<precision>& second)
      {
        vector<precision> dest(first);
        assert(dest.size() == first.size());
        dest += second;
        return dest;
      }

      template<typename precision>
      vector<vector<precision>> operator+(const vector<vector<precision>>& first,
                                                  const vector<vector<precision>>& second)
      {
        vector<vector<precision>> dest(first);
        assert(dest.size() == first.size());
        for (size_t index = 0; index < dest.size(); ++index) {
          dest[index] += second[index];
        }
        return dest;
      }

      template<typename precision>
      vector<vector<precision>> operator+(const vector<vector<precision>>& first,
                                                  const vector<precision>& second)
      {
        vector<vector<precision>> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] += second;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<precision> operator+(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] += value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>> operator+(const vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<vector<precision>> dest(vec);
        assert(dest->size() == vec->size());
        for (size_t index = 0; index < dest.size(); ++index) {
          dest[index] += value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<precision> operator+(const ValueT& value, const vector<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec + value;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>> operator+(const ValueT& value, const vector<vector<precision>>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec + value;
      }
      // OPERATORS: ADD
      ////

      ////
      // OPERATORS: INPLACE-ADD
      template<typename precision>
      vector<precision>& operator+=(vector<precision>& first,
                                               const vector<precision>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] += second[i];
        }
        return first;
      }

      template<typename precision>
      vector<vector<precision>>& operator+=(vector<vector<precision>>& first,
                                                    const vector<vector<precision>>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] += second[i];
        }
        return first;
      }

      template<typename precision>
      vector<vector<precision>>& operator+=(vector<vector<precision>>& first,
                                                    const vector<precision>& second)
      {
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] += second;
        }
        return first;
      }

      template<typename precision, typename ValueT>
      vector<precision>& operator+=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] += value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator+=(vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] += value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-ADD
      ////

      ////
      // OPERATORS: MINUS
      template<typename precision>
      vector<precision> operator-(const vector<precision>& first,
                                             const vector<precision>& second)
      {
        vector<precision> dest(first);
        assert(dest.size() == first.size());
        dest -= second;
        return dest;
      }

      template<typename precision>
      vector<vector<precision>> operator-(const vector<vector<precision>>& first,
                                                  const vector<vector<precision>>& second)
      {
        vector<vector<precision>> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= second[i];
        }
        return dest;
      }
      template<typename precision>
      vector<vector<precision>> operator-(const vector<vector<precision>>& first,
                                                  const vector<precision>& second)
      {
        vector<vector<precision>> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= second;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<precision> operator-(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>> operator-(const vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<vector<precision>> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<precision> operator-(const ValueT& value, const vector<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec - value;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>> operator-(const ValueT& value, const vector<vector<precision>>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec - value;
      }
      // OPERATORS: MINUS
      ////

      ////
      // OPERATORS: INPLACE-MINUS
      template<typename precision>
      vector<precision>& operator-=(vector<precision>& first, const vector<precision>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] -= second[i];
        }
        return first;
      }

      template<typename precision>
      vector<vector<precision>>& operator-=(vector<vector<precision>>& first,
                                                    const vector<vector<precision>>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] -= second[i];
        }
        return first;
      }

      template<typename precision>
      vector<vector<precision>>& operator-=(vector<vector<precision>>& first,
                                                    const vector<precision>& second)
      {
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] -= second;
        }
        return first;
      }

      template<typename precision, typename ValueT>
      vector<precision>& operator-=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] -= value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator-=(vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] -= value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-MINUS
      ////

      ////
      // OPERATORS: MUL
      template<typename precision, typename ValueT>
      vector<precision> operator*(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] *= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>> operator*(const vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<vector<precision>> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] *= value;
        }
        return dest;
      }

      template<typename precision>
      vector<vector<precision>> operator*(const vector<vector<precision>>& vec,
                                                  const vector<precision>& values)
      {
        assert(vec.size() == values.size());
        vector<vector<precision>> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] *= values[i];
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<precision> operator*(const ValueT& value, const vector<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec * value;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>> operator*(const ValueT& value, const vector<vector<precision>>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec * value;
      }
      // OPERATORS: MUL
      ////

      ////
      // OPERATORS: INPLACE-MUL
      template<typename precision, typename ValueT>
      vector<precision>& operator*=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] *= value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator*=(vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] *= value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-MUL
      ////

      ////
      // OPERATORS: DIV
      template<typename precision, typename ValueT>
      vector<precision>  operator/(const vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] /= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator/(const vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        vector<vector<precision>> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] /= value;
        }
        return dest;
      }

      template<typename precision>
      vector<vector<precision>> operator/(const vector<vector<precision>>& vec,
                                                  const vector<precision>& values)
      {
        assert(vec.size() == values.size());
        vector<vector<precision>> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] /= values[i];
        }
        return dest;
      }
      // OPERATORS: DIV
      ////

      ////
      // OPERATORS: INPLACE-DIV
      template<typename precision, typename ValueT>
      vector<precision>& operator/=(vector<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] /= value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator/=(vector<vector<precision>>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] /= value;
        }
        return vec;
      }
      // OPERATORS: INPLACE-DIV
      ////
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst
