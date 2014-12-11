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
      static ParticleCloudComponent<precision> cloud_component_factory(const size_t num_particles, const size_t dim)
      {
        ParticleCloudComponent<precision> out(num_particles);
        for (size_t p = 0; p < out.size(); ++p) {
          out[p].resize(dim);
        }
        return out;
      }


      template<typename precision>
      static ParticleComponent<precision> cross_prod(const ParticleComponent<precision>& first,
                                                     const ParticleComponent<precision>& second)
      {
        assert(first.size() == 3 && first.size() == second.size());
        ParticleComponent<precision> result(3);
        result[0] = first[1] * second[2] - first[2] * second[1];
        result[1] = first[2] * second[0] - first[0] * second[2];
        result[2] = first[0] * second[1] - first[1] * second[0];
        return result;
      }

      template<typename precision>
      static ParticleCloudComponent<precision> cross_prod(const ParticleCloudComponent<precision>& first,
                                                          const ParticleCloudComponent<precision>& second)
      {
        assert(first.size() == second.size());
        ParticleCloudComponent<precision> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] = cross_prod(first[i], second[i]);
        }
        return dest;
      }

      template<typename precision>
      static ParticleCloudComponent<precision> cross_prod(const ParticleCloudComponent<precision>& first,
                                                          const ParticleComponent<precision>& second)
      {
        ParticleCloudComponent<precision> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] = cross_prod(first[i], second);
        }
        return dest;
      }


      template<typename precision>
      static ParticleComponent<precision> cmp_wise_mul(const ParticleComponent<precision>& first,
                                                       const ParticleComponent<precision>& second)
      {
        ParticleComponent<precision> out(first);
        assert(out.size() == first.size());
        for (size_t i = 0; i < out.size(); ++i) {
          out[i] = first[i] * second[i];
        }
        return out;
      }

      template<typename precision>
      static ParticleComponent<precision> cmp_wise_div(const ParticleComponent<precision>& first,
                                                       const ParticleComponent<precision>& second)
      {
        ParticleComponent<precision> out(first);
        assert(out.size() == first.size());
        for (size_t i = 0; i < out.size(); ++i) {
          out[i] = first[i] / second[i];
        }
        return out;
      }


      template<typename precision>
      static precision distance(const Particle<precision>& first, const Particle<precision>& second)
      {
        assert(first.DIM() == second.DIM());
        ParticleComponent<precision> diff = first.pos() - second.pos();
        precision dist = precision(0.0);
        for (size_t i = 0; i < first.DIM(); ++i) {
          dist += diff[i] * diff[i];
        }
        return sqrt(dist);
      }

      template<typename precision>
      static precision distance(const shared_ptr<Particle<precision>> first, const shared_ptr<Particle<precision>> second)
      {
        return distance(*(first.get()), *(second.get()));
      }

      template<typename precision>
      static vector<precision> distance_to_reference(const ParticleCloud<precision>& cloud,
                                                     const Particle<precision>&      reference)
      {
        vector<precision> distances(cloud.size());
        for (size_t i = 0; i < distances.size(); ++i) {
          distances[i] = distance(*cloud[i], reference);
        }
        return distances;
      }

      template<typename precision>
      static vector<precision> distance_to_reference(const shared_ptr<ParticleCloud<precision>>& cloud,
                                                     const shared_ptr<Particle<precision>>&      reference)
      {
        return distance_to_reference(*cloud, *reference);
      }


      template<typename precision>
      static precision norm0(const ParticleComponent<precision>& data)
      {
        precision norm = precision(0.0);
        for (auto elem : data) {
          norm += elem * elem;
        }
        return sqrt(norm);
      }

      template<typename precision>
      static vector<precision> norm0(const ParticleCloudComponent<precision>& data)
      {
        vector<precision> norm(data.size(), precision(0.0));
        for (size_t p = 0; p < data.size(); ++p) {
          norm[p] = norm0(data[p]);
        }
        return norm;
      }


      template<typename precision>
      static void zero(ParticleCloudComponent<precision>& data)
      {
        for (auto p : data) {
          fill(p.begin(), p.end(), precision(0.0));
        }
      }

      template<typename precision>
      static void zero(shared_ptr<ParticleCloudComponent<precision>>& data)
      {
        zero(*data.get());
      }


      ////
      // OPERATORS: ADD
      template<typename precision>
      ParticleComponent<precision> operator+(const ParticleComponent<precision>& first,
                                             const ParticleComponent<precision>& second)
      {
        ParticleComponent<precision> dest(first);
        assert(dest.size() == first.size());
        dest += second;
        return dest;
      }

      template<typename precision>
      ParticleCloudComponent<precision> operator+(const ParticleCloudComponent<precision>& first,
                                                  const ParticleCloudComponent<precision>& second)
      {
        ParticleCloudComponent<precision> dest(first);
        assert(dest.size() == first.size());
        for (size_t index = 0; index < dest.size(); ++index) {
          dest[index] += second[index];
        }
        return dest;
      }

      template<typename precision>
      ParticleCloudComponent<precision> operator+(const ParticleCloudComponent<precision>& first,
                                                  const ParticleComponent<precision>& second)
      {
        ParticleCloudComponent<precision> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] += second;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleComponent<precision> operator+(const ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] += value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision> operator+(const ParticleCloudComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleCloudComponent<precision> dest(vec);
        assert(dest->size() == vec->size());
        for (size_t index = 0; index < dest.size(); ++index) {
          dest[index] += value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleComponent<precision> operator+(const ValueT& value, const ParticleComponent<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec + value;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision> operator+(const ValueT& value, const ParticleCloudComponent<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec + value;
      }
      // OPERATORS: ADD
      ////

      ////
      // OPERATORS: INPLACE-ADD
      template<typename precision>
      ParticleComponent<precision>& operator+=(ParticleComponent<precision>& first,
                                               const ParticleComponent<precision>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] += second[i];
        }
        return first;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& operator+=(ParticleCloudComponent<precision>& first,
                                                    const ParticleCloudComponent<precision>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] += second[i];
        }
        return first;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& operator+=(ParticleCloudComponent<precision>& first,
                                                    const ParticleComponent<precision>& second)
      {
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] += second;
        }
        return first;
      }

      template<typename precision, typename ValueT>
      ParticleComponent<precision>& operator+=(ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] += value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision>& operator+=(ParticleCloudComponent<precision>& vec, const ValueT& value)
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
      ParticleComponent<precision> operator-(const ParticleComponent<precision>& first,
                                             const ParticleComponent<precision>& second)
      {
        ParticleComponent<precision> dest(first);
        assert(dest.size() == first.size());
        dest -= second;
        return dest;
      }

      template<typename precision>
      ParticleCloudComponent<precision> operator-(const ParticleCloudComponent<precision>& first,
                                                  const ParticleCloudComponent<precision>& second)
      {
        ParticleCloudComponent<precision> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= second[i];
        }
        return dest;
      }
      template<typename precision>
      ParticleCloudComponent<precision> operator-(const ParticleCloudComponent<precision>& first,
                                                  const ParticleComponent<precision>& second)
      {
        ParticleCloudComponent<precision> dest(first);
        assert(dest.size() == first.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= second;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleComponent<precision> operator-(const ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision> operator-(const ParticleCloudComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleCloudComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] -= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleComponent<precision> operator-(const ValueT& value, const ParticleComponent<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec - value;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision> operator-(const ValueT& value, const ParticleCloudComponent<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec - value;
      }
      // OPERATORS: MINUS
      ////

      ////
      // OPERATORS: INPLACE-MINUS
      template<typename precision>
      ParticleComponent<precision>& operator-=(ParticleComponent<precision>& first, const ParticleComponent<precision>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] -= second[i];
        }
        return first;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& operator-=(ParticleCloudComponent<precision>& first,
                                                    const ParticleCloudComponent<precision>& second)
      {
        assert(first.size() == second.size());
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] -= second[i];
        }
        return first;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& operator-=(ParticleCloudComponent<precision>& first,
                                                    const ParticleComponent<precision>& second)
      {
        for (size_t i = 0; i < first.size(); ++i) {
          first[i] -= second;
        }
        return first;
      }

      template<typename precision, typename ValueT>
      ParticleComponent<precision>& operator-=(ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] -= value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision>& operator-=(ParticleCloudComponent<precision>& vec, const ValueT& value)
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
      ParticleComponent<precision> operator*(const ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] *= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision> operator*(const ParticleCloudComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleCloudComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] *= value;
        }
        return dest;
      }

      template<typename precision>
      ParticleCloudComponent<precision> operator*(const ParticleCloudComponent<precision>& vec,
                                                  const AttributeValues<precision>& values)
      {
        assert(vec.size() == values.size());
        ParticleCloudComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] *= values[i];
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleComponent<precision> operator*(const ValueT& value, const ParticleComponent<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec * value;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision> operator*(const ValueT& value, const ParticleCloudComponent<precision>& vec)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        return vec * value;
      }
      // OPERATORS: MUL
      ////

      ////
      // OPERATORS: INPLACE-MUL
      template<typename precision, typename ValueT>
      ParticleComponent<precision>& operator*=(ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] *= value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision>& operator*=(ParticleCloudComponent<precision>& vec, const ValueT& value)
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
      ParticleComponent<precision>  operator/(const ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] /= value;
        }
        return dest;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision>  operator/(const ParticleCloudComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        ParticleCloudComponent<precision> dest(vec);
        assert(dest.size() == vec.size());
        for (size_t i = 0; i < dest.size(); ++i) {
          dest[i] /= value;
        }
        return dest;
      }

      template<typename precision>
      ParticleCloudComponent<precision> operator/(const ParticleCloudComponent<precision>& vec,
                                                  const AttributeValues<precision>& values)
      {
        assert(vec.size() == values.size());
        ParticleCloudComponent<precision> dest(vec);
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
      ParticleComponent<precision>& operator/=(ParticleComponent<precision>& vec, const ValueT& value)
      {
        static_assert(is_arithmetic<ValueT>::value, "");
        for (size_t i = 0; i < vec.size(); ++i) {
          vec[i] /= value;
        }
        return vec;
      }

      template<typename precision, typename ValueT>
      ParticleCloudComponent<precision>& operator/=(ParticleCloudComponent<precision>& vec, const ValueT& value)
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
