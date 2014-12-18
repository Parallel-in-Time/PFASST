#ifndef _EXAMPLES__BORIS__PARTICLE_UTIL_HPP_
#define _EXAMPLES__BORIS__PARTICLE_UTIL_HPP_

#include <memory>
#include <vector>
using namespace std;


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      template<typename precision>
      static vector<vector<precision>> cloud_component_factory(const size_t num_particles, const size_t dim);


      template<typename precision>
      static vector<precision>      cross_prod(const vector<precision>&      first,
                                                          const vector<precision>&      second);
      template<typename precision>
      static vector<vector<precision>> cross_prod(const vector<vector<precision>>& first,
                                                          const vector<vector<precision>>& second);
      template<typename precision>
      static vector<vector<precision>> cross_prod(const vector<vector<precision>>& first,
                                                          const vector<precision>&      second);

      template<typename precision>
      static vector<precision> cmp_wise_mul(const vector<precision>& first,
                                                       const vector<precision>& second);

      template<typename precision>
      static vector<precision> cmp_wise_div(const vector<precision>& first,
                                                       const vector<precision>& second);
      
      
      template<typename precision>
      static precision max(const vector<precision>& data);


      template<typename precision>
      static precision norm0(const vector<precision>& data);

      template<typename precision>
      static vector<precision> norm0(const vector<vector<precision>>& data);


      template<typename precision>
      static void zero(vector<precision>& data);

      template<typename precision>
      static void zero(vector<vector<precision>>& data);

      template<typename precision>
      static void zero(shared_ptr<vector<vector<precision>>>& data);


      template<typename precision>
      vector<precision>       operator+ (const vector<precision>&      first,
                                                    const vector<precision>&      second);
      template<typename precision>
      vector<vector<precision>>  operator+ (const vector<vector<precision>>& first,
                                                    const vector<vector<precision>>& second);
      template<typename precision>
      vector<vector<precision>>  operator+ (const vector<vector<precision>>& first,
                                                    const vector<precision>&      second);
      template<typename precision, typename ValueT>
      vector<precision>       operator+ (const vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator+ (const vector<vector<precision>>& vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<precision>       operator+ (const ValueT&                            value,
                                                    const vector<precision>&      vec   );
      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator+ (const ValueT&                            value,
                                                    const vector<vector<precision>>& vec   );

      template<typename precision>
      vector<precision>&      operator+=(      vector<precision>&      first,
                                                    const vector<precision>&      second);
      template<typename precision>
      vector<vector<precision>>& operator+=(      vector<vector<precision>>& first,
                                                    const vector<vector<precision>>& second);
      template<typename precision>
      vector<vector<precision>>& operator+=(      vector<vector<precision>>& first,
                                                    const vector<precision>&      second);
      template<typename precision, typename ValueT>
      vector<precision>&      operator+=(      vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator+=(      vector<vector<precision>>& vec,
                                                    const ValueT&                            value );

      template<typename precision>
      vector<precision>       operator- (const vector<precision>&      first,
                                                    const vector<precision>&      second);
      template<typename precision>
      vector<vector<precision>>  operator- (const vector<vector<precision>>& first,
                                                    const vector<vector<precision>>& second);
      template<typename precision>
      vector<vector<precision>>  operator- (const vector<vector<precision>>& first,
                                                    const vector<precision>&      second);
      template<typename precision, typename ValueT>
      vector<precision>       operator- (const vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator- (const vector<vector<precision>>& vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<precision>       operator- (const ValueT&                            value,
                                                    const vector<precision>&      vec   );
      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator- (const ValueT&                            value,
                                                    const vector<vector<precision>>& vec   );

      template<typename precision>
      vector<precision>&      operator-=(      vector<precision>&      first,
                                                    const vector<precision>&      second);
      template<typename precision>
      vector<vector<precision>>& operator-=(      vector<vector<precision>>& first,
                                                    const vector<vector<precision>>& second);
      template<typename precision>
      vector<vector<precision>>& operator-=(      vector<vector<precision>>& first,
                                                    const vector<precision>&      second);
      template<typename precision, typename ValueT>
      vector<precision>&      operator-=(      vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator-=(      vector<vector<precision>>& vec,
                                                    const ValueT&                            value );

      template<typename precision, typename ValueT>
      vector<precision>       operator* (const vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator* (const vector<vector<precision>>& vec,
                                                    const ValueT&                            value );
      template<typename precision>
      vector<vector<precision>>  operator* (const vector<vector<precision>>& vec,
                                                    const vector<precision>&        value );
      template<typename precision, typename ValueT>
      vector<precision>       operator* (const ValueT&                            value,
                                                    const vector<precision>&      vec   );
      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator* (const ValueT&                            value,
                                                    const vector<vector<precision>>& vec   );

      template<typename precision, typename ValueT>
      vector<precision>&      operator*=(      vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator*=(      vector<vector<precision>>& vec,
                                                    const ValueT&                            value );

      template<typename precision, typename ValueT>
      vector<precision>       operator/ (const vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>  operator/ (const vector<vector<precision>>& vec,
                                                    const ValueT&                            value );
      template<typename precision>
      vector<vector<precision>>  operator/ (const vector<vector<precision>>& vec,
                                                    const vector<precision>&        value );

      template<typename precision, typename ValueT>
      vector<precision>&      operator/=(      vector<precision>&      vec,
                                                    const ValueT&                            value );
      template<typename precision, typename ValueT>
      vector<vector<precision>>& operator/=(      vector<vector<precision>>& vec,
                                                    const ValueT&                            value );
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#include "particle_util_impl.hpp"

#endif  // _EXAMPLES__BORIS__PARTICLE_UTIL_HPP_
