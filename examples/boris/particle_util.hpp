/**
 * @ingroup BorisFiles
 * @file examples/boris/particle_util.hpp
 */
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
      /**
       * @defgroup BorisUtilities Utilities
       * @ingroup Boris
       */

      //! @{
      // TODO: check whether this function is still required
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> cloud_component_factory(const size_t num_particles, const size_t dim);

      //! @ingroup BorisUtilities
      template<typename precision>
      inline static void zero(vector<precision>& data);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static void zero(shared_ptr<vector<precision>>& data);
      //! @}

      //! @{
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> cross_prod(const vector<precision>& first,
                                                 const vector<precision>& second);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static void cross_prod_1part(typename vector<precision>::const_iterator __first,
                                          typename vector<precision>::const_iterator __second,
                                          typename vector<precision>::iterator __result);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> cross_prod_npart(const vector<precision>& first,
                                                       const vector<precision>& second);


      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> kronecker(const vector<precision>& first,
                                                const vector<precision>& second);


      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> cmp_wise_mul(const vector<precision>& first,
                                                   const vector<precision>& second);

      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> cmp_wise_div(const vector<precision>& first,
                                                   const vector<precision>& second);


      //! @ingroup BorisUtilities
      template<typename precision>
      inline static precision max(const vector<precision>& data);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static precision max_abs(const vector<precision>& data);


      //! @ingroup BorisUtilities
      template<typename precision>
      inline static precision norm_sq(const vector<precision>& data);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static precision norm_sq(typename vector<precision>::const_iterator __first,
                                      typename vector<precision>::const_iterator __second);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> norm_sq_npart(const vector<precision>& data, const size_t npart);

      //! @ingroup BorisUtilities
      template<typename precision>
      inline static precision norm0(const vector<precision>& data);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static precision norm0(typename vector<precision>::const_iterator __first,
                                    typename vector<precision>::const_iterator __second);
      //! @ingroup BorisUtilities
      template<typename precision>
      inline static vector<precision> norm0_npart(const vector<precision>& data, const size_t npart);
      //! @}


      //! @{
      //! @ingroup BorisUtilities
      template<typename precision>
      inline vector<precision>  operator+ (const vector<precision>& first,
                                           const vector<precision>& second);
      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>  operator+ (const vector<precision>& vec,
                                           const ValueT&            value );
      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>  operator+ (const ValueT&            value,
                                           const vector<precision>& vec   );

      //! @ingroup BorisUtilities
      template<typename precision>
      inline vector<precision>& operator+=(      vector<precision>& first,
                                           const vector<precision>& second);
      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>& operator+=(      vector<precision>& vec,
                                           const ValueT&            value );

      //! @ingroup BorisUtilities
      template<typename precision>
      inline vector<precision>  operator- (const vector<precision>& first,
                                           const vector<precision>& second);
      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>  operator- (const vector<precision>& vec,
                                           const ValueT&            value );

      //! @ingroup BorisUtilities
      template<typename precision>
      inline vector<precision>& operator-=(      vector<precision>& first,
                                           const vector<precision>& second);
      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>& operator-=(      vector<precision>& vec,
                                           const ValueT&            value );

      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>  operator* (const vector<precision>& vec,
                                           const ValueT&            value );
      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>  operator* (const ValueT&            value,
                                           const vector<precision>& vec   );
      //! @ingroup BorisUtilities
      template<typename precision>
      inline vector<precision>  operator* (const vector<precision>& vec,
                                           const vector<precision>& values);

      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>& operator*=(      vector<precision>& vec,
                                           const ValueT&            value );

      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>  operator/ (const vector<precision>& vec,
                                           const ValueT&            value );
      //! @ingroup BorisUtilities
      template<typename precision>
      inline vector<precision>  operator/ (const vector<precision>& vec,
                                           const vector<precision>& values);

      //! @ingroup BorisUtilities
      template<typename precision, typename ValueT>
      inline vector<precision>& operator/=(      vector<precision>& vec,
                                           const ValueT&            value );
      //! @}
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#include "particle_util_impl.hpp"

#endif  // _EXAMPLES__BORIS__PARTICLE_UTIL_HPP_
