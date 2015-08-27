#ifndef _PFASST__ENCAP__TRAITS_HPP_
#define _PFASST__ENCAP__TRAITS_HPP_

#include <vector>
using namespace std;


namespace pfasst
{
  /**
   * Type Traits for encapsulation of user data types.
   *
   * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
   * @tparam SpacialPrecision the spacial data precision
   * @tparam DataT            the actual data type encapsulated
   *
   * @ingroup Traits
   */
  template<
    class TimePrecision,
    class SpacialPrecision,
    class DataT,
    class... Ts
  >
  struct encap_traits
  {
    //! public member type for the time precision
    typedef TimePrecision    time_type;

    //! public member type for the spacial precision
    typedef SpacialPrecision spacial_type;

    //! public member type for the encapsulated data type
    typedef DataT            data_type;
  };


  /**
   * Spacialized Type Traits for encapsulation of std::vector.
   *
   * @tparam TimePrecision    the time precision, e.g. precision of the integration nodes
   * @tparam SpacialPrecision the spacial data precision
   *
   * @ingroup Traits
   */
  template<
    class TimePrecision,
    class SpacialPrecision
  >
  struct vector_encap_traits
    : public encap_traits<TimePrecision, SpacialPrecision, vector<SpacialPrecision>>
  {
    typedef TimePrecision        time_type;
    typedef SpacialPrecision     spacial_type;
    typedef vector<spacial_type> data_type;
  };

}  // ::pfasst

#endif  // _PFASST__ENCAP__TRAITS_HPP_
