#ifndef _PFASST__ENCAP__TRAITS_HPP_
#define _PFASST__ENCAP__TRAITS_HPP_

#include <vector>
using namespace std;


namespace pfasst
{
  template<
    class TimePrecision,
    class SpacialPrecision,
    class DataT,
    class... Ts
  >
  struct encap_traits
  {
    typedef TimePrecision    time_type;
    typedef SpacialPrecision spacial_type;
    typedef DataT            data_type;
  };


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
