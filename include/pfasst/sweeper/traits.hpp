#ifndef _PFASST__SWEEPER__TRAITS_HPP_
#define _PFASST__SWEEPER__TRAITS_HPP_

#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  template<
    class EncapsulationTraits,
    class... Ts
  >
  struct sweeper_traits
  {
    typedef          EncapsulationTraits                       encap_traits;
    typedef          encap::Encapsulation<EncapsulationTraits> encap_type;
    typedef typename encap_traits::time_type                   time_type;
    typedef typename encap_traits::spacial_type                spacial_type;
  };
}  // ::pfasst

#endif  // _PFASST__SWEEPER__TRAITS_HPP_
