#ifndef _PFASST__TRANSFER__TRAITS_HPP_
#define _PFASST__TRANSFER__TRAITS_HPP_

#include <type_traits>
using namespace std;


namespace pfasst
{
  template<
    class CoarseSweeper,
    class FineSweeper,
    class... Ts
  >
  struct transfer_traits
  {
    typedef          CoarseSweeper                       coarse_sweeper_type;
    typedef typename conditional<
                       is_void<coarse_sweeper_type>::value,
                       void,
                       typename coarse_sweeper_type::traits
                     >::type                             coarse_sweeper_traits;
    typedef typename conditional<
                       is_void<coarse_sweeper_traits>::value,
                       void,
                       typename coarse_sweeper_traits::encap_traits
                     >::type                             coarse_encap_traits;
    typedef typename conditional<
                       is_void<coarse_sweeper_traits>::value,
                       void,
                       typename coarse_sweeper_traits::encap_type
                     >::type                             coarse_encap_type;
    typedef typename conditional<
                       is_void<coarse_sweeper_traits>::value,
                       void,
                       typename coarse_sweeper_traits::time_type
                     >::type                             coarse_time_type;
    typedef typename conditional<
                       is_void<coarse_sweeper_traits>::value,
                       void,
                       typename coarse_sweeper_traits::spacial_type
                     >::type                             coarse_spacial_type;

    typedef          FineSweeper                         fine_sweeper_type;
    typedef typename fine_sweeper_type::traits           fine_sweeper_traits;
    typedef typename fine_sweeper_traits::encap_traits   fine_encap_traits;
    typedef typename fine_sweeper_traits::encap_type     fine_encap_type;
    typedef typename fine_sweeper_traits::time_type      fine_time_type;
    typedef typename fine_sweeper_traits::spacial_type   fine_spacial_type;

    static const size_t num_levels = conditional<
                                       is_void<coarse_sweeper_type>::value,
                                       integral_constant<size_t, 1>,
                                       integral_constant<size_t, 2>
                                     >::type::value;
  };
}  // ::pfasst

#endif  // _PFASST__TRANSFER__TRAITS_HPP_
