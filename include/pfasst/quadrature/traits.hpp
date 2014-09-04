#ifndef _PFASST__QUADRATURE__TRAITS_HPP_
#define _PFASST__QUADRATURE__TRAITS_HPP_

#include <type_traits>

using namespace std;

namespace pfasst
{
  namespace quadrature
  {
    enum class QuadratureType : int {
        GaussLegendre   =  0
      , GaussLobatto    =  1
      , GaussRadau      =  2
      , ClenshawCurtis  =  3
      , Uniform         =  4
      , UNDEFINED       = -1
    };

    typedef integral_constant<QuadratureType, QuadratureType::GaussLegendre> gauss_legendre;
    typedef integral_constant<QuadratureType, QuadratureType::GaussLobatto> gauss_lobatto;
    typedef integral_constant<QuadratureType, QuadratureType::GaussRadau> gauss_radau;
    typedef integral_constant<QuadratureType, QuadratureType::ClenshawCurtis> clenshaw_curtis;
    typedef integral_constant<QuadratureType, QuadratureType::Uniform> uniform;
    typedef integral_constant<QuadratureType, QuadratureType::UNDEFINED> undefined;

    template<typename QuadratureT>
    struct quadrature_traits
    {
      typedef pfasst::quadrature::undefined integral_constant;
      static const bool left_is_node = false;
      static const bool right_is_node = false;
    };
  }
}

#endif  // _PFASST__QUADRATURE__TRAITS_HPP_
