/**
 * @file pfasst/quadrature/gauss_legendre.hpp
 * @since v0.3.0
 */
#ifndef _PFASST__QUADRATURE__GAUSS_LEGENDRE_HPP_
#define _PFASST__QUADRATURE__GAUSS_LEGENDRE_HPP_

#include "pfasst/quadrature/interface.hpp"


namespace pfasst
{
  namespace quadrature
  {
    /**
     * Quadrature handler for Gauss-Legendre quadrature.
     *
     * @tparam scalar precision of quadrature (i.e. `double`)
     *
     * @since v0.3.0
     */
    template<typename precision = pfasst::time_precision>
    class GaussLegendre
      : public IQuadrature<precision>
    {
      protected:
        //! @{
        static const bool LEFT_IS_NODE = false;
        static const bool RIGHT_IS_NODE = false;
        //! @}

      public:
        //! @{
        explicit GaussLegendre(const size_t num_nodes);
        GaussLegendre() = default;
        virtual ~GaussLegendre() = default;
        //! @}

        //! @{
        virtual bool left_is_node() const override;
        virtual bool right_is_node() const override;
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override;
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#include "pfasst/quadrature/gauss_legendre_impl.hpp"

#endif  // _PFASST__QUADRATURE__GAUSS_LEGENDRE_HPP_
