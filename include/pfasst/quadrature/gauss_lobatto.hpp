#ifndef _PFASST__QUADRATURE__GAUSS_LOBATTO_HPP_
#define _PFASST__QUADRATURE__GAUSS_LOBATTO_HPP_

#include "pfasst/quadrature/interface.hpp"


namespace pfasst
{
  namespace quadrature
  {
    template<typename precision = pfasst::time_precision>
    class GaussLobatto
      : public IQuadrature<precision>
    {
      protected:
        //! @{
        static const bool LEFT_IS_NODE = true;
        static const bool RIGHT_IS_NODE = true;
        //! @}

      public:
        //! @{
        explicit GaussLobatto(const size_t num_nodes);
        GaussLobatto() = default;
        virtual ~GaussLobatto() = default;
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

#include "pfasst/quadrature/gauss_lobatto_impl.hpp"

#endif  // _PFASST__QUADRATURE__GAUSS_LOBATTO_HPP_
