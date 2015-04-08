/**
 * @file pfasst/quadrature/clenshaw_curtis.hpp
 * @since v0.3.0
 */
#ifndef _PFASST__QUADRATURE__CLENSHAW_CURTIS_HPP_
#define _PFASST__QUADRATURE__CLENSHAW_CURTIS_HPP_

#include <cassert>
#include <vector>
using namespace std;

#include "pfasst/interfaces.hpp"
#include "pfasst/quadrature/interface.hpp"


namespace pfasst
{
  namespace quadrature
  {
    /**
     * Quadrature handler for Clenshaw-Curtis quadrature.
     *
     * @tparam scalar precision of quadrature (i.e. `double`)
     *
     * @since v0.3.0
     */
    template<typename precision = pfasst::time_precision>
    class ClenshawCurtis
      : public IQuadrature<precision>
    {
      protected:
        //! @{
        static const bool LEFT_IS_NODE = true;
        static const bool RIGHT_IS_NODE = true;
        //! @}

      public:
        //! @{
        /**
         * @throws invalid_argument if less than two nodes are requested
         */
        explicit ClenshawCurtis(const size_t num_nodes);
        ClenshawCurtis() = default;
        virtual ~ClenshawCurtis() = default;
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

#include "pfasst/quadrature/clenshaw_curtis_impl.hpp"

#endif  // _PFASST__QUADRATURE__CLENSHAW_CURTIS_HPP_
