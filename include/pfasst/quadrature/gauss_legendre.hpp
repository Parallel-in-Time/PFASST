#ifndef _PFASST__QUADRATURE__GAUSS_LEGENDRE_HPP_
#define _PFASST__QUADRATURE__GAUSS_LEGENDRE_HPP_

#include <cassert>
#include <vector>

#include "../interfaces.hpp"
#include "polynomial.hpp"
#include "interface.hpp"

using namespace std;


namespace pfasst
{
  namespace quadrature
  {
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
        explicit GaussLegendre(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          this->compute_nodes();
          this->compute_weights();
          this->compute_delta_nodes();
        }

        GaussLegendre() = default;

        virtual ~GaussLegendre() = default;
        //! @}

        //! @{
        virtual bool left_is_node() const { return LEFT_IS_NODE; }

        virtual bool right_is_node() const { return RIGHT_IS_NODE; }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override
        {
          this->nodes = vector<precision>(this->num_nodes, precision(0.0));
          auto roots = Polynomial<precision>::legendre(this->num_nodes).roots();

          for (size_t j = 0; j < this->num_nodes; j++) {
            this->nodes[j] = 0.5 * (1.0 + roots[j]);
          }
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__GAUSS_LEGENDRE_HPP_
