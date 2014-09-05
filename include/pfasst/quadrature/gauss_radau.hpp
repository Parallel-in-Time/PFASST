#ifndef _PFASST__QUADRATURE__GAUSS_RADAU_HPP_
#define _PFASST__QUADRATURE__GAUSS_RADAU_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "../interfaces.hpp"
#include "polynomial.hpp"
#include "interface.hpp"
#include "traits.hpp"

template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename scalar>
using Index = typename Matrix<scalar>::Index;

using namespace std;


namespace pfasst
{
  namespace quadrature
  {
    template<typename precision = pfasst::time_precision>
    class GaussRadau
      : public IQuadrature<precision>
    {
      protected:
        //! @{
        static const bool LEFT_IS_NODE = false;
        static const bool RIGHT_IS_NODE = true;
        //! @}

      public:
        //! @{
        GaussRadau(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->num_nodes < 2) {
            throw invalid_argument("Gauss-Radau quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
          this->compute_delta_nodes();
        }

        GaussRadau()
          : IQuadrature<precision>()
        {}

        GaussRadau(const GaussRadau<precision>& other)
          : IQuadrature<precision>(other)
        {}

        GaussRadau(GaussRadau<precision>&& other)
          : GaussRadau<precision>()
        {
          swap(*this, other);
        }

        virtual ~GaussRadau()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return LEFT_IS_NODE; }

        virtual bool right_is_node() const
        { return RIGHT_IS_NODE; }
        //! @}

        //! @{
        GaussRadau<precision>& operator=(GaussRadau<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override
        {
          this->nodes = vector<precision>(this->num_nodes, precision(0.0));
          auto l   = Polynomial<precision>::legendre(this->num_nodes);
          auto lm1 = Polynomial<precision>::legendre(this->num_nodes - 1);

          for (size_t i = 0; i < this->num_nodes; i++) {
            l[i] += lm1[i];
          }
          auto roots = l.roots();
          for (size_t j = 1; j < this->num_nodes; j++) {
            this->nodes[j - 1] = 0.5 * (1.0 - roots[this->num_nodes - j]);
          }
          this->nodes.back() = 1.0;
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__GAUSS_RADAU_HPP_
