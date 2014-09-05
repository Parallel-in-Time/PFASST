#ifndef _PFASST__QUADRATURE__CLENSHAW_CURTIS_HPP_
#define _PFASST__QUADRATURE__CLENSHAW_CURTIS_HPP_

#include <cassert>
#include <vector>

#include <boost/math/constants/constants.hpp>
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
using namespace boost::math::constants;


namespace pfasst
{
  namespace quadrature
  {
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
        ClenshawCurtis(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->num_nodes < 2) {
            throw invalid_argument("Clenshaw-Curtis quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
          this->compute_delta_nodes();
        }

        ClenshawCurtis()
          : IQuadrature<precision>()
        {}

        ClenshawCurtis(const ClenshawCurtis<precision>& other)
          : IQuadrature<precision>(other)
        {}

        ClenshawCurtis(ClenshawCurtis<precision>&& other)
          : ClenshawCurtis<precision>()
        {
          swap(*this, other);
        }

        virtual ~ClenshawCurtis()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return LEFT_IS_NODE; }

        virtual bool right_is_node() const
        { return RIGHT_IS_NODE; }
        //! @}

        //! @{
        ClenshawCurtis<precision>& operator=(ClenshawCurtis<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes()
        {
          this->nodes = vector<precision>(this->num_nodes, precision(0.0));
          auto roots = Polynomial<precision>::legendre(this->num_nodes).roots();

          for (size_t j = 0; j < this->num_nodes; j++) {
            this->nodes[j] = 0.5 * (1.0 - cos(j * pi<precision>() / (this->num_nodes - 1)));
          }
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__CLENSHAW_CURTIS_HPP_
