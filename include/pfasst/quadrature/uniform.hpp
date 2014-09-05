#ifndef _PFASST__QUADRATURE__UNIFORM_HPP_
#define _PFASST__QUADRATURE__UNIFORM_HPP_

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
    class Uniform
      : public IQuadrature<precision>
    {
      protected:
        //! @{
        static const bool LEFT_IS_NODE = true;
        static const bool RIGHT_IS_NODE = true;
        //! @}

      public:
        //! @{
        Uniform(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->num_nodes < 2) {
            throw invalid_argument("Uniform quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
          this->compute_delta_nodes();
        }

        Uniform()
          : IQuadrature<precision>()
        {}

        Uniform(const Uniform<precision>& other)
          : IQuadrature<precision>(other)
        {}

        Uniform(Uniform<precision>&& other)
          : Uniform<precision>()
        {
          swap(*this, other);
        }

        virtual ~Uniform()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return LEFT_IS_NODE; }

        virtual bool right_is_node() const
        { return RIGHT_IS_NODE; }
        //! @}

        //! @{
        Uniform<precision>& operator=(Uniform<precision> other)
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
          for (size_t j = 0; j < this->num_nodes; j++) {
            this->nodes[j] = precision(j) / (this->num_nodes - 1);
          }
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__UNIFORM_HPP_
