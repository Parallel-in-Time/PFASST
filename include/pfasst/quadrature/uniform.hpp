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
    class Uniform;

    template<typename precision>
    struct quadrature_traits<Uniform<precision>>
    {
      typedef pfasst::quadrature::uniform integral_constant;
      static const bool left_is_node = true;
      static const bool right_is_node = true;
    };

    template<typename precision>
    class Uniform
      : public IQuadrature<precision>
    {
      public:
        //! @{
        Uniform(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->m_num_nodes < 2) {
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
        { return quadrature_traits<Uniform<precision>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<Uniform<precision>>::right_is_node; }
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
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          for (size_t j = 0; j < this->m_num_nodes; j++) {
            this->m_nodes[j] = precision(j) / (this->m_num_nodes - 1);
          }
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__UNIFORM_HPP_
