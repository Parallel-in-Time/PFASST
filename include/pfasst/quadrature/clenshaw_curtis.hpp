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
    class ClenshawCurtis;

    template<>
    struct quadrature_traits<ClenshawCurtis<>>
    {
      typedef pfasst::quadrature::clenshaw_curtis integral_constant;
      static const bool left_is_node = true;
      static const bool right_is_node = true;
    };

    template<typename precision>
    class ClenshawCurtis
      : public IQuadrature<precision>
    {
      public:
        //! @{
        ClenshawCurtis(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->m_num_nodes < 2) {
            throw invalid_argument("Clenshaw-Curtis quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
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
        { return quadrature_traits<ClenshawCurtis<>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<ClenshawCurtis<>>::right_is_node; }
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
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          auto roots = Polynomial<precision>::legendre(this->m_num_nodes).roots();

          for (size_t j = 0; j < this->m_num_nodes; j++) {
            this->m_nodes[j] = 0.5 * (1.0 - cos(j * pi<precision>() / (this->m_num_nodes - 1)));
          }
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__CLENSHAW_CURTIS_HPP_
