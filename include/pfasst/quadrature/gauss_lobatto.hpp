#ifndef _PFASST__QUADRATURE__GAUSS_LOBATTO_HPP_
#define _PFASST__QUADRATURE__GAUSS_LOBATTO_HPP_

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
    class GaussLobatto;

    template<typename precision>
    struct quadrature_traits<GaussLobatto<precision>>
    {
      typedef pfasst::quadrature::gauss_lobatto integral_constant;
      static const bool left_is_node = true;
      static const bool right_is_node = true;
    };

    template<typename precision>
    class GaussLobatto
      : public IQuadrature<precision>
    {
      public:
        //! @{
        GaussLobatto(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->m_num_nodes < 2) {
            throw invalid_argument("Gauss-Lobatto quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
        }

        GaussLobatto()
          : IQuadrature<precision>()
        {}

        GaussLobatto(const GaussLobatto<precision>& other)
          : IQuadrature<precision>(other)
        {}

        GaussLobatto(GaussLobatto<precision>&& other)
          : GaussLobatto<precision>()
        {
          swap(*this, other);
        }

        virtual ~GaussLobatto()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return quadrature_traits<GaussLobatto<precision>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<GaussLobatto<precision>>::right_is_node; }
        //! @}

        //! @{
        GaussLobatto<precision>& operator=(GaussLobatto<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override
        {
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          auto roots = Polynomial<precision>::legendre(this->m_num_nodes - 1).differentiate().roots();

          for (size_t j = 0; j < this->m_num_nodes - 2; j++) {
            this->m_nodes[j + 1] = 0.5 * (1.0 + roots[j]);
          }
          this->m_nodes.front() = 0.0;
          this->m_nodes.back() = 1.0;
        }
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__GAUSS_LOBATTO_HPP_
