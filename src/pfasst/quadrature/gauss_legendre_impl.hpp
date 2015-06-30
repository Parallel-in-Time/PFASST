#include "pfasst/quadrature/gauss_legendre.hpp"

#include <vector>
using namespace std;

#include "pfasst/quadrature/polynomial.hpp"


namespace pfasst
{
  namespace quadrature
  {
    template<typename precision>
    GaussLegendre<precision>::GaussLegendre(const size_t num_nodes)
      : IQuadrature<precision>(num_nodes)
    {
      this->compute_nodes();
      this->compute_weights();
    }

    template<typename precision>
    bool GaussLegendre<precision>::left_is_node() const
    {
      return LEFT_IS_NODE;
    }

    template<typename precision>
    bool GaussLegendre<precision>::right_is_node() const
    {
      return RIGHT_IS_NODE;
    }

    template<typename precision>
    void GaussLegendre<precision>::compute_nodes()
    {
      this->nodes = vector<precision>(this->num_nodes, precision(0.0));
      auto roots = Polynomial<precision>::legendre(this->num_nodes).roots();

      for (size_t j = 0; j < this->num_nodes; j++) {
        this->nodes[j] = 0.5 * (1.0 + roots[j]);
      }
    }
  }  // ::pfasst::quadrature
}  // ::pfasst
