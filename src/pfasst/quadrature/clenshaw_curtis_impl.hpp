#include "pfasst/quadrature/clenshaw_curtis.hpp"

#include <stdexcept>
using namespace std;

#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;

#include "pfasst/quadrature/polynomial.hpp"


namespace pfasst
{
  namespace quadrature
  {
    template<typename precision>
    ClenshawCurtis<precision>::ClenshawCurtis(const size_t num_nodes)
      : IQuadrature<precision>(num_nodes)
    {
      if (this->num_nodes < 2) {
        throw invalid_argument("Clenshaw-Curtis quadrature requires at least two quadrature nodes.");
      }
      this->compute_nodes();
      this->compute_weights();
    }

    template<typename precision>
    bool ClenshawCurtis<precision>::left_is_node() const
    {
      return LEFT_IS_NODE;
    }

    template<typename precision>
    bool ClenshawCurtis<precision>::right_is_node() const
    {
      return RIGHT_IS_NODE;
    }

    template<typename precision>
    void ClenshawCurtis<precision>::compute_nodes()
    {
      this->nodes = vector<precision>(this->num_nodes, precision(0.0));
      auto roots = Polynomial<precision>::legendre(this->num_nodes).roots();

      for (size_t j = 0; j < this->num_nodes; j++) {
        this->nodes[j] = 0.5 * (1.0 - cos(j * pi<precision>() / (this->num_nodes - 1)));
      }
    }
  }  // ::pfasst::quadrature
}  // ::pfasst
