#include "pfasst/quadrature/uniform.hpp"

#include <stdexcept>
#include <vector>
using namespace std;


namespace pfasst
{
  namespace quadrature
  {
    template<typename precision>
    Uniform<precision>::Uniform(const size_t num_nodes)
      : IQuadrature<precision>(num_nodes)
    {
      if (this->num_nodes < 2) {
        throw invalid_argument("Uniform quadrature requires at least two quadrature nodes.");
      }
      this->compute_nodes();
      this->compute_weights();
    }

    template<typename precision>
    bool Uniform<precision>::left_is_node() const
    {
      return LEFT_IS_NODE;
    }

    template<typename precision>
    bool Uniform<precision>::right_is_node() const
    {
      return RIGHT_IS_NODE;
    }

    template<typename precision>
    void Uniform<precision>::compute_nodes()
    {
      this->nodes = vector<precision>(this->num_nodes, precision(0.0));
      for (size_t j = 0; j < this->num_nodes; j++) {
        this->nodes[j] = precision(j) / (this->num_nodes - 1);
      }
    }
  }  // ::pfasst::quadrature
}  // ::pfasst
