#include "pfasst/quadrature/gauss_radau.hpp"

#include <stdexcept>
#include <vector>
using namespace std;

#include "pfasst/logging.hpp"
#include "pfasst/quadrature/polynomial.hpp"


namespace pfasst
{
  namespace quadrature
  {
    template<typename precision>
    GaussRadau<precision>::GaussRadau(const size_t num_nodes)
      : IQuadrature<precision>(num_nodes)
    {
      if (this->num_nodes < 2) {
        throw invalid_argument("Gauss-Radau quadrature requires at least two quadrature nodes.");
      }
      this->compute_nodes();
      this->compute_weights();

      CLOG(INFO, "QUAD") << "Gauss-Radau (right) on [0.0, 1.0] with " << num_nodes
                         << " nodes at " << this->get_nodes();
      CLOG(DEBUG, "QUAD") << "Q:" << endl << this->get_q_mat();
    }

    template<typename precision>
    bool GaussRadau<precision>::left_is_node() const
    {
      return LEFT_IS_NODE;
    }

    template<typename precision>
    bool GaussRadau<precision>::right_is_node() const
    {
      return RIGHT_IS_NODE;
    }

    template<typename precision>
    void GaussRadau<precision>::compute_nodes()
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

    template<typename precision>
    void GaussRadau<precision>::compute_weights()
    {
      const size_t n = this->get_num_nodes();

      this->q_mat = Matrix<precision>::Zero(n + 1, n + 1);
      this->q_mat.block(1, 1, n, n) << compute_q_matrix(this->get_nodes());

      this->s_mat = compute_s_matrix(this->get_q_mat());

      this->q_vec = compute_q_vec(this->get_nodes());
      this->q_vec.insert(this->q_vec.begin(), precision(0.0));

      this->b_mat = Matrix<precision>(1, n + 1);
      for (size_t i = 0; i < n + 1; i++) {
        this->b_mat(0, i) = this->q_vec[i];
      }
    }
  }  // ::pfasst::quadrature
}  // ::pfasst
