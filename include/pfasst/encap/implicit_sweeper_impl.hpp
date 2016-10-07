#include "pfasst/encap/implicit_sweeper.hpp"

#include <cassert>

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {
    template<typename scalar>
    using lu_pair = pair< Matrix<scalar>, Matrix<scalar> >;

    /**
     * LU (without pivoting) decomposition.
     */
    template<typename scalar>
    static lu_pair<scalar> lu_decomposition(const Matrix<scalar>& A)
    {
      assert(A.rows() == A.cols());

      auto n = A.rows();

      Matrix<scalar> L = Matrix<scalar>::Zero(n, n);
      Matrix<scalar> U = Matrix<scalar>::Zero(n, n);

      if (A.rows() == 1) {

        L(0, 0) = 1.0;
        U(0, 0) = A(0,0);

      } else {

        // first row of U is first row of A
        auto U12 = A.block(0, 1, 1, n-1);

        // first column of L is first column of A / a11
        auto L21 = A.block(1, 0, n-1, 1) / A(0, 0);

        // remove first row and column and recurse
        auto A22  = A.block(1, 1, n-1, n-1);
        Matrix<scalar> tmp = A22 - L21 * U12;
        auto LU22 = lu_decomposition(tmp);

        L(0, 0) = 1.0;
        U(0, 0) = A(0, 0);
        L.block(1, 0, n-1, 1) = L21;
        U.block(0, 1, 1, n-1) = U12;
        L.block(1, 1, n-1, n-1) = get<0>(LU22);
        U.block(1, 1, n-1, n-1) = get<1>(LU22);

      }

      return lu_pair<scalar>(L, U);
    }

    /**
     * Augment nodes: nodes <- [t0] + dt * nodes
     */
    template<typename time>
    vector<time> augment(time t0, time dt, vector<time> const & nodes)
    {
      vector<time> t(1 + nodes.size());
      t[0] = t0;
      for (size_t m = 0; m < nodes.size(); m++) {
        t[m+1] = t0 + dt * nodes[m];
      }
      return t;
    }

    /*
     * Implementations
     */

    template<typename time>
    void ImplicitSweeper<time>::set_end_state()
    {
      if (this->quadrature->right_is_node()) {
        this->end_state->copy(this->state.back());
      } else {
        vector<shared_ptr<Encapsulation<time>>> dst = { this->end_state };
        dst[0]->copy(this->start_state);
        dst[0]->mat_apply(dst, this->get_controller()->get_step_size(), this->quadrature->get_b_mat(), this->fs_impl, false);
      }
    }

    template<typename time>
    void ImplicitSweeper<time>::setup(bool coarse)
    {
      pfasst::encap::EncapSweeper<time>::setup(coarse);

      auto const nodes = this->quadrature->get_nodes();
      auto const num_nodes = this->quadrature->get_num_nodes();

      if (this->quadrature->left_is_node()) {
        ML_CLOG(INFO, "Sweeper", "implicit sweeper shouldn't include left endpoint");
        throw ValueError("implicit sweeper shouldn't include left endpoint");
      }

      for (size_t m = 0; m < num_nodes; m++) {
        this->s_integrals.push_back(this->get_factory()->create(pfasst::encap::solution));
        this->fs_impl.push_back(this->get_factory()->create(pfasst::encap::function));
      }

      Matrix<time> QT = this->quadrature->get_q_mat().transpose();
      auto lu = lu_decomposition(QT);
      auto L = get<0>(lu);
      auto U = get<1>(lu);
      this->q_tilde = U.transpose();

      ML_CLOG(DEBUG, "Sweeper", "Q':" << endl << QT);
      ML_CLOG(DEBUG, "Sweeper", "L:" << endl << L);
      ML_CLOG(DEBUG, "Sweeper", "U:" << endl << U);
      ML_CLOG(DEBUG, "Sweeper", "LU:" << endl << L * U);
      ML_CLOG(DEBUG, "Sweeper", "q_tilde:" << endl << this->q_tilde);
    }

    template<typename time>
    void ImplicitSweeper<time>::predict(bool initial)
    {
      UNUSED(initial);

      auto const dt = this->get_controller()->get_step_size();
      auto const t  = this->get_controller()->get_time();

      ML_CLOG(INFO, "Sweeper", "predicting step " << this->get_controller()->get_step() + 1
                               << " (t=" << t << ", dt=" << dt << ")");

      auto const anodes = augment(t, dt, this->quadrature->get_nodes());
      for (size_t m = 0; m < anodes.size() - 1; ++m) {
        this->impl_solve(this->fs_impl[m], this->state[m], anodes[m], anodes[m+1] - anodes[m],
                         m == 0 ? this->get_start_state() : this->state[m-1]);
      }

      this->set_end_state();
    }

    template<typename time>
    void ImplicitSweeper<time>::sweep()
    {
      auto const dt = this->get_controller()->get_step_size();
      auto const t  = this->get_controller()->get_time();

      ML_CLOG(INFO, "Sweeper", "sweeping on step " << this->get_controller()->get_step() + 1
                               << " in iteration " << this->get_controller()->get_iteration()
                               << " (dt=" << dt << ")");

      this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->quadrature->get_s_mat(), this->fs_impl, true);
      if (this->fas_corrections.size() > 0) {
        for (size_t m = 0; m < this->s_integrals.size(); m++) {
          this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m]);
        }
      }

      for (size_t m = 0; m < this->s_integrals.size(); m++) {
        for (size_t n = 0; n < m; n++) {
          this->s_integrals[m]->saxpy(-dt*this->q_tilde(m, n), this->fs_impl[n]);
        }
      }

      shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

      auto const anodes = augment(t, dt, this->quadrature->get_nodes());
      for (size_t m = 0; m < anodes.size() - 1; ++m) {
        auto const ds = anodes[m+1] - anodes[m];
        rhs->copy(m == 0 ? this->get_start_state() : this->state[m-1]);
        rhs->saxpy(1.0, this->s_integrals[m]);
        rhs->saxpy(-ds, this->fs_impl[m]);
        for (size_t n = 0; n < m; n++) {
          rhs->saxpy(dt*this->q_tilde(m, n), this->fs_impl[n]);
        }
        this->impl_solve(this->fs_impl[m], this->state[m], anodes[m], ds, rhs);
      }
      this->set_end_state();
    }

    template<typename time>
    void ImplicitSweeper<time>::advance()
    {
      this->start_state->copy(this->end_state);
    }

    template<typename time>
    void ImplicitSweeper<time>::reevaluate(bool initial_only)
    {
      if (initial_only) {
        return;
      }
      auto const dt = this->get_controller()->get_step_size();
      auto const t0 = this->get_controller()->get_time();
      auto const nodes = this->quadrature->get_nodes();
      for (size_t m = 0; m < nodes.size(); m++) {
        this->f_impl_eval(this->fs_impl[m], this->state[m], t0 + dt * nodes[m]);
      }
    }

    template<typename time>
    void ImplicitSweeper<time>::integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
    {
      dst[0]->mat_apply(dst, dt, this->quadrature->get_q_mat(), this->fs_impl, true);
    }

  }
}
