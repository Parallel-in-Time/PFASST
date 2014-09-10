
#ifndef _PFASST_ENCAP_IMEX_SWEEPER_HPP_
#define _PFASST_ENCAP_IMEX_SWEEPER_HPP_

#include <cstdlib>
#include <cassert>
#include <vector>
#include <memory>

#include "../globals.hpp"
#include "../quadrature.hpp"
#include "encapsulation.hpp"
#include "encap_sweeper.hpp"
#include "vector.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {
    using pfasst::encap::Encapsulation;

    /**
     * Semi-implicit IMEX sweeper.
     *
     * Given an ODE \\( \\frac{\\partial}{\\partial t}u(t) = F(t,u) \\) where the function of the
     * right hand side \\( F(t,u) \\) can be split into a non-stiff and a stiff part.
     * To reduce complexity and computational effort one would want to solve the non-stiff part
     * explicitly and the stiff part implicitly.
     * Therefore, we define the splitting \\( F(t,u) = F_{expl}(t,u) + F_{impl}(t,u) \\).
     *
     * This sweeper requires three interfaces to be implement for such ODEs: two routines to
     * evaluate the explicit \\( F_{\\rm expl} \\) and implicit \\( F_{\\rm impl} \\) pieces for a
     * given state, and one that solves (perhaps with an external solver) the backward-Euler
     * equation \\( U^{n+1} - \\Delta t F_{\\rm impl}(U^{n+1}) = RHS \\) for \\( U^{n+1} \\).
     *
     * @tparam time precision type of the time dimension
     */
    template<typename time = time_precision>
    class IMEXSweeper
      : public pfasst::encap::EncapSweeper<time>
    {
      protected:
        //! @{
        /**
         * solution values \\( u(\\tau) \\) at all time nodes \\( \\tau \\in [0, M-1] \\) of the
         * current iteration
         */
        vector<shared_ptr<Encapsulation<time>>> u_state;

        /**
         * solution values \\( u(t) \\) at all time nodes \\( t \\in [0, M-1] \\) of the
         * previous iteration
         */
        vector<shared_ptr<Encapsulation<time>>> previous_u_state;

        /**
         * node-to-node integrated values of \\( F(t,u) \\) at all time nodes \\( t \\in
         * [0, M-1] \\) of the current iteration
         */
        vector<shared_ptr<Encapsulation<time>>> s_integrals;

        /**
         * FAS corrections \\( \\tau_t \\) at all time nodes \\( t \\in [0, M-1] \\) of the current
         * iteration
         */
        vector<shared_ptr<Encapsulation<time>>> fas_corrections;

        /**
         * values of the explicit part of the right hand side \\( F_{expl}(t,u) \\) at all time
         * nodes \\( t \\in [0, M-1] \\) of the current iteration
         */
        vector<shared_ptr<Encapsulation<time>>> fs_expl;

        shared_ptr<Encapsulation<time>> fs_expl_start;
        // shared_ptr<Encapsulation<time>> fs_expl_end;

        /**
         * values of the implicit part of the right hand side \\( F_{impl}(t,u) \\) at all time
         * nodes \\( t \\in [0, M-1] \\) of the current iteration
         */
        vector<shared_ptr<Encapsulation<time>>> fs_impl;

        // shared_ptr<Encapsulation<time>> fs_impl_start;
        // shared_ptr<Encapsulation<time>> fs_impl_end;
        //! @}

        //! @{
        /**
         * Quadrature matrix containing weights for node-to-node integration for the explicit part.
         *
         * @see IMEXSweeper::setup(bool) for a short description
         */
        Matrix<time> s_mat_expl;

        /**
         * Quadrature matrix containing weights for node-to-node integration of the implicit part.
         *
         * @see IMEXSweeper::setup(bool) for a short description
         */
        Matrix<time> s_mat_impl;

        /**
         * Compact quadrature matrix.
         */
        Matrix<time> s_mat_cmpt;
        //! @}

        /**
        * Set end state to \\( U_0 + \\int F_{expl} + F_{expl} \\).
        */
        void integrate_end_state(time dt)
        {
          vector<shared_ptr<Encapsulation<time>>> dst = { this->u_end };
          dst[0]->copy(this->u_start);
          dst[0]->mat_apply(dst, dt, this->quad->get_b_mat(), this->fs_expl, false);
          dst[0]->mat_apply(dst, dt, this->quad->get_b_mat(), this->fs_impl, false);
        }

      public:
        //! @{
        IMEXSweeper() = default;

        virtual ~IMEXSweeper() = default;
        //! @}

        //! @{
        virtual void set_state(shared_ptr<const Encapsulation<time>> u0, size_t m) override
        {
          this->u_state[m]->copy(u0);
        }

        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const override
        {
          return this->u_state[m];
        }

        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const override
        {
          return this->fas_corrections[m];
        }

        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const override
        {
          return this->previous_u_state[m];
        }
        //! @}

        //! @{
        /**
         * @copydoc ISweeper::setup(bool)
         *
         * To reduce computational overhead, we precompute the partial integration matrices
         * IMEXSweeper::s_mat_expl \\( (\\tilde{s}^{expl})_{m,j} \\) and IMEXSweeper::s_mat_impl
         * \\( (\\tilde{s}^{impl})_{m,j} \\) by incorporating known values of the SDC sweep
         * equation.
         *
         * Let \\( F = F_{impl} + F_{expl} \\), \\( f = F_{impl} \\), \\( g = F_{expl} \\) and
         * \\( \\Delta t_m = t_{m+1} - t_m \\).
         * @f{eqnarray*}{
         *   u_{m+1}^{k+1} &=& u_m^k + \Delta t_m \left( f_{m+1}^{k+1} - f_{m+1}^k \right)
         *                     + \Delta t_m \left( g_m^{k+1} - g_m^k \right)
         *                     + \sum_{j=1}^M s_{m,j} F_j^k \\
         *                 &=& u_m^k + \Delta t_m f_{m+1}^{k+1} + \Delta t_m g_m^{k+1}
         *                     + \sum_{j=1}^M s_{m,j} F_j^k
         *                     - \Delta t_m \left( f_{m+1}^k - g_m^k \right) \\
         *                 &=& u_m^k + \Delta t_m f_{m+1}^{k+1} + \Delta t_m g_m^{k+1}
         *                     + \sum_{j=1}^M s_{m,j} f_j^k - \Delta t_m f_{m+1}^k
         *                     + \sum_{j=1}^M s_{m,j} g_j^k - \Delta t_m g_m^k \\
         *                 &=& u_m^k + \Delta t_m f_{m+1}^{k+1} + \Delta t_m g_m^{k+1}
         *                     + \sum_{j=1}^M \tilde{s}_{m,j}^{impl} f_j^k
         *                     + \sum_{j=1}^M \tilde{s}_{m,j}^{expl} g_j^k
         * @f}
         * with
         * @f[
         *  \tilde{s}_{m,j}^{impl} =
         *    \begin{cases}
         *      s_{m,j} - \Delta t_m &\mbox{if } j \equiv m+1 \\
         *      s_{m,j} &\mbox{else}
         *    \end{cases} \\
         *  \tilde{s}_{m,j}^{expl} =
         *    \begin{cases}
         *      s_{m,j} - \Delta t_m &\mbox{if } j \equiv m \\
         *      s_{m,j} &\mbox{else}
         *    \end{cases} \\
         * @f]
         *
         * This procedure can be derived from the so called \\( Q_{\\Delta} \\) notation.
         */
        virtual void setup(bool coarse) override
        {
          const auto nodes = this->quad->get_nodes();
          const size_t num_nodes = this->quad->get_num_nodes();
          const auto delta_nodes = this->quad->get_delta_nodes();

          auto s_mat = this->quad->get_s_mat();

          if (this->quad->left_is_node()) {
            this->s_mat_cmpt = s_mat.block(1, 0, num_nodes-1, num_nodes);
          } else {
            this->s_mat_cmpt = s_mat;
          }

          this->s_mat_expl = s_mat_cmpt;
          this->s_mat_impl = s_mat_cmpt;
          // for (size_t m = 1; m < s_mat_cmpt.rows(); ++m) {
          //   this->s_mat_expl(m, m - 1) -= nodes[m] - nodes[m-1];
          //   this->s_mat_impl(m, m)     -= nodes[m] - nodes[m-1];
          //   // this->s_mat_expl(m, m - 1) -= delta_nodes[m];
          //   // this->s_mat_impl(m, m)     -= delta_nodes[m];
          // }

         for (size_t m = 0; m < nodes.size() - 1; m++) {
            time ds = nodes[m + 1] - nodes[m];
            this->s_mat_expl(m, m)     -= ds;
            this->s_mat_impl(m, m + 1) -= ds;
          }

          this->u_start = this->get_factory()->create(pfasst::encap::solution);
          this->u_end = this->get_factory()->create(pfasst::encap::solution);
          for (size_t m = 0; m < num_nodes; m++) {
            this->u_state.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->previous_u_state.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
            this->fs_expl.push_back(this->get_factory()->create(pfasst::encap::function));
            this->fs_impl.push_back(this->get_factory()->create(pfasst::encap::function));

          }

          auto num_non_zero = this->quad->left_is_node() ? num_nodes - 1 : num_nodes;
          for (size_t m = 0; m < num_non_zero; m++) {
            this->s_integrals.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->fas_corrections.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
          }

          if (!this->quad->left_is_node()) {
            this->fs_expl_start = this->get_factory()->create(pfasst::encap::function);
          }
          // this->fs_expl_end = this->get_factory()->create(pfasst::encap::function);
          // this->fs_impl_start = this->get_factory()->create(pfasst::encap::function);
          // this->fs_impl_end = this->get_factory()->create(pfasst::encap::function);

          assert(this->u_state.size() == this->quad->get_num_nodes());
          assert(this->fs_expl.size() == this->quad->get_num_nodes());
          // assert(this->fs_impl.size() == this->quad->get_num_nodes());
          // assert(this->s_integrals.size() == this->quad->get_num_nodes());
        }

        /**
         * Compute low-order provisional solution.
         *
         * This performs forward/backward Euler steps across the nodes to compute a low-order provisional solution.
         *
         * @param[in] initial if `true` the explicit and implicit part of the right hand side of the
         *     ODE get evaluated with the initial value
         */
        virtual void predict(bool initial) override
        {
          const auto   nodes       = this->quad->get_nodes();
          const size_t num_nodes   = this->quad->get_num_nodes();
          const auto   delta_nodes = this->quad->get_delta_nodes();

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          if (initial) {
            if (this->quad->left_is_node()) {
              this->u_state[0]->copy(this->u_start);
              this->f_impl_eval(this->fs_impl[0], this->u_start, t);
              this->f_expl_eval(this->fs_expl[0], this->u_start, t);
            } else {
              this->f_expl_eval(this->fs_expl_start, this->u_start, t);
            }
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          // step to first node if necessary
          if (!this->quad->left_is_node()) {
            time ds = dt * nodes[0];
            rhs->copy(this->u_start);
            rhs->saxpy(ds, this->fs_expl_start);
            this->impl_solve(this->fs_impl[0], this->u_state[0], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[0], this->u_state[0], t + ds);
          }

          // step across all nodes
          for (size_t m = 0; m < num_nodes - 1; ++m) {
            time ds = dt * (nodes[m+1] - nodes[m]);
            rhs->copy(this->u_state[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            this->impl_solve(this->fs_impl[m + 1], this->u_state[m + 1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m + 1], this->u_state[m + 1], t + ds);
            t += ds;
          }

          // set end state
          if (this->quad->right_is_node()) {
            this->u_end->copy(this->u_state.back());
          } else {
            this->integrate_end_state(dt);
          }
        }

        virtual void sweep() override
        {
          const auto   nodes  = this->quad->get_nodes();
          const size_t num_nodes = nodes.size();
          const auto delta_nodes = this->quad->get_delta_nodes();

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->s_mat_expl, this->fs_expl, true);
          this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->s_mat_impl, this->fs_impl, false);
          if (this->fas_corrections.size() > 0) {
            for (size_t m = 0; m < num_nodes; m++) {
              this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m]);
            }
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          // step to first node if necessary
          if (!this->quad->left_is_node()) {
            time ds = dt * nodes[0];
            rhs->copy(this->u_start);
            rhs->saxpy(ds, this->fs_expl_start);
            rhs->saxpy(1.0, this->s_integrals[0]);
            this->impl_solve(this->fs_impl[0], this->u_state[0], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[0], this->u_state[0], t + ds);
          }

          // step across all nodes
          for (size_t m = 0; m < num_nodes - 1; ++m) {
            time ds = dt * (nodes[m+1] - nodes[m]);
            rhs->copy(this->u_state[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            rhs->saxpy(1.0, this->s_integrals[m]);
            this->impl_solve(this->fs_impl[m + 1], this->u_state[m + 1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m + 1], this->u_state[m + 1], t + ds);
            t += ds;
          }

          // set end state
          if (this->quad->right_is_node()) {
            this->u_end->copy(this->u_state.back());
          } else {
            this->integrate_end_state(dt);
          }
        }

        virtual void advance() override
        {
          if (this->quad->left_is_node() && this->quad->right_is_node()) {
            this->u_state[0]->copy(this->u_end);
            this->fs_expl.front()->copy(this->fs_expl.back());
            this->fs_impl.front()->copy(this->fs_impl.back());
          } else {
            throw NotImplementedYet("imex advance");
          }
        }

        virtual void save(bool initial_only) override
        {
          if (initial_only) {
            this->previous_u_state[0]->copy(u_state[0]);
          } else {
            for (size_t m = 0; m < this->previous_u_state.size(); m++) {
              this->previous_u_state[m]->copy(u_state[m]);
            }
          }
        }

        /**
         * @copybrief EncapSweeper::evaluate()
         *
        * If the node `m` is virtual (not proper), we can save some
        * implicit/explicit evaluations.
        */
        virtual void evaluate(size_t m) override
        {
          time t0 = this->get_controller()->get_time();
          time dt = this->get_controller()->get_time_step();
          time t =  t0 + dt * this->quad->get_nodes()[m];
          if (m == 0) {
            this->f_expl_eval(this->fs_expl.front(), this->u_state.front(), t);
          }
          this->f_impl_eval(this->fs_impl[m], this->u_state[m], t);
        }

        /**
         * @copybrief EncapSweeper::integrate()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst integrated values; will get zeroed out beforehand
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override
        {
          dst[0]->mat_apply(dst, dt, this->quad->get_s_mat(), this->fs_expl, true);
          dst[0]->mat_apply(dst, dt, this->quad->get_s_mat(), this->fs_impl, false);
        }
        //! @}

        //! @{
        /**
         * Evaluates the explicit part of the right hand side of the ODE at the given time.
         *
         * @param[in,out] f_expl_encap Encapsulation to store the evaluated right hand side
         * @param[in] u_encap Encapsulation storing the solution values to use for computing the
         *     explicit part of the right hand side of the ODE
         * @param[in] t time point of the evaluation
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void f_expl_eval(shared_ptr<Encapsulation<time>> f_expl_encap,
                                 shared_ptr<Encapsulation<time>> u_encap,
                                 time t)
        {
          UNUSED(f_expl_encap); UNUSED(u_encap); UNUSED(t);
          throw NotImplementedYet("imex (f_expl_eval)");
        }

        /**
         * Evaluates the implicit part of the right hand side of the ODE at the given time.
         *
         * This is typically called to compute the implicit part of the right hand side at the first
         * collocation node, and on all nodes after restriction or interpolation.
         *
         * @param[in,out] f_impl_encap Encapsulation to store the evaluated right hand side
         * @param[in] u_encap Encapsulation storing the solution values to use for computing the
         *     implicit part of the right hand side of the ODE
         * @param[in] t time point of the evaluation
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void f_impl_eval(shared_ptr<Encapsulation<time>> f_impl_encap,
                                 shared_ptr<Encapsulation<time>> u_encap,
                                 time t)
        {
          UNUSED(f_impl_encap); UNUSED(u_encap); UNUSED(t);
          throw NotImplementedYet("imex (f_impl_eval)");
        }

        /**
         * Solves \\( U - \\Delta t f_{\\rm impl}(U) = RHS \\) for \\( U \\).
         *
         * During an IMEX SDC sweep, the correction equation is evolved using a forward-Euler
         * stepper for the explicit piece, and a backward-Euler stepper for the implicit piece.
         * This routine (implemented by the user) performs the solve required to perform one
         * backward-Euler sub-step, and also returns \\( f_{\\rm impl}(U) \\).
         *
         * @param[in,out] f_encap Encapsulation to store the evaluated right hand side
         * @param[in,out] u_encap Encapsulation to store the solution of the backward-Euler sub-step
         * @param[in] t time point (of \\( RHS \\))
         * @param[in] dt sub-step size to the previous time point (\\( \\Delta t \\))
         * @param[in] rhs_encap Encapsulation storing \\( RHS \\)
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void impl_solve(shared_ptr<Encapsulation<time>> f_encap,
                                shared_ptr<Encapsulation<time>> u_encap,
                                time t, time dt,
                                shared_ptr<Encapsulation<time>> rhs_encap)
        {
          UNUSED(f_encap); UNUSED(u_encap); UNUSED(t); UNUSED(dt); UNUSED(rhs_encap);
          throw NotImplementedYet("imex (impl_solve)");
        }
        //! @}
    };

  }  // ::pfasst::encap
}  // ::pfasst

#endif
