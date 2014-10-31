
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
     * This IMEX sweeper is for ODEs of the form \\( \\dot{U} = F_{\\rm expl}(t,U) + F_{\\rm impl}(t, U)\\).  To reduce
     * complexity and computational effort the non-stiff part is treated explicitly and the stiff part implicitly.
     *
     * This sweeper requires three interfaces to be implemented: two routines to evaluate the explicit \\( F_{\\rm expl}
     * \\) and implicit \\( F_{\\rm impl} \\) pieces for a given state, and one routine that solves (perhaps with an
     * external solver) the backward-Euler equation \\( U^{n+1} - \\Delta t F_{\\rm impl}(U^{n+1}) = RHS \\) for \\(
     * U^{n+1} \\).
     *
     * @tparam time precision type of the time dimension
     */
    template<typename time = time_precision>
    class IMEXSweeper
      : public pfasst::encap::EncapSweeper<time>
    {
      protected:
        vector<time> aug_nodes;


        //! @{
        /**
         * Solution values \\( U \\) at all time nodes of the current iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> state;

        /**
         * Solution values \\( U \\) at all time nodes of the previous iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> saved_state;

        /**
         * Node-to-node integrals of \\( F(t,u) \\) at all time nodes of the current iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> s_integrals;

        /**
         * FAS corrections \\( \\tau \\) at all time nodes of the current iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fas_corrections;

        /**
         * Values of the explicit part of the right hand side \\( F_{expl}(t,u) \\) at all time nodes of the current
         * iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fs_expl;

        /**
         * Values of the implicit part of the right hand side \\( F_{impl}(t,u) \\) at all time nodes of the current
         * iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fs_impl;
        //! @}

        //! @{
        /**
         * Augmented solution and function values.
         *
         * For quadrature rules that do not include the left end point (0.0), we add a virtual (augmented) node.
         */
        vector<shared_ptr<Encapsulation<time>>> aug_state;
        vector<shared_ptr<Encapsulation<time>>> aug_fs_impl;
        vector<shared_ptr<Encapsulation<time>>> aug_fs_expl;
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
         * Augmented quadrature matrix.
         */
        Matrix<time> s_mat_augm;

        /**
         * Compact quadrature matrix.
         */
        Matrix<time> q_mat_cmpt;
        //! @}

        /**
        * Set end state to \\( U_0 + \\int F_{expl} + F_{expl} \\).
        */
        void integrate_end_state(time dt)
        {
          vector<shared_ptr<Encapsulation<time>>> dst = { this->end_state };
          dst[0]->copy(this->start_state);
          dst[0]->mat_apply(dst, dt, this->quad->get_b_mat(), this->fs_expl, false);
          dst[0]->mat_apply(dst, dt, this->quad->get_b_mat(), this->fs_impl, false);
        }

      public:
        //! @{
        IMEXSweeper() = default;

        virtual ~IMEXSweeper() = default;
        //! @}

        //! @{
        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const override
        {
          return this->state[m];
        }

        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const override
        {
          return this->fas_corrections[m];
        }

        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const override
        {
          return this->saved_state[m];
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
          //
          // nodes
          //
          auto const nodes = this->quad->get_nodes();
          auto const num_nodes = this->quad->get_num_nodes();

          this->aug_nodes = nodes;
          if (!this->quad->left_is_node()) {
            this->aug_nodes.insert(this->aug_nodes.begin(), 0.0);
          }
          auto const num_aug_nodes = aug_nodes.size();

          //
          // quadrature matrices
          //
          auto const s_mat = this->quad->get_s_mat();
          auto const q_mat = this->quad->get_q_mat();
          this->s_mat_augm = Matrix<time>::Zero(num_aug_nodes-1, num_aug_nodes);
          if (this->quad->left_is_node()) {
            this->s_mat_augm.block(0, 0, num_nodes-1, num_nodes) = s_mat.block(1, 0, num_nodes-1, num_nodes);
            this->q_mat_cmpt = q_mat.block(1, 0, num_nodes-1, num_nodes);
          } else {
            this->s_mat_augm.block(0, 1, num_nodes, num_nodes) = s_mat.block(0, 0, num_nodes, num_nodes);
            this->q_mat_cmpt = q_mat;
          }

          this->s_mat_expl = s_mat_augm;
          this->s_mat_impl = s_mat_augm;
          for (size_t m = 0; m < aug_nodes.size() - 1; m++) {
            time ds = aug_nodes[m + 1] - aug_nodes[m];
            this->s_mat_expl(m, m)     -= ds;
            this->s_mat_impl(m, m + 1) -= ds;
          }

          //
          // allocate
          //
          this->start_state = this->get_factory()->create(pfasst::encap::solution);
          this->end_state = this->get_factory()->create(pfasst::encap::solution);
          for (size_t m = 0; m < num_nodes; m++) {
            this->state.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->saved_state.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
            this->fs_expl.push_back(this->get_factory()->create(pfasst::encap::function));
            this->fs_impl.push_back(this->get_factory()->create(pfasst::encap::function));
          }

          for (size_t m = 0; m < num_aug_nodes-1; m++) {
            this->s_integrals.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->fas_corrections.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
          }

          //
          // augment
          //
          this->aug_state = this->state;
          this->aug_fs_expl = this->fs_expl;
          this->aug_fs_impl = this->fs_impl;

          if (!this->quad->left_is_node()) {
            this->aug_state.insert(this->aug_state.begin(), this->start_state);
            this->aug_fs_expl.insert(this->aug_fs_expl.begin(), this->get_factory()->create(pfasst::encap::function));
            // XXX: this should really be a nullptr... but then we need to some fancy mat_apply magic
            this->aug_fs_impl.insert(this->aug_fs_impl.begin(), this->get_factory()->create(pfasst::encap::function));
          }

          assert(this->state.size() == this->quad->get_num_nodes());
          assert(this->fs_expl.size() == this->quad->get_num_nodes());
          assert(this->fs_impl.size() == this->quad->get_num_nodes());
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
          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          if (initial) {
            this->aug_state[0]->copy(this->start_state);
            this->f_expl_eval(this->aug_fs_expl[0], this->aug_state[0], t);
            if (this->quad->left_is_node()) {
              this->f_impl_eval(this->aug_fs_impl[0], this->aug_state[0], t);
            }
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          // step across all nodes
          for (size_t m = 0; m < this->aug_nodes.size() - 1; ++m) {
            time ds = dt * (this->aug_nodes[m+1] - this->aug_nodes[m]);
            rhs->copy(this->aug_state[m]);
            rhs->saxpy(ds, this->aug_fs_expl[m]);
            this->impl_solve(this->aug_fs_impl[m + 1], this->aug_state[m + 1], t, ds, rhs);
            this->f_expl_eval(this->aug_fs_expl[m + 1], this->aug_state[m + 1], t + ds);
            t += ds;
          }

          // set end state
          if (this->quad->right_is_node()) {
            this->end_state->copy(this->aug_state.back());
          } else {
            this->integrate_end_state(dt);
          }
        }

        /**
         * Perform one SDC sweep/iteration.
         *
         * This computes a high-order solution from the previous iteration's function values and
         * corrects it using forward/backward Euler steps across the nodes.
         */
        virtual void sweep() override
        {
          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->s_mat_expl, this->aug_fs_expl, true);
          this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->s_mat_impl, this->aug_fs_impl, false);
          if (this->fas_corrections.size() > 0) {
            for (size_t m = 0; m < this->s_integrals.size(); m++) {
              this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m]);
            }
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          // step across all nodes
          for (size_t m = 0; m < this->aug_nodes.size() - 1; ++m) {
            time ds = dt * (aug_nodes[m+1] - aug_nodes[m]);
            rhs->copy(this->aug_state[m]);
            rhs->saxpy(ds, this->aug_fs_expl[m]);
            rhs->saxpy(1.0, this->s_integrals[m]);
            this->impl_solve(this->aug_fs_impl[m + 1], this->aug_state[m + 1], t, ds, rhs);
            this->f_expl_eval(this->aug_fs_expl[m + 1], this->aug_state[m + 1], t + ds);
            t += ds;
          }

          // set end state
          if (this->quad->right_is_node()) {
            this->end_state->copy(this->aug_state.back());
          } else {
            this->integrate_end_state(dt);
          }
        }

        /**
         * Advance the end solution to start solution.
         */
        virtual void advance() override
        {
          // solutions
          this->start_state->copy(this->end_state);
          if (this->quad->left_is_node()) {
            this->aug_state[0]->copy(this->end_state);
          }
          // functions
          if (this->quad->left_is_node() && this->quad->right_is_node()) {
            this->aug_fs_expl[0]->copy(this->aug_fs_expl.back());
            this->aug_fs_impl[0]->copy(this->aug_fs_impl.back());
          } else if (this->quad->right_is_node()) {
            this->aug_fs_expl[0]->copy(this->fs_expl.back());
          } else {
            time t0 = this->get_controller()->get_time();
            time dt = this->get_controller()->get_time_step();
            this->f_expl_eval(this->aug_fs_expl[0], this->start_state, t0+dt);
          }
        }

        /**
         * Save current solution states.
         */
        virtual void save(bool initial_only) override
        {
          if (initial_only) {
            this->saved_state[0]->copy(state[0]);
          } else {
            for (size_t m = 0; m < this->saved_state.size(); m++) {
              this->saved_state[m]->copy(state[m]);
            }
          }
        }

        /**
         * @copybrief EncapSweeper::evaluate()
         */
        virtual void evaluate(size_t m) override
        {
          time t0 = this->get_controller()->get_time();
          time dt = this->get_controller()->get_time_step();
          time t =  t0 + dt * this->quad->get_nodes()[m];
          this->f_expl_eval(this->fs_expl[m], this->state[m], t);
          this->f_impl_eval(this->fs_impl[m], this->state[m], t);
        }

        /**
         * @copybrief EncapSweeper::integrate()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst integrated values; will get zeroed out beforehand
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override
        {
          dst[0]->mat_apply(dst, dt, this->q_mat_cmpt, this->fs_expl, true);
          dst[0]->mat_apply(dst, dt, this->q_mat_cmpt, this->fs_impl, false);
        }

        /**
         * @copybrief EncapSweeper::residual()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst residuals
         */
        virtual void residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override
        {
          for (size_t m=0; m<this->quad->get_num_nodes(); m++) {
            dst[m]->copy(this->get_start_state());
            dst[m]->saxpy(-1.0, this->get_state(m));
            // XXX: add tau corrections
          }
          dst[0]->mat_apply(dst, dt, this->q_mat_cmpt, this->fs_expl, false);
          dst[0]->mat_apply(dst, dt, this->q_mat_cmpt, this->fs_impl, false);
        }
        //! @}

        //! @{
        /**
         * Evaluate the explicit part of the ODE.
         *
         * @param[in,out] f_expl_encap Encapsulation to store the explicit function evaluation.
         * @param[in] u_encap Encapsulation that stores the solution state at which to evaluate the
         *     explicit part of the ODE.
         * @param[in] t Time point of the evaluation.
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
         * Evaluate the implicit part of the ODE.
         *
         * This is typically called to compute the implicit part of the right hand side at the first
         * collocation node, and on all nodes after restriction or interpolation.
         *
         * @param[in,out] f_impl_encap Encapsulation to store the implicit function evaluation.
         * @param[in] u_encap Encapsulation storing the solution state at which to evaluate the
         *     implicit part of the ODE.
         * @param[in] t Time point of the evaluation.
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
         * Solve \\( U - \\Delta t F_{\\rm impl}(U) = RHS \\) for \\( U \\).
         *
         * During an IMEX SDC sweep, the correction equation is evolved using a forward-Euler
         * stepper for the explicit piece, and a backward-Euler stepper for the implicit piece.
         * This routine (implemented by the user) performs the solve required to perform one
         * backward-Euler sub-step, and also returns \\( F_{\\rm impl}(U) \\).
         *
         * @param[in,out] f_encap Encapsulation to store the evaluated implicit piece.
         * @param[in,out] u_encap Encapsulation to store the solution of the backward-Euler sub-step.
         * @param[in] t time point (of \\( RHS \\)).
         * @param[in] dt sub-step size to the previous time point (\\( \\Delta t \\)).
         * @param[in] rhs_encap Encapsulation that stores \\( RHS \\).
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
