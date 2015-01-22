
#ifndef _PFASST_ENCAP_IMEX_SWEEPER_HPP_
#define _PFASST_ENCAP_IMEX_SWEEPER_HPP_

#include <cstdlib>
#include <cassert>
#include <vector>
#include <memory>

#include "../globals.hpp"
#include "../logging.hpp"
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

        //! @{
        /**
         * Node-to-node integrals of \\( F(t,u) \\) at all time nodes of the current iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> s_integrals;

        /**
         * Values of the explicit part of the right hand side \\( F_{expl}(t,u) \\) at all time nodes of the current
         * iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fs_expl;
        shared_ptr<Encapsulation<time>> fs_expl_start;

        /**
         * Values of the implicit part of the right hand side \\( F_{impl}(t,u) \\) at all time nodes of the current
         * iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fs_impl;
        //! @}

        /**
        * Set end state to \\( U_0 + \\int F_{expl} + F_{expl} \\).
        */
        void integrate_end_state(time dt)
        {
          vector<shared_ptr<Encapsulation<time>>> dst = { this->end_state };
          dst[0]->copy(this->start_state);
          dst[0]->mat_apply(dst, dt, this->quadrature->get_b_mat(), this->fs_expl, false);
          dst[0]->mat_apply(dst, dt, this->quadrature->get_b_mat(), this->fs_impl, false);
        }

      public:
        //! @{
        IMEXSweeper() = default;

        virtual ~IMEXSweeper() = default;
        //! @}

        //! @{
        /**
         * @copydoc ISweeper::setup(bool)
         */
        virtual void setup(bool coarse) override
        {
          pfasst::encap::EncapSweeper<time>::setup(coarse);

          auto const num_nodes = this->quadrature->get_num_nodes();
          auto const num_s_integrals = this->quadrature->left_is_node() ? num_nodes - 1 : num_nodes;

          for (size_t m = 0; m < num_nodes; m++) {
            this->fs_expl.push_back(this->get_factory()->create(pfasst::encap::function));
            this->fs_impl.push_back(this->get_factory()->create(pfasst::encap::function));
          }
          for (size_t m = 0; m < num_s_integrals; m++) {
            this->s_integrals.push_back(this->get_factory()->create(pfasst::encap::solution));
          }

          if (! this->quadrature->left_is_node()) {
            this->fs_expl_start = this->get_factory()->create(pfasst::encap::function);
          }
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
          if (this->quadrature->left_is_node()) {
            this->predict_with_left(initial);
          } else {
            this->predict_without_left(initial);
          }

          if (this->quadrature->right_is_node()) {
            this->end_state->copy(this->state.back());
          } else {
            this->integrate_end_state(this->get_controller()->get_time_step());
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
          if (this->quadrature->left_is_node()) {
            this->sweep_with_left();
          } else {
            this->sweep_without_left();
          }

          if (this->quadrature->right_is_node()) {
            this->end_state->copy(this->state.back());
          } else {
            this->integrate_end_state(this->get_controller()->get_time_step());
          }
        }

        /**
         * Advance the end solution to start solution.
         */
        virtual void advance() override
        {
          this->start_state->copy(this->end_state);
          if (this->quadrature->left_is_node() && this->quadrature->right_is_node()) {
            this->state[0]->copy(this->start_state);
            this->fs_expl.front()->copy(this->fs_expl.back());
            this->fs_impl.front()->copy(this->fs_impl.back());
          }
        }

        /**
         * @copybrief EncapSweeper::reevaluate()
         */
        virtual void reevaluate(bool initial_only) override
        {
          time t0 = this->get_controller()->get_time();
          time dt = this->get_controller()->get_time_step();
          if (initial_only) {
            if (this->quadrature->left_is_node()) {
              this->f_expl_eval(this->fs_expl[0], this->state[0], t0);
              this->f_impl_eval(this->fs_impl[0], this->state[0], t0);
            } else {
              throw NotImplementedYet("reevaluate");
            }
          } else {
            for (size_t m = 0; m < this->quadrature->get_num_nodes(); m++) {
              time t =  t0 + dt * this->quadrature->get_nodes()[m];
              this->f_expl_eval(this->fs_expl[m], this->state[m], t);
              this->f_impl_eval(this->fs_impl[m], this->state[m], t);
            }
          }
        }

        /**
         * @copybrief EncapSweeper::integrate()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst integrated values; will get zeroed out beforehand
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override
        {
          auto const q_mat = this->quadrature->get_q_mat();
          dst[0]->mat_apply(dst, dt, q_mat, this->fs_expl, true);
          dst[0]->mat_apply(dst, dt, q_mat, this->fs_impl, false);
        }

        /**
         * @copybrief EncapSweeper::residual()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst residuals
         */
        virtual void residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override
        {
          for (size_t m = 0; m < this->quadrature->get_num_nodes(); m++) {
            dst[m]->copy(this->start_state);
            dst[m]->saxpy(-1.0, this->state[m]);
          }
          if (this->fas_corrections.size() > 0) {
            // XXX: build a matrix and use mat_apply to do this
            for (size_t m = 0; m < this->quadrature->get_num_nodes(); m++) {
              for (size_t n = 0; n <= m; n++) {
                dst[m]->saxpy(1.0, this->fas_corrections[n]);
              }
            }
          }
          dst[0]->mat_apply(dst, dt, this->quadrature->get_q_mat(), this->fs_expl, false);
          dst[0]->mat_apply(dst, dt, this->quadrature->get_q_mat(), this->fs_impl, false);
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


      private:
        void predict_with_left(bool initial)
        {
          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();
          CLOG(INFO, "Sweeper") << "predicting step " << this->get_controller()->get_step() + 1
                                << " (t=" << t << ", dt=" << dt << ")";

          if (initial) {
            this->state[0]->copy(this->start_state);
            this->f_expl_eval(this->fs_expl[0], this->state[0], t);
            this->f_impl_eval(this->fs_impl[0], this->state[0], t);
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          auto const nodes = this->quadrature->get_nodes();

          // step across all nodes
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            time ds = dt * (nodes[m+1] - nodes[m]);
            rhs->copy(this->state[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            this->impl_solve(this->fs_impl[m + 1], this->state[m + 1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m + 1], this->state[m + 1], t + ds);
            t += ds;
          }
        }

        void predict_without_left(bool initial)
        {
          UNUSED(initial);
          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();
          CLOG(INFO, "Sweeper") << "predicting step " << this->get_controller()->get_step() + 1
                                << " (t=" << t << ", dt=" << dt << ")";
          time ds;

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          auto const nodes = this->quadrature->get_nodes();

          // step to first node
          ds = dt * nodes[0];
          this->f_expl_eval(this->fs_expl_start, this->start_state, t);
          rhs->copy(this->start_state);
          rhs->saxpy(ds, this->fs_expl_start);
          this->impl_solve(this->fs_impl[0], this->state[0], t, ds, rhs);
          this->f_expl_eval(this->fs_expl[0], this->state[0], t + ds);

          // step across all nodes
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            rhs->copy(this->state[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            this->impl_solve(this->fs_impl[m+1], this->state[m+1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m+1], this->state[m+1], t + ds);
            t += ds;
          }
        }

        void sweep_with_left()
        {
          auto const nodes = this->quadrature->get_nodes();
          auto const dt    = this->get_controller()->get_time_step();
          auto const s_mat = this->quadrature->get_s_mat().block(1, 0, nodes.size()-1, nodes.size());
          CLOG(INFO, "Sweeper") << "sweeping on step " << this->get_controller()->get_step() + 1
                                << " in iteration " << this->get_controller()->get_iteration() << " (dt=" << dt << ")";
          time ds;

          this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_expl, true);
          this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_impl, false);
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            this->s_integrals[m]->saxpy(-ds, this->fs_expl[m]);
            this->s_integrals[m]->saxpy(-ds, this->fs_impl[m+1]);
          }
          if (this->fas_corrections.size() > 0) {
            for (size_t m = 0; m < this->s_integrals.size(); m++) {
              this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m+1]);
            }
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          // step across all nodes
          auto t = this->get_controller()->get_time();
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            rhs->copy(this->state[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            rhs->saxpy(1.0, this->s_integrals[m]);
            this->impl_solve(this->fs_impl[m + 1], this->state[m + 1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m + 1], this->state[m + 1], t + ds);
            t += ds;
          }
        }

        void sweep_without_left()
        {
          auto const nodes = this->quadrature->get_nodes();
          auto const dt    = this->get_controller()->get_time_step();
          auto const s_mat = this->quadrature->get_s_mat();
          CLOG(INFO, "Sweeper") << "sweeping on step " << this->get_controller()->get_step() + 1
                                << " in iteration " << this->get_controller()->get_iteration() << " (dt=" << dt << ")";
          time ds;

          this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_expl, true);
          this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_impl, false);
          ds = dt * nodes[0];
          this->s_integrals[0]->saxpy(-ds, this->fs_expl_start);
          this->s_integrals[0]->saxpy(-ds, this->fs_impl[0]);
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            this->s_integrals[m+1]->saxpy(-ds, this->fs_expl[m]);
            this->s_integrals[m+1]->saxpy(-ds, this->fs_impl[m+1]);
          }
          if (this->fas_corrections.size() > 0) {
            for (size_t m = 0; m < this->s_integrals.size(); m++) {
              this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m]);
            }
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          // step to first node
          auto t = this->get_controller()->get_time();
          ds = dt * nodes[0];
          this->f_expl_eval(this->fs_expl_start, this->start_state, t);
          rhs->copy(this->start_state);
          rhs->saxpy(ds, this->fs_expl_start);
          rhs->saxpy(1.0, this->s_integrals[0]);
          this->impl_solve(this->fs_impl[0], this->state[0], t, ds, rhs);
          this->f_expl_eval(this->fs_expl[0], this->state[0], t + ds);

          // step across all nodes
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            rhs->copy(this->state[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            rhs->saxpy(1.0, this->s_integrals[m+1]);
            this->impl_solve(this->fs_impl[m+1], this->state[m+1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m+1], this->state[m+1], t + ds);
            t += ds;
          }
        }

    };

  }  // ::pfasst::encap
}  // ::pfasst

#endif
