
#ifndef _PFASST_ENCAP_IMPLICIT_SWEEPER_HPP_
#define _PFASST_ENCAP_IMPLICIT_SWEEPER_HPP_

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
     * Implicit sweeper.
     *
     * @tparam time precision type of the time dimension
     */
    template<typename time = time_precision>
    class ImplicitSweeper
      : public pfasst::encap::EncapSweeper<time>
    {
      protected:
        //! @{
        /**
         * Node-to-node integrals of \\( F(t,u) \\) at all time nodes of the current iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> s_integrals;

        /**
         * FAS corrections \\( \\tau \\) at all time nodes of the current iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fas_corrections;

        /**
         * Values of the implicit part of the right hand side \\( F_{impl}(t,u) \\) at all time nodes of the current
         * iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fs_impl;
        //! @}

        /**
        * Set end state to \\( U_0 + \\int F_{expl} + F_{expl} \\).
        */
        void set_end_state(time dt)
        {
          if (this->quadrature->right_is_node()) {
            this->end_state->copy(this->state.back());
          } else {
            vector<shared_ptr<Encapsulation<time>>> dst = { this->end_state };
            dst[0]->copy(this->start_state);
            dst[0]->mat_apply(dst, dt, this->quadrature->get_b_mat(), this->fs_impl, false);
          }
        }

      public:
        //! @{
        ImplicitSweeper() = default;
        virtual ~ImplicitSweeper() = default;
        //! @}

        //! @{
        /**
         * @copydoc ISweeper::setup(bool)
         *
         */
        virtual void setup(bool coarse) override
        {
          pfasst::encap::EncapSweeper<time>::setup(coarse);

          auto const nodes = this->quadrature->get_nodes();
          auto const num_nodes = this->quadrature->get_num_nodes();

          if (this->quadrature->left_is_node()) {
            throw ValueError("implicit sweeper shouldn't include left endpoint");
          }

          for (size_t m = 0; m < num_nodes; m++) {
            this->s_integrals.push_back(this->get_factory()->create(pfasst::encap::solution));
            this->fs_impl.push_back(this->get_factory()->create(pfasst::encap::function));
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
          UNUSED(initial);

          CLOG(INFO, "Sweeper") << "predicting step " << this->get_controller()->get_step() + 1
                                << " (t=" << this->get_controller()->get_time()
                                << ", dt=" << this->get_controller()->get_time_step() << ")";

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();
          time ds;

          auto const nodes = this->quadrature->get_nodes();

          ds = dt * nodes[0];
          this->impl_solve(this->fs_impl[0], this->state[0], t, ds, this->get_start_state());

          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            this->impl_solve(this->fs_impl[m+1], this->state[m+1], t, ds, this->state[m]);
            t += ds;
          }

          this->set_end_state(dt);
        }

        /**
         * Perform one SDC sweep/iteration.
         *
         * This computes a high-order solution from the previous iteration's function values and
         * corrects it using forward/backward Euler steps across the nodes.
         */
        virtual void sweep() override
        {
          CLOG(INFO, "Sweeper") << "sweeping on step " << this->get_controller()->get_step() + 1
                                << " in iteration " << this->get_controller()->get_iteration()
                                << " (dt=" << this->get_controller()->get_time_step() << ")";

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();
          time ds;
          auto const nodes = this->quadrature->get_nodes();


          this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->quadrature->get_s_mat(), this->fs_impl, true);

          ds = dt * nodes[0];
          this->s_integrals[0]->saxpy(-ds, this->fs_impl[0]);
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            this->s_integrals[m+1]->saxpy(-ds, this->fs_impl[m+1]);
          }
          if (this->fas_corrections.size() > 0) {
            for (size_t m = 0; m < this->s_integrals.size(); m++) {
              this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m]);
            }
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          // step to first node
          ds = dt * nodes[0];
          rhs->copy(this->start_state);
          rhs->saxpy(1.0, this->s_integrals[0]);
          this->impl_solve(this->fs_impl[0], this->state[0], t, ds, rhs);

          // step across all nodes
          for (size_t m = 0; m < nodes.size() - 1; ++m) {
            ds = dt * (nodes[m+1] - nodes[m]);
            rhs->copy(this->state[m]);
            rhs->saxpy(1.0, this->s_integrals[m+1]);
            this->impl_solve(this->fs_impl[m+1], this->state[m+1], t, ds, rhs);
            t += ds;
          }

          this->set_end_state(dt);
        }

        /**
         * Advance the end solution to start solution.
         */
        virtual void advance() override
        {
          this->start_state->copy(this->end_state);
        }

        /**
         * @copybrief EncapSweeper::evaluate()
         */
        virtual void reevaluate(bool initial_only) override
        {
          if (initial_only) {
            return;
          }

          time t0 = this->get_controller()->get_time();
          time dt = this->get_controller()->get_time_step();
          for (size_t m = 0; m < this->quadrature->get_num_nodes(); m++) {
            time t =  t0 + dt * this->quadrature->get_nodes()[m];
            this->f_impl_eval(this->fs_impl[m], this->state[m], t);
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
          dst[0]->mat_apply(dst, dt, this->quadrature->get_q_mat(), this->fs_impl, true);
        }
        //! @}

        //! @{
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
         * During an Implicit SDC sweep, the correction equation is evolved using a forward-Euler
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
