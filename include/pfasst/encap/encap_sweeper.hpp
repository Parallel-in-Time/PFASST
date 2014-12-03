/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAP_ENCAP_SWEEPER_HPP_
#define _PFASST_ENCAP_ENCAP_SWEEPER_HPP_

#include <algorithm>
#include <cstdlib>
#include <vector>
#include <memory>

#include "../globals.hpp"
#include "../interfaces.hpp"
#include "../quadrature.hpp"
#include "encapsulation.hpp"

namespace pfasst
{
  namespace encap
  {
    using namespace pfasst::quadrature;

    template<typename time = time_precision>
    class EncapSweeper
      : public ISweeper<time>
    {
      protected:
        //! @{
        shared_ptr<IQuadrature<time>> quadrature;
        shared_ptr<EncapFactory<time>> factory;
        shared_ptr<Encapsulation<time>> start_state;
        shared_ptr<Encapsulation<time>> end_state;
        vector<shared_ptr<Encapsulation<time>>> residuals;
        //! @}

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
         * FAS corrections \\( \\tau \\) at all time nodes of the current iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fas_corrections;
        //! @}

        int residual_norm_order;
        time abs_residual_tol, rel_residual_tol;

      public:

        EncapSweeper()
          : quadrature(nullptr), abs_residual_tol(0.0), rel_residual_tol(0.0)
        {}

        //! @{
        /**
         * Retrieve solution values of current iteration at time node index `m`.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const
        {
          return this->state[m];
        }

        /**
         * Retrieve FAS correction of current iteration at time node index `m`.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const
        {
          return this->fas_corrections[m];
        }

        /**
         * Retrieve solution values of previous iteration at time node index `m`.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const
        {
          return this->saved_state[m];
        }
        //! @}

        virtual void setup(bool coarse) override
        {
          auto const nodes = this->quadrature->get_nodes();
          auto const num_nodes = this->quadrature->get_num_nodes();

          this->start_state = this->get_factory()->create(pfasst::encap::solution);
          this->end_state = this->get_factory()->create(pfasst::encap::solution);

          for (size_t m = 0; m < num_nodes; m++) {
            this->state.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->saved_state.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
          }

          if (coarse) {
            for (size_t m = 0; m < num_nodes; m++) {
              this->fas_corrections.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
          }
        }

        //! @{
        virtual void spread() override
        {
          for (size_t m = 1; m < this->quadrature->get_num_nodes(); m++) {
            this->state[m]->copy(this->state[0]);
          }
        }
        //! @}

        /**
         * Save current solution states.
         */
        virtual void save(bool initial_only) override
        {
          // XXX: if !left_is_node, this is a problem...
          if (initial_only) {
            this->saved_state[0]->copy(state[0]);
          } else {
            for (size_t m = 0; m < this->saved_state.size(); m++) {
              this->saved_state[m]->copy(state[m]);
            }
          }
        }

        //! @{
        void set_quadrature(shared_ptr<IQuadrature<time>> quadrature)
        {
          this->quadrature = quadrature;
        }

        shared_ptr<const IQuadrature<time>> get_quadrature() const
        {
          return this->quadrature;
        }

        shared_ptr<Encapsulation<time>> get_start_state() const
        {
          return this->start_state;
        }

        const vector<time> get_nodes() const
        {
          return this->quadrature->get_nodes();
        }

        void set_factory(shared_ptr<EncapFactory<time>> factory)
        {
          this->factory = factory;
        }

        virtual shared_ptr<EncapFactory<time>> get_factory() const
        {
          return factory;
        }

        virtual shared_ptr<Encapsulation<time>> get_end_state()
        {
          return this->end_state;
        }
        //! @}

        //! @{
        /**
         * @copydoc ISweeper::advance()
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void advance() override
        {
          throw NotImplementedYet("sweeper");
        }

        /**
         * Re-evaluate function values.
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void reevaluate(bool initial_only=false)
        {
          throw NotImplementedYet("sweeper");
        }

        /**
         * integrates values of right hand side at all time nodes \\( t \\in [0,M-1] \\)
         * simultaneously
         *
         * @param[in] dt width of the time interval to integrate
         * @param[in,out] dst integrated values
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
        {
          UNUSED(dt); UNUSED(dst);
          throw NotImplementedYet("sweeper");
        }
        //! @}

        //! @{
        /**
         * Set residual tolerances for convergence checking.
         */
        void set_residual_tolerances(time abs_residual_tol, time rel_residual_tol, int order=0)
        {
          this->abs_residual_tol = abs_residual_tol;
          this->rel_residual_tol = rel_residual_tol;
          this->residual_norm_order = order;
        }

        /**
         * Compute residual at each SDC node (including FAS corrections).
         */
        virtual void residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
        {
          throw NotImplementedYet("residual");
        }

        /**
         * Return convergence status.
         *
         * This is used by controllers to shortcircuit iterations.
         */
        virtual bool converged() override
        {
          if (this->abs_residual_tol > 0.0 || this->rel_residual_tol > 0.0) {
            if (this->residuals.size() == 0) {
              for (auto x: this->get_nodes()) {
                this->residuals.push_back(this->get_factory()->create(pfasst::encap::solution));
              }
            }
            this->residual(this->get_controller()->get_time_step(), this->residuals);
            vector<time> rnorms;
            for (auto r: this->residuals) {
              rnorms.push_back(r->norm0());
            }
            auto rmax = *std::max_element(rnorms.begin(), rnorms.end());
            if (rmax < this->abs_residual_tol) {
              return true;
            }
            // XXX: check rel norms too
          }
          return false;
        }
        //! @}

        //! @{
        virtual void post(ICommunicator* comm, int tag) override
        {
          this->start_state->post(comm, tag);
        }

        virtual void send(ICommunicator* comm, int tag, bool blocking) override
        {
          this->end_state->send(comm, tag, blocking);
        }

        virtual void recv(ICommunicator* comm, int tag, bool blocking) override
        {
          this->start_state->recv(comm, tag, blocking);
          if (this->quadrature->left_is_node()) {
            this->state[0]->copy(this->start_state);
          }
        }

        virtual void broadcast(ICommunicator* comm) override
        {
          if (comm->rank() == comm->size() - 1) {
            this->start_state->copy(this->end_state);
          }
          this->start_state->broadcast(comm);
        }
        //! @}
    };


    template<typename time>
    EncapSweeper<time>& as_encap_sweeper(shared_ptr<ISweeper<time>> x)
    {
      shared_ptr<EncapSweeper<time>> y = dynamic_pointer_cast<EncapSweeper<time>>(x);
      assert(y);
      return *y.get();
    }


    template<typename time>
    const EncapSweeper<time>& as_encap_sweeper(shared_ptr<const ISweeper<time>> x)
    {
      shared_ptr<const EncapSweeper<time>> y = dynamic_pointer_cast<const EncapSweeper<time>>(x);
      assert(y);
      return *y.get();
    }

  }  // ::pfasst::encap
}  // ::pfasst

#endif
