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
    template<typename time = time_precision>
    class EncapSweeper
      : public ISweeper<time>
    {
      protected:
        //! @{
        quadrature::IQuadrature<time>* quad;
        shared_ptr<EncapFactory<time>> factory;
        shared_ptr<Encapsulation<time>> start_state;
        shared_ptr<Encapsulation<time>> end_state;
        vector<shared_ptr<Encapsulation<time>>> residuals;
        //! @}

        int residual_norm_order;
        time abs_residual_tol, rel_residual_tol;

      public:
        //! @{
        EncapSweeper()
          : quad(nullptr), abs_residual_tol(0.0), rel_residual_tol(0.0)
        {}

        virtual ~EncapSweeper()
        {
          if (this->quad) delete this->quad;
        }
        //! @}

        //! @{
        virtual void spread() override
        {
          for (size_t m = 1; m < this->quad->get_num_nodes(); m++) {
            //            this->get_state(m)->copy(this->start_state);
            this->get_state(m)->copy(this->get_state(0));
          }
        }
        //! @}

        //! @{
        void set_quadrature(quadrature::IQuadrature<time>* quad)
        {
          this->quad = quad;
        }

        const quadrature::IQuadrature<time>* get_quadrature() const
        {
          return this->quad;
        }

        shared_ptr<Encapsulation<time>> get_start_state() const
        {
          return this->start_state;
        }

        const vector<time> get_nodes() const
        {
          return this->quad->get_nodes();
        }

        void set_factory(shared_ptr<EncapFactory<time>> factory)
        {
          this->factory = factory;
        }

        virtual shared_ptr<EncapFactory<time>> get_factory() const
        {
          return factory;
        }

        /**
         * retrieve solution values of current iteration at time node index `m`
         *
         * @param[in] m 0-based index of time node
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const
        {
          UNUSED(m);
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        /**
         * retrieves FAS correction of current iteration at time node index `m`
         *
         * @param[in] m 0-based index of time node
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const
        {
          UNUSED(m);
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        /**
         * retrieves solution values of previous iteration at time node index `m`
         *
         * @param[in] m 0-based index of time node
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const
        {
          UNUSED(m);
          throw NotImplementedYet("sweeper");
          return NULL;
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
         * evaluates the right hand side at given time node
         *
         * This evaluates the right hand side at the given time node with index `m` as returned by
         * pfasst::encap::EncapSweeper::get_nodes:
         *
         * @param[in] m index of the time node to evaluate at
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void evaluate(size_t m)
        {
          UNUSED(m);
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
          this->get_state(0)->post(comm, tag);
        }

        virtual void send(ICommunicator* comm, int tag, bool blocking) override
        {
          this->get_state(this->get_nodes().size() - 1)->send(comm, tag, blocking);
        }

        virtual void recv(ICommunicator* comm, int tag, bool blocking) override
        {
          this->get_state(0)->recv(comm, tag, blocking);
        }

        virtual void broadcast(ICommunicator* comm) override
        {
          if (comm->rank() == comm->size() - 1) {
            this->get_state(0)->copy(this->get_state(this->get_nodes().size() - 1));
          }
          this->get_state(0)->broadcast(comm);
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
