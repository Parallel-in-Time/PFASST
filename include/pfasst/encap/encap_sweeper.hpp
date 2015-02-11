/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAP_ENCAP_SWEEPER_HPP_
#define _PFASST_ENCAP_ENCAP_SWEEPER_HPP_

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/interfaces.hpp"
#include "pfasst/quadrature.hpp"
#include "pfasst/encap/encapsulation.hpp"


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
        EncapSweeper();

        //! @{
        /**
         * Retrieve solution values of current iteration at time node index `m`.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const;

        /**
         * Retrieve FAS correction of current iteration at time node index `m`.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const;

        /**
         * Retrieve solution values of previous iteration at time node index `m`.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const;
        //! @}

        //! @{
        virtual void set_options() override;
        virtual void setup(bool coarse) override;
        //! @}

        //! @{
        virtual void spread() override;

        /**
         * Save current solution states.
         */
        virtual void save(bool initial_only) override;
        //! @}

        //! @{
        virtual void set_quadrature(shared_ptr<IQuadrature<time>> quadrature);
        virtual shared_ptr<const IQuadrature<time>> get_quadrature() const;
        virtual shared_ptr<Encapsulation<time>> get_start_state() const;
        virtual const vector<time> get_nodes() const;
        virtual void set_factory(shared_ptr<EncapFactory<time>> factory);
        virtual shared_ptr<EncapFactory<time>> get_factory() const;
        virtual shared_ptr<Encapsulation<time>> get_end_state();
        //! @}

        //! @{
        /**
         * @copydoc ISweeper::advance()
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void advance() override;

        /**
         * Re-evaluate function values.
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void reevaluate(bool initial_only = false);

        /**
         * integrates values of right hand side at all time nodes \\( t \\in [0,M-1] \\)
         * simultaneously
         *
         * @param[in] dt width of the time interval to integrate
         * @param[in,out] dst integrated values
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const;
        //! @}

        //! @{
        /**
         * Set residual tolerances for convergence checking.
         */
        void set_residual_tolerances(time abs_residual_tol, time rel_residual_tol, int order = 0);

        /**
         * Compute residual at each SDC node (including FAS corrections).
         */
        virtual void residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const;

        /**
         * Return convergence status.
         *
         * This is used by controllers to shortcircuit iterations.
         */
        virtual bool converged() override;
        //! @}

        //! @{
        virtual void post(ICommunicator* comm, int tag) override;
        virtual void send(ICommunicator* comm, int tag, bool blocking) override;
        virtual void recv(ICommunicator* comm, int tag, bool blocking) override;
        virtual void broadcast(ICommunicator* comm) override;
        //! @}
    };


    template<typename time>
    EncapSweeper<time>& as_encap_sweeper(shared_ptr<ISweeper<time>> x);

    template<typename time>
    const EncapSweeper<time>& as_encap_sweeper(shared_ptr<const ISweeper<time>> x);
  }  // ::pfasst::encap
}  // ::pfasst

#include "pfasst/encap/encap_sweeper_impl.hpp"

#endif
