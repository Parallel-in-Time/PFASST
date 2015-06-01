/**
 * @file pfasst/encap/encap_sweeper.hpp
 * @since v0.1.0
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


    /**
     * Host based encapsulated base sweeper.
     *
     * @tparam time time precision; defaults to pfasst::time_precision
     */
    template<typename time = time_precision>
    class EncapSweeper
      : public ISweeper<time>
    {
      protected:
        //! @{
        //! Quadrature rule used by this sweeper.
        shared_ptr<IQuadrature<time>> quadrature;

        //! Encapsulation data structure factory.
        shared_ptr<EncapFactory<time>> factory;

        //! Separate start state, i.e. initial condition for the sweeper's current time step.
        shared_ptr<Encapsulation<time>> start_state;

        //! Current solution at \\( T_{end} \\).
        shared_ptr<Encapsulation<time>> end_state;

        /**
         * Place for the residuals at the different time nodes.
         *
         * The index of the vector corresponds to the index of the quadrature nodes, i.e.
         * `residuals.size() == quadrature->get_num_nodes()`.
         */
        vector<shared_ptr<Encapsulation<time>>> residuals;
        //! @}

        //! @{
        /**
         * Solution values \\( U \\) at all time nodes of the current iteration.
         *
         * The index of the vector corresponds to the index of the quadrature nodes, i.e.
         * `state.size() == quadrature->get_num_nodes()`.
         */
        vector<shared_ptr<Encapsulation<time>>> state;

        /**
         * Solution values \\( U \\) at all time nodes of the previous iteration.
         *
         * The index of the vector corresponds to the index of the quadrature nodes, i.e.
         * `saved_state.size() == quadrature->get_num_nodes()`.
         */
        vector<shared_ptr<Encapsulation<time>>> saved_state;

        /**
         * FAS corrections \\( \\tau \\) at all time nodes of the current iteration.
         *
         * The index of the vector corresponds to the index of the quadrature nodes, i.e.
         * `fas_corrections.size() == quadrature->get_num_nodes()`.
         */
        vector<shared_ptr<Encapsulation<time>>> fas_corrections;
        //! @}

        //! @{
        //! @todo Write documentation for this member.
        int residual_norm_order;

        /**
         * Tolerance for absolute residual.
         *
         * The absolute residual is the residual between the very first iteration and the current
         * state.
         */
        time abs_residual_tol;

        /**
         * Tolerance for relative residual.
         *
         * The relative residual is the residual between the previous iteration and current state.
         */
        time rel_residual_tol;
        //! @}

        string FORMAT_STR;

      public:
        EncapSweeper();

        //! @{
        /**
         * Retrieve solution values of current iteration at time node index @p m.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const;

        /**
         * Retrieve FAS correction of current iteration at time node index @p m.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const;

        /**
         * Retrieve solution values of previous iteration at time node index @p m.
         *
         * @param[in] m 0-based index of time node
         */
        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const;
        //! @}

        //! @{
        /**
         * @copybrief ISweeper::set_options()
         */
        virtual void set_options() override;

        /**
         * @copybrief ISweeper::setup()
         */
        virtual void setup(bool coarse) override;
        //! @}

        //! @{
        /**
         * @copybrief ISweeper::spread()
         */
        virtual void spread() override;

        /**
         * @copybrief ISweeper::save()
         */
        virtual void save(bool initial_only) override;
        //! @}

        //! @{
        virtual void set_quadrature(shared_ptr<IQuadrature<time>> quadrature);
        virtual shared_ptr<const IQuadrature<time>> get_quadrature() const;

        virtual const vector<time> get_nodes() const;

        virtual void set_factory(shared_ptr<EncapFactory<time>> factory);
        virtual shared_ptr<EncapFactory<time>> get_factory() const;

        virtual shared_ptr<Encapsulation<time>> get_start_state() const;
        virtual shared_ptr<Encapsulation<time>> get_end_state();
        //! @}

        //! @{
        /**
         * @copybrief ISweeper::advance()
         */
        virtual void advance() override;

        /**
         * Re-evaluate function values.
         *
         * @param[in] initial_only whether the right hand side should only be evaluated at the
         *   initial time point
         */
        virtual void reevaluate(bool initial_only = false);

        /**
         * Integrates values of right hand side at all time nodes \\( t \\in [0,M-1] \\)
         * simultaneously.
         *
         * @param[in]     dt  width of the time interval to integrate
         * @param[in,out] dst integrated values
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const;
        //! @}

        //! @{
        /**
         * Set residual tolerances for convergence checking.
         *
         * @param[in] abs_residual_tol tolerance for the absolute residual
         * @param[in] rel_residual_tol tolerance for the relative residual
         * @param[in] order
         */
        void set_residual_tolerances(time abs_residual_tol, time rel_residual_tol, int order = 0);

        /**
         * Compute residual at each SDC node (including FAS corrections).
         *
         * @param[in]     dt  width of the time interval to compute the residual for
         * @param[in,out] dst place to store the residuals at time nodes
         */
        virtual void residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const;

        /**
         * @copybrief ISweeper::converged()
         */
        virtual bool converged() override;
        //! @}

        //! @{
        /**
         * @copybrief ISweeper::post()
         */
        virtual void post(ICommunicator* comm, int tag) override;

        /**
         * @copybrief ISweeper::send()
         */
        virtual void send(ICommunicator* comm, int tag, bool blocking) override;

        /**
         * @copybrief ISweeper::recv()
         */
        virtual void recv(ICommunicator* comm, int tag, bool blocking) override;

        /**
         * @copybrief ISweeper::broadcast()
         */
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
