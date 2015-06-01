#ifndef _PFASST_ENCAP_IMEX_SWEEPER_HPP_
#define _PFASST_ENCAP_IMEX_SWEEPER_HPP_

#include <memory>
#include <vector>

#include "pfasst/encap/encapsulation.hpp"
#include "pfasst/encap/encap_sweeper.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {
    using pfasst::encap::Encapsulation;

    /**
     * Semi-implicit IMEX sweeper.
     *
     * This IMEX sweeper is for ODEs of the form
     * \\( \\dot{U} = F_{\\rm expl}(t,U) + F_{\\rm impl}(t, U) \\).
     * To reduce complexity and computational effort the non-stiff part is treated explicitly and
     * the stiff part implicitly.
     *
     * This sweeper requires three interfaces to be implemented: two routines to evaluate the
     * explicit \\( F_{\\rm expl} \\) and implicit \\( F_{\\rm impl} \\) pieces for a given state,
     * and one routine that solves (perhaps with an external solver) the backward-Euler equation
     * \\( U^{n+1} - \\Delta t F_{\\rm impl}(U^{n+1}) = RHS \\) for \\( U^{n+1} \\).
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
        virtual void integrate_end_state(time dt);

      public:
        //! @{
        IMEXSweeper() = default;
        virtual ~IMEXSweeper() = default;
        //! @}

        //! @{
        /**
         * @copydoc ISweeper::setup(bool)
         */
        virtual void setup(bool coarse) override;

        /**
         * Compute low-order provisional solution.
         *
         * This performs forward/backward Euler steps across the nodes to compute a low-order provisional solution.
         *
         * @param[in] initial if `true` the explicit and implicit part of the right hand side of the
         *     ODE get evaluated with the initial value
         */
        virtual void predict(bool initial) override;

        /**
         * Perform one SDC sweep/iteration.
         *
         * This computes a high-order solution from the previous iteration's function values and
         * corrects it using forward/backward Euler steps across the nodes.
         */
        virtual void sweep() override;

        /**
         * Advance the end solution to start solution.
         */
        virtual void advance() override;

        /**
         * @copybrief EncapSweeper::reevaluate()
         */
        virtual void reevaluate(bool initial_only) override;

        /**
         * @copybrief EncapSweeper::integrate()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst integrated values; will get zeroed out beforehand
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override;

        /**
         * @copybrief EncapSweeper::residual()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst residuals
         */
        virtual void residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override;
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
                                 time t);

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
                                 time t);

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
                                shared_ptr<Encapsulation<time>> rhs_encap);
        //! @}

      private:
        virtual void predict_with_left(bool initial);
        virtual void predict_without_left(bool initial);
        virtual void sweep_with_left();
        virtual void sweep_without_left();
    };
  }  // ::pfasst::encap
}  // ::pfasst

#include "pfasst/encap/imex_sweeper_impl.hpp"

#endif
