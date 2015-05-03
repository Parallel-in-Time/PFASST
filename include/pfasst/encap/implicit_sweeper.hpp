#ifndef _PFASST_ENCAP_IMPLICIT_SWEEPER_HPP_
#define _PFASST_ENCAP_IMPLICIT_SWEEPER_HPP_

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
         * Values of the implicit part of the right hand side \\( F_{impl}(t,u) \\) at all time nodes of the current
         * iteration.
         */
        vector<shared_ptr<Encapsulation<time>>> fs_impl;
        //! @}

        Matrix<time> q_tilde;

        void set_end_state();

      public:
        //! @{
        ImplicitSweeper() = default;
        virtual ~ImplicitSweeper() = default;
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
         * @copybrief EncapSweeper::evaluate()
         */
        virtual void reevaluate(bool initial_only) override;

        /**
         * @copybrief EncapSweeper::integrate()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst integrated values; will get zeroed out beforehand
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override;
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
          throw NotImplementedYet("implicit (f_impl_eval)");
        }

        /**
         * Solve \\( U - \\Delta t F_{\\rm impl}(U) = RHS \\) for \\( U \\).
         *
         * During an implicit SDC sweep, the correction equation is evolved using a backward-Euler
         * stepper.  This routine (implemented by the user) performs the solve required to perform
         * one backward-Euler sub-step, and also returns \\( F_{\\rm impl}(U) \\).
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
          throw NotImplementedYet("implicit (impl_solve)");
        }
        //! @}
    };

  }  // ::pfasst::encap
}  // ::pfasst

#include "pfasst/encap/implicit_sweeper_impl.hpp"

#endif
