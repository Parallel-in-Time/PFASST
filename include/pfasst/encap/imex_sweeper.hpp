
#ifndef _PFASST_ENCAP_IMEX_SWEEPER_HPP_
#define _PFASST_ENCAP_IMEX_SWEEPER_HPP_

#include <cstdlib>
#include <cassert>
#include <vector>
#include <memory>

#include <boost/numeric/ublas/matrix.hpp>

#include "../globals.hpp"
#include "encapsulation.hpp"
#include "encap_sweeper.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {
    using pfasst::encap::Encapsulation;

    /**
     * interface for an semi-implicit sweeper
     *
     * Given an ODE \\( \\frac{\\partial}{\\partial t}u(t) = F(t,u) \\) where the function of the
     * right hand side \\( F(t,u) \\) can be split into a non-stiff and a stiff part.
     * To reduce complexity and computational efford one would want to solve the non-stiff part
     * explicitly and the stiff part implicitly.
     * Therefore, we define the splitting \\( F(t,u) = F_{expl}(t,u) + F_{impl}(t,u) \\).
     *
     * This sweeper provides an interface for such ODEs were the implicit part can be computed by
     * an external implicit solver without actually evaluating \\( F_{impl}(t,u) \\), which is
     * possibly very expensive.
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
        vector<shared_ptr<Encapsulation<time>>> us;

        /**
         * solution values \\( u(t) \\) at all time nodes \\( t \\in [0, M-1] \\) of the
         * previous iteration
         */
        vector<shared_ptr<Encapsulation<time>>> previous_us;

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

        /**
         * values of the explicit part of the right hand side \\( F_{impl}(t,u) \\) at all time
         * nodes \\( t \\in [0, M-1] \\) of the current iteration
         */
        vector<shared_ptr<Encapsulation<time>>> fs_impl;
        //! @}

        //! @{
        /**
         * quadrature matrix containing weights for node-to-node integration
         */
        matrix<time> s_mat;

        /**
         * quadrature matrix containing weights for node-to-node integration of explicit part
         *
         * @see IMEXSweeper::setup(bool) for a short description
         */
        matrix<time> s_mat_expl;

        /**
         * quadrature matrix containing weights for node-to-node integration of implicit part
         *
         * @see IMEXSweeper::setup(bool) for a short description
         */
        matrix<time> s_mat_impl;
        //! @}

      public:
        //! @{
        virtual ~IMEXSweeper()
        {}
        //! @}

        //! @{
        virtual void set_state(shared_ptr<const Encapsulation<time>> u0, size_t m) override
        {
          this->us[m]->copy(u0);
        }

        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const override
        {
          return this->us[m];
        }

        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const override
        {
          return this->fas_corrections[m];
        }

        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const override
        {
          return this->previous_us[m];
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
          auto nodes = this->get_nodes();
          assert(nodes.size() >= 1);

          this->s_mat = compute_quadrature(nodes, nodes, 's');

          this->s_mat_expl = this->s_mat;
          this->s_mat_impl = this->s_mat;
          for (size_t m = 0; m < nodes.size() - 1; m++) {
            time ds = nodes[m + 1] - nodes[m];
            this->s_mat_expl(m, m)     -= ds;
            this->s_mat_impl(m, m + 1) -= ds;
          }

          for (size_t m = 0; m < nodes.size(); m++) {
            this->us.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->previous_us.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
            this->fs_expl.push_back(this->get_factory()->create(pfasst::encap::function));
            this->fs_impl.push_back(this->get_factory()->create(pfasst::encap::function));
          }

          for (size_t m = 0; m < nodes.size() - 1; m++) {
            this->s_integrals.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->fas_corrections.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
          }
        }

        virtual void predict(bool initial)
        {
          const auto   nodes  = this->get_nodes();
          const size_t nnodes = nodes.size();
          assert(nnodes >= 1);

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          if (initial) {
            this->f_expl_eval(this->fs_expl[0], this->us[0], t);
            this->f_impl_eval(this->fs_impl[0], this->us[0], t);
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          for (size_t m = 0; m < nnodes - 1; m++) {
            time ds = dt * (nodes[m + 1] - nodes[m]);
            rhs->copy(this->us[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            this->impl_solve(this->fs_impl[m + 1], this->us[m + 1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m + 1], this->us[m + 1], t + ds);

            t += ds;
          }
        }

        virtual void sweep()
        {
          const auto   nodes  = this->get_nodes();
          const size_t nnodes = nodes.size();
          assert(nnodes >= 1);

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          // integrate
          this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->s_mat_expl, this->fs_expl, true);
          this->s_integrals[0]->mat_apply(this->s_integrals, dt, this->s_mat_impl, this->fs_impl, false);
          if (this->fas_corrections.size() > 0) {
            for (size_t m = 0; m < nnodes - 1; m++) {
              this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m]);
            }
          }

          // sweep
          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          for (size_t m = 0; m < nnodes - 1; m++) {
            time ds = dt * (nodes[m + 1] - nodes[m]);

            rhs->copy(this->us[m]);
            rhs->saxpy(ds, this->fs_expl[m]);
            rhs->saxpy(1.0, this->s_integrals[m]);
            this->impl_solve(this->fs_impl[m + 1], this->us[m + 1], t, ds, rhs);
            this->f_expl_eval(this->fs_expl[m + 1], this->us[m + 1], t + ds);

            t += ds;
          }
        }

        virtual void advance() override
        {
          this->us[0]->copy(this->us.back());
          this->fs_expl[0]->copy(this->fs_expl.back());
          this->fs_impl[0]->copy(this->fs_impl.back());
        }

        virtual void save(bool initial_only) override
        {
          if (initial_only) {
            this->previous_us[0]->copy(us[0]);
          } else {
            for (size_t m = 0; m < this->previous_us.size(); m++) {
              this->previous_us[m]->copy(us[m]);
            }
          }
        }

        virtual void evaluate(size_t m) override
        {
          time t = this->get_nodes()[m]; // XXX
          this->f_expl_eval(this->fs_expl[m], this->us[m], t);
          this->f_impl_eval(this->fs_impl[m], this->us[m], t);
        }

        /**
         * @copybrief EncapSweeper::integrate()
         *
         * @param[in] dt width of time interval to integrate over
         * @param[in,out] dst integrated values; will get zeroed out beforehand
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override
        {
          dst[0]->mat_apply(dst, dt, this->s_mat, this->fs_expl, true);
          dst[0]->mat_apply(dst, dt, this->s_mat, this->fs_impl, false);
        }
        //! @}

        //! @{
        /**
         * Evaluates the explicit part of the right hand side at the given time.
         *
         * @param[in,out] f_expl_encap Encapsulation to store the evaluated right hand side
         * @param[in] u_encap Encapsulation storing the solution values to use for computing the
         *     explicit part of the right hand side
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
         * Evaluates the implicit part of the right hand side at the given time.
         *
         * This is typically called to compute the implicit part of the right hand side at the first
         * collocation node, and on all nodes after restriction or interpolation.
         *
         * @param[in,out] f_impl_encap Encapsulation to store the evaluated right hand side
         * @param[in] u_encap Encapsulation storing the solution values to use for computing the
         *     implicit part of the right hand side
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
         * Solves \\( U - \\Delta t f_{\\rm impl}(U) = RHS \\) for \\(U\\).
         *
         * During an IMEX SDC sweep, the correction equation is evolved using a forward-Euler
         * stepper for the explicit piece, and a backward-Euler stepper for the implicit piece.
         * This routine (implemented by the user) performs the solve required to perform one
         * backward-Euler sub-step, and also returns \\(f_{\\rm impl}(U)\\).
         *
         * @param[in,out] f_encap Encapsulation to store the evaluated right hand side
         * @param[in,out] u_encap Encapsulation to store the solution of the backward-Euler sub-step
         * @param[in] t time point (of \\(RHS\\))
         * @param[in] dt sub-step size to the previous time point (\\(\\Delta t \\))
         * @param[in] rhs_encap Encapsulation storing \\(RHS\\)
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void impl_solve(shared_ptr<Encapsulation<time>> f_encap,
                                shared_ptr<Encapsulation<time>> u_encap,
                                time t, time dt,
                                shared_ptr<Encapsulation<time>> rhs_encap)
        {
          UNUSED(f_encap); UNUSED(u_encap); UNUSED(t); UNUSED(dt); UNUSED(rhs_encap);
          throw NotImplementedYet("imex (f2comp)");
        }
        //! @}
    };

  }  // ::pfasst::encap
}  // ::pfasst

#endif
