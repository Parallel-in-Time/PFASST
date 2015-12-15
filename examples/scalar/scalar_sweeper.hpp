/**
 * @defgroup ScalarFiles Files
 * @ingroup Scalar
 *
 * This directory contains a simple implementations of scalar ODE solver using the PFASST framework.
 *
 * The SDC sweeper is defined in `scalar_sweeper.hpp`.
 * As this simple equation does not require any special handling, all the SDC magic is derived and
 * taken from pfasst::encap::IMEXSweeper::sweep().
 *
 * @file examples/scalar/scalar_sweeper.hpp
 * @since v0.2.0
 */
#ifndef _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_
#define _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_

#include <complex>
#include <vector>
using namespace std;

#include <pfasst/logging.hpp>
#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/vector.hpp>


namespace pfasst
{
  namespace examples
  {
    /**
     * @defgroup Scalar Scalar
     * @ingroup Examples
     *
     * This directory contains a simple implementations of scalar ODE solver using the PFASST framework.
     *
     * The SDC sweeper is defined in `scalar_sweeper.hpp`.
     * As this simple equation does not require any special handling, all the SDC magic is derived and
     * taken from pfasst::encap::IMEXSweeper::sweep().
     */
    namespace scalar
    {
      /**
       * Sweeper for scalar test equation.
       *
       * \\[ u' = \\lambda*u \\quad\\text{ , } u(0) = u_0 \\]
       * with complex lambda using an IMEX scheme.
       *
       * @ingroup Scalar
       */
      template<typename time = pfasst::time_precision>
      class ScalarSweeper
        : public encap::IMEXSweeper<time>
      {
        private:
          typedef encap::Encapsulation<time> encap_type;

          //! define a type for a complex PFASST vector encapsulation.
          typedef encap::VectorEncapsulation<complex<double>> complex_vector_type;

           //! parameter lambda and initial value \\( u_0 \\)
          complex<double> lambda, u0;

          //! the complex unit \\( i = \\sqrt{-1} \\)
          const complex<double> i_complex = complex<double>(0, 1);

           //! error at the final time. For the scalar example, an analytical solution is known.
          double error;

           //! counters for how often `f_expl_eval`, `f_impl_eval` and `impl_solve` are called.
          size_t n_f_expl_eval, n_f_impl_eval, n_impl_solve;

        public:
          /**
           * generic constructor; initialize all function call counters with zero.
           *
           * @param[in] lambda coefficient in test equation
           * @param[in] u0 initial value at \\( t=0 \\)
           */
          ScalarSweeper(const complex<double>& lambda, const complex<double>& u0)
            :   lambda(lambda)
              , u0(u0)
              , error(0.0)
              , n_f_expl_eval(0)
              , n_f_impl_eval(0)
              , n_impl_solve(0)
          {}

          /**
           * upon destruction, report final error and number of function calls
           */
          virtual ~ScalarSweeper()
          {
            ML_LOG(INFO, "Final error:                   " << this->error);
            ML_LOG(INFO, "Number of explicit evaluations:" << this->n_f_expl_eval);
            ML_LOG(INFO, "Number of implicit evaluations:" << this->n_f_impl_eval);
            ML_LOG(INFO, "Number of implicit solves:     " << this->n_impl_solve);
          }

          /**
           * compute error between last state and exact solution at time tand print it to cout
           *
           * @param[in] t Time
           */
          void echo_error(time t)
          {
            auto& qend = encap::as_vector<complex<double>, time>(this->get_end_state());

            complex_vector_type qex(qend.size());

            this->exact(qex, t);
            double max_err = abs(qend[0] - qex[0]) / abs(qex[0]);
            ML_LOG(INFO, "err:" << max_err);
            this->error = max_err;
          }

          /**
           * returns error, but does not update it!
           */
          double get_errors()
          {
            return this->error;
          }

          /**
           * post prediction step, update of error.
           */
          void post_predict() override
          {
            time t  = this->get_controller()->get_time();
            time dt = this->get_controller()->get_step_size();
            this->echo_error(t + dt);
          }

          /**
           * post sweep, update error.
           */
          void post_sweep() override
          {
            time t  = this->get_controller()->get_time();
            time dt = this->get_controller()->get_step_size();
            this->echo_error(t + dt);
          }

          /**
           * computes the exact solution \\( u_0 \\exp \\left( \\lambda*t \\right) \\) at a given time t.
           *
           * @param[in] t Time
           */
          void exact(complex_vector_type& q, time t)
          {
            q[0] = this->u0 * exp(this->lambda * t);
          }

          void exact(shared_ptr<encap_type> q_encap, time t)
          {
            auto& q = encap::as_vector<complex<double>, time>(q_encap);
            this->exact(q, t);
          }

          /**
           * evaluate the explicit part of the right hand side: Multiply with \\( \\text{imag}(\\lambda) \\)
           */
          void f_expl_eval(shared_ptr<encap_type> f_encap,
                           shared_ptr<encap_type> q_encap, time t) override
          {
            UNUSED(t);
            auto& f = encap::as_vector<complex<double>, time>(f_encap);
            auto& q = encap::as_vector<complex<double>, time>(q_encap);

            // f_expl = multiply with imaginary part of lambda
            f[0] = this->i_complex * imag(this->lambda) * q[0];

            this->n_f_expl_eval++;
          }

          /**
           * evaluate the implicit part of the right hand side: Multiply with \\( \\text{real}(\\lambda) \\)
           */
          void f_impl_eval(shared_ptr<encap_type> f_encap,
                           shared_ptr<encap_type> q_encap, time t) override
          {
            UNUSED(t);
            auto& f = encap::as_vector<complex<double>, time>(f_encap);
            auto& q = encap::as_vector<complex<double>, time>(q_encap);

            // f_impl = multiply with real part of lambda
            f[0] = real(this->lambda) * q[0];

            this->n_f_impl_eval++;
          }

          /**
           * for given \\( b \\), solve
           * \\( \\left( \\mathbb{I}_d - \\Delta t \\text{real}(\\lambda) \\right) u = b \\)
           * for \\( u \\) and set f_encap to \\( \\text{real}(\\lambda) u \\)
           */
          void impl_solve(shared_ptr<encap_type> f_encap,
                          shared_ptr<encap_type> q_encap, time t, time dt,
                          shared_ptr<encap_type> rhs_encap) override
          {
            UNUSED(t);
            auto& f = encap::as_vector<complex<double>, time>(f_encap);
            auto& q = encap::as_vector<complex<double>, time>(q_encap);
            auto& rhs = encap::as_vector<complex<double>, time>(rhs_encap);

            // invert f_impl = multiply with inverse of real part of lambda
            double inv = 1.0 / (1.0 - double(dt) * real(this->lambda));
            q[0] = inv * rhs[0];

            // compute f_impl_eval of q[0], i.e. multiply with real(lambda)
            f[0] = real(this->lambda) * q[0];

            this->n_impl_solve++;
          }
      };
    }  // ::pfasst::examples::scalar
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_
