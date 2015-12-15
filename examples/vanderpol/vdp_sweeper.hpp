/**
 * @defgroup VanDerPolFiles Files
 * @ingroup VanDerPol
 *
 * @file examples/vanderpol/vdp_sweeper.hpp
 * @since v0.3.0
 */
#ifndef _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_
#define _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#include <pfasst/logging.hpp>
#include <pfasst/encap/implicit_sweeper.hpp>
#include <pfasst/encap/vector.hpp>


namespace pfasst
{
  namespace examples
  {
    /**
     * @defgroup VanDerPol Van-der-Pol Oscillator
     * @ingroup Examples
     *
     * This directory contains an implementations of solver for the Van-der-Pol oscillator using the
     * PFASST framework.
     *
     * The SDC sweeper is defined in vdp_sweeper.hpp.
     * As this simple equation does not require any special handling, all the SDC magic is derived and
     * taken from pfasst::encap::IMEXSweeper::sweep().
     *
     */
    namespace vdp
    {
      /**
       * Sweeper for the van der Pol oscillator.
       *
       * \\[
       * x' = y\\\\
       * y' = \\nu*(1 - x^2)*y - x
       * \\]
       *
       * Note that an analytical solution is available only for \\( \\nu=0 \\), where the vdP simplifies
       * to the standard linear oscillator. Hence, actual errors are computed only for \\( \\nu=0 \\).
       *
       * @ingroup VanDerPol
       */
      template<typename time = pfasst::time_precision>
      class VdpSweeper
        : public encap::ImplicitSweeper<time>
      {
        private:
          typedef encap::Encapsulation<time> encap_type;

          //! Define a type for a complex PFASST vector encapsulation.
          typedef encap::VectorEncapsulation<double> real_vector_type;

          //! Parameter \\( \\nu \\) in the van-der-Pol oscillator
          double nu;

          //! Starting values
          double x0, y0;

          //! Maximum number of iterations for the nonlinear Newton solver
          size_t newton_maxit;

          //! Tolerance for the nonlinear Newton iteration
          double newton_tol;

          //! Counters for how often `f_impl_eval` and `impl_solve` are called.
          size_t n_f_impl_eval, n_impl_solve, n_newton_iter;

          //! Output file
          fstream output_file;

          //! error: For \\( \\nu=0 \\), vdP reduces to linear oscillator, otherwise no analytic solution is available
          double error;

        public:
          /**
           * generic constructor; initialize all function call counters with zero.
           *
           * Also open file to write solution.
           *
           * @param[in] coefficient nu of van-der-Pol oscillator
           */
          VdpSweeper(double nu, double x0, double y0)
            :   nu(nu)
              , x0(x0)
              , y0(y0)
              , newton_maxit(50)
              , newton_tol(1e-12)
              , n_f_impl_eval(0)
              , n_impl_solve(0)
              , n_newton_iter(0)
          {
            this->output_file.open("./vanderpol.txt", ios_base::out);

            // The file should also contain the initial value, so perform a first write here
            this->output_file << x0 << "    " << y0 << endl;
          }

          /**
           * upon destruction, report final error and number of function calls and close output file.
           */
          virtual ~VdpSweeper()
          {
            ML_LOG(INFO, "Number of implicit evaluations:" << this->n_f_impl_eval);
            ML_LOG(INFO, "Number of implicit solves:     " << this->n_impl_solve);
            this->output_file.close();
          }

          /**
           * computes and prints out error, which is meaningful only for \\( \\nu=0 \\).
           *
           * Also writes the current solution into an ascii file
           *
           * @param[in] t Time
           */
          void echo_error(time t)
          {
            auto& qend = encap::as_vector<double, time>(this->get_end_state());

            real_vector_type qex(qend.size());

            this->exact(qex, t);

            if (this->nu==0)
            {
              double max_err = max(abs(qend[0] - qex[0])/abs(qex[0]) , abs(qend[1]-qex[1])/abs(qex[1]) );
              ML_LOG(INFO, "error:" << max_err);
              this->error = max_err;
            }
            this->output_file << qend[0] << "    " << qend[1] << endl;
          }

          /**
           * returns error, but does not update it!
           */
          double get_errors()
          {
            return this->error;
          }

          /**
           * post prediction step, update error.
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

          // I don't know why Doxygen is eating up so many backslashes here ...
          /**
           * if \\\\( \\nu=0 \\\\), this computes the exact solution at time @p t.
           *
           * For other values of \\\\( \\nu \\\\), it just returns the inital value, because there
           * is no analytical solution to evaluate.
           *
           * @param[in] t Time
           */
          void exact(real_vector_type& q, time t)
          {
            if (this->nu==0)
            {
              /**
               * For \\( \\nu=0 \\) and given initial value \\( x0 \\), \\( y0 \\), the analytic
               * solution reads
               *
               * \\[
               * x(t) =  y0*sin(t) + x0*cos(t)\\\\
               * y(t) = -x0*sin(t) + y0*cos(t)
               * \\]
               *
               * Note that for \\( t=0 \\), we recover \\( x(0) = x0 \\), \\( y(0) = y0 \\).
               */
              q[0] = this->y0*sin(t) + this->x0*cos(t);
              q[1] = -this->x0*sin(t) + this->y0*cos(t);
            }
            else
            {
              // For the nonlinear case, there is no analytic solution for the vdP oscillator, so
              // simply return the initial value
              q[0] = this->x0;
              q[1] = this->y0;
            }
          }

          void exact(shared_ptr<encap_type> q_encap, time t)
          {
            auto& q = encap::as_vector<double, time>(q_encap);
            this->exact(q, t);
          }

          /**
           * Evaluate the implicit part of the right hand side.
           */
          void f_impl_eval(shared_ptr<encap_type> f_encap,
                           shared_ptr<encap_type> q_encap, time t) override
          {
            UNUSED(t);
            auto& f = encap::as_vector<double, time>(f_encap);
            auto& q = encap::as_vector<double, time>(q_encap);

            // x' = y
            f[0] = q[1];
            // y' = nu(1-x^2)*y - x
            f[1] = this->nu*(1.0-q[0]*q[0])*q[1] - q[0];

            this->n_f_impl_eval++;
          }

          // I don't know why Doxygen is eating up so many backslashes here ...
          /**
           * for given \\\\( b \\\\), solve \\\\( \\left( u - \\Delta t f(u) \\right) u = b \\\\) for
           * \\\\( u \\\\) and set @p f_encap to \\\\( f(u) \\\\).
           */
          void impl_solve(shared_ptr<encap_type> f_encap,
                          shared_ptr<encap_type> q_encap, time t, time dt,
                          shared_ptr<encap_type> rhs_encap) override
          {
            UNUSED(t);
            auto& f = encap::as_vector<double, time>(f_encap);
            auto& q = encap::as_vector<double, time>(q_encap);
            auto& rhs = encap::as_vector<double, time>(rhs_encap);

            // I don't know why Doxygen is eating up so many backslashes here ...
            /**
             * Solves the nonlinear equation
             *
             * \\\\[
             * P(q) := q - dt*f(q) - rhs = 0
             * \\\\]
             *
             * using Newton's method.
             * For the vdP oscillator, it is
             *
             * \\\\[
             * P(q) = \\left( x - dt*y - rhs_0 ; y - dt*\\left( \\nu*(1-x^2)*y -x \\right) - rhs_1 \\right)
             * \\\\]
             *
             * The Newton update reads
             *
             * \\\\[
             * q_new = q + inv(P'(q))*(-P(q))
             * \\\\]
             *
             * Note that \\\\( P'(q) = Id - dt*f'(q) \\\\).
             *
             * The inverse of \\\\( P' \\\\) has been computed symbolically, so that \\\\( P'(q)*x \\\\)
             * can be evaluated directly.
             */

            // initialize residual
            double residual = this->newton_tol + 1.0;
            size_t iter     = 0.0;

            // Initial value for q is just rhs: For small dt, P is approximately the identity
            q[0] = rhs[0];
            q[1] = rhs[1];

            // NEWTON ITERATION: q_new = q + inv(J(q))*(-f(q))
            do {
              // Store -P(q) in f0, f1
              double f0 = -( q[0] - dt*q[1] - rhs[0] );
              double f1 = -( q[1] - dt*( this->nu*(1-q[0]*q[0])*q[1]-q[0]) - rhs[1] );

              /**
               * The Jacobian of the right hand side of the van der Pol oscillator for \\( q=[x;y] \\)
               * is
               *
               * \\[
               * Df(q) = \\left( 0  1  ; -2*\\nu*x*y - 1  \\nu*(1-x^2) \\right)
               * \\]
               *
               * Because here we need to invert \\( Id - dt*f(q) = b \\), the corresponding Jacobian
               * is
               *
               * \\[
               * J(q) = \\left( 1 -dt ; dt*2*\\nu*x*y+1  1 + dt*\\nu(x^2-1) \\right)
               * \\]
               *
               * which can be inverted analytically
               *
               * \\[
               * inv(J(q)) = (1/c)*\\left( dt*x^2-dt+1  dt ; -2*dt*\\nu*x*y-dt  1\\right) =: [a dt; b 1.0]
               * \\]
               *
               * with \\( c:=2*\\nu*x*y*dt^2 + dt^2 + dt*x^2 - dt + 1 \\)
               */
              double a = dt * q[0] * q[0] - dt + 1.0;
              double b = -2.0 * dt * this->nu * q[0] * q[1] - dt;
              double c = 2.0 * this->nu * q[0] * q[1] * dt * dt + dt * dt + dt * q[0] * q[0] - dt + 1.0;

              // Compute inv(J(q))*(-f(q)) and store it in f
              f[0] = (1.0/c)*( a*f0 + dt*f1 );
              f[1] = (1.0/c)*( b*f0 + f1 );

              // Perform Newton update q_new = q + inv(J(q))*(-f(q))
              q[0] += f[0];
              q[1] += f[1];

              // Compute relative residual; note that f = q_new - q is the update
              residual = fmax( abs(f[0]), abs(f[1]) ) / fmax( abs(q[0]), abs(q[1]) );

              iter++;
              this->n_newton_iter++;

            } while ( (iter<this->newton_maxit) && (residual>this->newton_tol) );

            // Say something, only if the residual tolerance has not been reach in the maximum
            // number of iterations
            if (residual >=this->newton_tol)
            {
              cout << "Newton failed to converge: res = " << scientific << residual << " -- n_iter = " << iter << " of maxit = " << this->newton_maxit << endl;
            }

            // Set f to f(q)
            f[0] = q[1];
            f[1] = this->nu*(1.0-q[0]*q[0])*q[1] - q[0];

            this->n_impl_solve++;
          }
      };
    }  // ::pfasst::examples::vdp
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_
