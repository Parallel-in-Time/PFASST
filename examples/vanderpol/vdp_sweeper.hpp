#ifndef _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_
#define _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_

#include <vector>

#include <pfasst/logging.hpp>
#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/vector.hpp>

#include <iostream>
#include <fstream>

using namespace std;

/**
 * Sweeper for the van der Pol oscillator
 *
 * x' = y
 * y' = nu*(1 - x^2)*y - x
 *
 * Although derived from the IMEX sweeper, the explicit part does not do anything
 * and a fully implicit Euler is used as a base method.
 *
 * Note that an analytical solution is available only for nu=0, where the vdP simplifies
 * to the standard linear oscillator. Hence, actual errors are computed only for nu=0
 */
template<typename time = pfasst::time_precision>
class VdpSweeper
  : public pfasst::encap::IMEXSweeper<time>
{
  private:
    typedef pfasst::encap::Encapsulation<time> encap_type;

    //! Define a type for a complex PFASST vector encapsulation.
    typedef pfasst::encap::VectorEncapsulation<double> real_vector_type;

    //! Parameter nu in the van-der-Pol oscillator
    double nu;

    //! Starting values
    double x0, y0;

    //! Maximum number of iterations for the nonlinear Newton solver
    size_t newton_maxit;

    //! Tolerance for the nonlinear Newton iteration
    double newton_tol;

    //! Counters for how often `f_expl_eval`, `f_impl_eval` and `impl_solve` are called.
    size_t n_f_expl_eval, n_f_impl_eval, n_impl_solve, n_newton_iter;

    //! Output file
    fstream output_file;

    //! error: For nu=0, vdP reduces to linear oscillator, otherwise no analytic solution is available
    double error;

  public:
    /**
     * Generic constructor; initialize all function call counters with zero.
     * Also open file to write solution.
     * @param[in] coefficient nu of van-der-Pol oscillator
     */
    VdpSweeper(double nu, double x0, double y0)
      :   nu(nu)
        , x0(x0)
        , y0(y0)
        , newton_maxit(50)
        , newton_tol(1e-12)
        , n_f_expl_eval(0)
        , n_f_impl_eval(0)
        , n_impl_solve(0)
        , n_newton_iter(0)
    {
      this->output_file.open("./vanderpol.txt", ios_base::out);

      // The file should also contain the initial value, so perform a first write here
      this->output_file << x0 << "    " << y0 << endl;
    }

    /**
     * Upon destruction, report final error and number of function calls and close output file.
     */
    virtual ~VdpSweeper()
    {
      LOG(INFO) << "Number of explicit evaluations:" << this->n_f_expl_eval;
      LOG(INFO) << "Number of implicit evaluations:" << this->n_f_impl_eval;
      LOG(INFO) << "Number of implicit solves:     " << this->n_impl_solve;
      this->output_file.close();
    }

    /**
     * Computes and prints out error, which is meaningful only for nu=0.
     * Also writes the current solution into an ascii file
     *
     * @param[in] t Time
     */
    void echo_error(time t)
    {
      auto& qend = pfasst::encap::as_vector<double, time>(this->get_end_state());

      real_vector_type qex(qend.size());

      this->exact(qex, t);

      if (this->nu==0)
      {
        double max_err = max(abs(qend[0] - qex[0])/abs(qex[0]) , abs(qend[1]-qex[1])/abs(qex[1]) );
        LOG(INFO) << "error:" << max_err;
        this->error = max_err;
      }
      this->output_file << qend[0] << "    " << qend[1] << endl;
    }

    /**
     * Returns error, but does not update it!
     */
    double get_errors()
    {
      return this->error;
    }

    /**
     * Prediction step and update of error. Uses predictor as provided by IMEXSweeper.
     */
    void predict(bool initial) override
    {
      pfasst::encap::IMEXSweeper<time>::predict(initial);
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_time_step();
      this->echo_error(t + dt);
    }

    /**
     * Perform a sweep and update error. Uses sweep as provided by IMEXSweeper.
     */
    void sweep() override
    {
      pfasst::encap::IMEXSweeper<time>::sweep();
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_time_step();
      this->echo_error(t + dt);
    }

    /**
     * If nu=0, this computes the exact solution at time t. For other values of nu, it just returns
     * the inital value, because there is no analytical solution to evaluate.
     *
     * @param[in] t Time
     */
    void exact(real_vector_type& q, time t)
    {
      if (this->nu==0)
      {
        /**
         * For nu=0 and given initial value x0, y0, the analytic solution reads
         * x(t) =  y0*sin(t) + x0*cos(t)
         * y(t) = -x0*sin(t) + y0*cos(t)
         *
         * Note that for t=0, we recover x(0) = x0, y(0) = y0.
         */
        q[0] = this->y0*sin(t) + this->x0*cos(t);
        q[1] = -this->x0*sin(t) + this->y0*cos(t);
      }
      else
      {
        // For the nonlinear case, there is no analytic solution for the vdP oscillator, so simply return the initial value
        q[0] = this->x0;
        q[1] = this->y0;
      }
    }

    void exact(shared_ptr<encap_type> q_encap, time t)
    {
      auto& q = pfasst::encap::as_vector<double, time>(q_encap);
      this->exact(q, t);
    }

    /**
     * Evaluate the explicit part of the right hand side.
     * Because we use a fully implicit Euler for the vdP oscillator, this returns zeros.
     */
    void f_expl_eval(shared_ptr<encap_type> f_encap,
                     shared_ptr<encap_type> q_encap, time t) override
    {
      UNUSED(t);
      UNUSED(q_encap);
      auto& f = pfasst::encap::as_vector<double, time>(f_encap);
      //auto& q = pfasst::encap::as_vector<double, time>(q_encap);

      // We use a fully implicit sweeper, so f_expl_eval returns zero
      f[0] = 0.0;
      f[1] = 0.0;
      this->n_f_expl_eval++;
    }

    /**
     * Evaluate the implicit part of the right hand side.
     * Because we use a fully implicit Euler for the vdP oscillator, here the full rhs of vdP is computed
     */
    void f_impl_eval(shared_ptr<encap_type> f_encap,
                     shared_ptr<encap_type> q_encap, time t) override
    {
      UNUSED(t);
      auto& f = pfasst::encap::as_vector<double, time>(f_encap);
      auto& q = pfasst::encap::as_vector<double, time>(q_encap);

      // x' = y
      f[0] = q[1];
      // y' = nu(1-x^2)*y - x
      f[1] = this->nu*(1.0-q[0]*q[0])*q[1] - q[0];

      this->n_f_impl_eval++;
    }

    /**
     * For given \\( b \\), solve
     * \\( \\left( u - \\Delta t f(u) \\right) u = b \\)
     * for \\( u \\) and set f_encap to \\( f(u) \\)
     */
    void impl_solve(shared_ptr<encap_type> f_encap,
                    shared_ptr<encap_type> q_encap, time t, time dt,
                    shared_ptr<encap_type> rhs_encap) override
    {
      UNUSED(t);
      auto& f = pfasst::encap::as_vector<double, time>(f_encap);
      auto& q = pfasst::encap::as_vector<double, time>(q_encap);
      auto& rhs = pfasst::encap::as_vector<double, time>(rhs_encap);

      /**
       * Solves the nonlinear equation
       *
       * P(q) := q - dt*f(q) - rhs = 0
       *
       * using Newton's method. For the vdP oscillator, it is
       *
       * P(q) = [ x - dt*y - rhs[0] ; y - dt*[ nu*(1-x^2)*y -x ] - rhs[1] ]
       *
       * The Newton update reads
       *
       * q_new = q + inv(P'(q))*(-P(q))
       *
       * Note that P'(q) = Id - dt*f'(q).
       *
       * The inverse of P' has been computed symbolically, so that P'(q)*x can be evaluated directly.
       */

      // initialize residual
      double residual = this->newton_tol + 1.0;
      size_t iter     = 0.0;

      // initialize temporary variables for use in loop
      double a, b, c, f0, f1;

      // Initial value for q is just rhs: For small dt, P is approximately the identity
      q[0] = rhs[0];
      q[1] = rhs[1];

      // NEWTON ITERATION: q_new = q + inv(J(q))*(-f(q))
      do {
        // Store -P(q) in f0, f1
        f0 = -( q[0] - dt*q[1] - rhs[0] );
        f1 = -( q[1] - dt*( this->nu*(1-q[0]*q[0])*q[1]-q[0]) - rhs[1] );

        /**
         * The Jacobian of the right hand side of the van der Pol oscillator for q=[x;y] is
         *
         * Df(q) = [ 0  1  ; -2*nu*x*y - 1  nu*(1-x^2) ]
         *
         * Because here we need to invert Id - dt*f(q) = b, the corresponding Jacobian is
         *
         * J(q) = [1 -dt ; dt*2*nu*x*y+1  1 + dt*nu(x^2-1) ]
         *
         * which can be inverted analytically
         *
         * inv(J(q)) = (1/c)*[ dt*x^2-dt+1  dt ; -2*dt*nu*x*y-dt  1] =: [a dt; b 1.0]
         *
         * with c:=2*nu*x*y*dt^2 + dt^2 + dt*x^2 - dt + 1
         */
        a = dt*q[0]*q[0]-dt+1.0;
        b = -2.0*dt*this->nu*q[0]*q[1]-dt;
        c = 2.0*this->nu*q[0]*q[1]*dt*dt + dt*dt + dt*q[0]*q[0] - dt + 1.0;

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

      // Say something, only if the residual tolerance has not been reach in the maximum number of iterations
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

#endif  // _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_
