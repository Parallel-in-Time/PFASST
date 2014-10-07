#ifndef _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_
#define _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_

#include <vector>

#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/vector.hpp>

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
 
    //! Maximum number of iterations for the nonlinear Newton solver
    size_t newton_maxit;
    
    //! Tolerance for the nonlinear Newton iteration
    double newton_tol;

     //! Counters for how often `f_expl_eval`, `f_impl_eval` and `impl_solve` are called.
    size_t n_f_expl_eval, n_f_impl_eval, n_impl_solve, n_newton_iter;

  public:
    /**
     * Generic constructor; initialize all function call counters with zero.
     * @param[in] coefficient nu of van-der-Pol oscillator
     */ 
    VdpSweeper(double nu)
      : nu(nu)
        , newton_maxit(25)
        , newton_tol(1e-12)
        , n_f_expl_eval(0)
        , n_f_impl_eval(0)
        , n_impl_solve(0)
        , n_newton_iter(0)
    {}

    /**
     * Upon destruction, report final error and number of function calls
     */
    virtual ~VdpSweeper()
    {
      cout << "Number of explicit evaluations: " << this->n_f_expl_eval << endl;
      cout << "Number of implicit evaluations: " << this->n_f_impl_eval << endl;
      cout << "Number of implicit solves:      " << this->n_impl_solve << endl;
    }

    /**
     * Compute error between last state and exact solution at time tand print it to cout
     *
     * @param[in] t Time
     */
    void echo_error(time t)
    {
      auto& qend = pfasst::encap::as_vector<double, time>(this->get_end_state());

      real_vector_type qex(qend.size());

      this->exact(qex, t);
      //double max_err = abs(qend[0] - qex[0]) / abs(qex[0]);
      cout << "x: " << qend[0] << " -- y: " << qend[1] << endl;

    }

    /**
     * Returns error, but does not update it!
     */
    double get_errors()
    {
      return -1.0;
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
     * Computes the exact solution \\( u_0 \\exp \\left( \\lambda*t \\right) \\) at a given time t.
     *
     * @param[in] t Time
     */
    void exact(real_vector_type& q, time t)
    {
      // TODO: There is no analytic solution for the vdP oscillator, so simply return the initial value
      UNUSED(t);
      q[0] = 2.0;
      q[1] = 0.0;
    }

    void exact(shared_ptr<encap_type> q_encap, time t)
    {
      auto& q = pfasst::encap::as_vector<double, time>(q_encap);
      this->exact(q, t);
    }

    /**
     * Evaluate the explicit part of the right hand side.
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
     * \\( \\left( \\mathbb{I}_d - \\Delta t \\text{real}(\\lambda) \\right) u = b \\) 
     * for \\( u \\) and set f_encap to \\( \\text{real}(\\lambda) u \\)
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
        do
        {
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
            
            // Compute norm of update q_new - q = inv(J(q))*(-f(q)) in maximum norm
            residual = fmax( abs(f[0]), abs(f[1]) );
            
            //std::cout << "Newton res: " << scientific << residual << std::endl;

            iter++;
            this->n_newton_iter++;
        
        } while ( (iter<this->newton_maxit) && (residual>this->newton_tol) );
      
        std::cout << "Newton: res = " << scientific << residual << " -- n_iter = " << iter << std::endl;
        
        // Set f to f(q)
        f[0] = q[1];
        f[1] = this->nu*(1.0-q[0]*q[0])*q[1] - q[0];
        
      this->n_impl_solve++;
    }
  
};

#endif  // _EXAMPLES__VANDERPOL__VDP_SWEEPER_HPP_
