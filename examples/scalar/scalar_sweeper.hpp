#ifndef _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_
#define _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_

#include <complex>
#include <vector>

#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/vector.hpp>

using namespace std;

/**
 * Sweeper for scalar test equation
 *
 * \\( u' = \\lambda*u \\quad\\text{ , } u(0) = u_0 \\)
 *
 * with complex lambda using an IMEX scheme. Derived from the generic imex_sweeper.
 *
 */
template<typename time = pfasst::time_precision>
class ScalarSweeper
  : public pfasst::encap::IMEXSweeper<time>
{
  private:
    typedef pfasst::encap::Encapsulation<time> encap_type;

    //! Define a type for a complex PFASST vector encapsulation.
    typedef pfasst::encap::VectorEncapsulation<complex<double>> complex_vector_type;

     //! Parameter lambda and initial value \\( u_0 \\)
    complex<double> lambda, u0;

    //! The complex unit \\( i = \\sqrt{-1} \\)
    const complex<double> i_complex = complex<double>(0, 1);

     //! Error at the final time. For the scalar example, an analytical solution is known.
    double error;

     //! Counters for how often `f_expl_eval`, `f_impl_eval` and `impl_solve` are called.
    size_t n_f_expl_eval, n_f_impl_eval, n_impl_solve;

  public:
    /**
     * Generic constructor; initialize all function call counters with zero.
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
     * Upon destruction, report final error and number of function calls
     */
    virtual ~ScalarSweeper()
    {
      cout << "Final error:                    " << scientific << this->error << endl;
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
      auto& qend = pfasst::encap::as_vector<complex<double>, time>(this->get_end_state());

      complex_vector_type qex(qend.size());

      this->exact(qex, t);
      double max_err = abs(qend[0] - qex[0]) / abs(qex[0]);
      cout << "err: " << scientific << max_err << endl;
      this->error = max_err;
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
     * Computes the exact solution \\( u_0 \\exp \\left( \\lambda*t \\right) \\) at a given time t.
     *
     * @param[in] t Time
     */
    void exact(complex_vector_type& q, time t)
    {
      q[0] = this->u0 * exp(this->lambda * t);
    }

    void exact(shared_ptr<encap_type> q_encap, time t)
    {
      auto& q = pfasst::encap::as_vector<complex<double>, time>(q_encap);
      this->exact(q, t);
    }

    /**
     * Evaluate the explicit part of the right hand side: Multiply with \\( \\text{imag}(\\lambda) \\)
     */
    void f_expl_eval(shared_ptr<encap_type> f_encap,
                     shared_ptr<encap_type> q_encap, time t) override
    {
      UNUSED(t);
      auto& f = pfasst::encap::as_vector<complex<double>, time>(f_encap);
      auto& q = pfasst::encap::as_vector<complex<double>, time>(q_encap);

      // f_expl = multiply with imaginary part of lambda
      f[0] = this->i_complex * imag(this->lambda) * q[0];

      this->n_f_expl_eval++;
    }

    /**
     * Evaluate the implicit part of the right hand side: Multiply with \\( \\text{real}(\\lambda) \\)
     */
    void f_impl_eval(shared_ptr<encap_type> f_encap,
                     shared_ptr<encap_type> q_encap, time t) override
    {
      UNUSED(t);
      auto& f = pfasst::encap::as_vector<complex<double>, time>(f_encap);
      auto& q = pfasst::encap::as_vector<complex<double>, time>(q_encap);

      // f_impl = multiply with real part of lambda
      f[0] = real(this->lambda) * q[0];

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
      auto& f = pfasst::encap::as_vector<complex<double>, time>(f_encap);
      auto& q = pfasst::encap::as_vector<complex<double>, time>(q_encap);
      auto& rhs = pfasst::encap::as_vector<complex<double>, time>(rhs_encap);

      // invert f_impl = multiply with inverse of real part of lambda
      double inv = 1.0 / (1.0 - double(dt) * real(this->lambda));
      q[0] = inv * rhs[0];

      // compute f_impl_eval of q[0], i.e. multiply with real(lambda)
      f[0] = real(this->lambda) * q[0];

      this->n_impl_solve++;
    }
};

#endif  // _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_
