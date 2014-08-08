/*
 * Sweeper for scalar test equation
 */

#ifndef _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_
#define _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_

#include <complex>
#include <vector>

#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/vector.hpp>

using namespace std;

template<typename time = pfasst::time_precision>
class ScalarSweeper
  : public pfasst::encap::IMEXSweeper<time>
{
  private:
    typedef pfasst::encap::Encapsulation<time> encap_type;
    typedef pfasst::encap::VectorEncapsulation<complex<double>> complex_vector_type;

    complex<double> lambda, u0;
    size_t n_f_expl_eval, n_f_impl_eval, n_impl_solve;
    const complex<double> i_complex = complex<double>(0, 1);
    double error;

  public:
    ScalarSweeper(complex<double> lambda, complex<double> u0)
      :   lambda(lambda)
        , u0(u0)
        , n_f_expl_eval(0)
        , n_f_impl_eval(0)
        , n_impl_solve(0)
        , error(0.0)
    {}

    virtual ~ScalarSweeper()
    {
      cout << "Final error:                    " << scientific << this->error << endl;
      cout << "Number of explicit evaluations: " << this->n_f_expl_eval << endl;
      cout << "Number of implicit evaluations: " << this->n_f_impl_eval << endl;
      cout << "Number of implicit solves:      " << this->n_impl_solve << endl;
    }

    void echo_error(time t)
    {
      auto& qend = pfasst::encap::as_vector<complex<double>, time>(this->get_end_state());

      complex_vector_type qex(qend.size());

      this->exact(qex, t);
      double max_err = abs(qend[0] - qex[0]) / abs(qex[0]);
      cout << "err: " << scientific << max_err << endl;
      this->error = max_err;
    }

    double get_errors()
    {
      return this->error;
    }

    void predict(bool initial) override
    {
      pfasst::encap::IMEXSweeper<time>::predict(initial);
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_time_step();
      this->echo_error(t + dt);
    }

    void sweep() override
    {
      pfasst::encap::IMEXSweeper<time>::sweep();
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_time_step();
      this->echo_error(t + dt);
    }

    void exact(complex_vector_type& q, time t)
    {
      q[0] = this->u0 * exp(this->lambda * t);
    }

    void exact(shared_ptr<encap_type> q_encap, time t)
    {
      auto& q = pfasst::encap::as_vector<complex<double>, time>(q_encap);
      this->exact(q, t);
    }

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
      f[0] = real(this->lambda) * q[0];

      this->n_impl_solve++;
    }
};

#endif  // _EXAMPLES__SCALAR__SCALAR_SWEEPER_HPP_
