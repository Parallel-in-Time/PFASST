/*
 * Sweeper for scalar test equation
 */

#ifndef _SCALAR_SWEEPER_HPP_
#define _SCALAR_SWEEPER_HPP_

#include <complex>
#include <vector>
#include <pfasst/encap/imex_sweeper.hpp>
#include <pfasst/encap/vector.hpp>

using namespace std;

template<typename time = pfasst::time_precision>
class ScalarSweeper : public pfasst::encap::IMEXSweeper<time>
{

    typedef pfasst::encap::Encapsulation<time> Encapsulation;
    typedef pfasst::encap::VectorEncapsulation<complex<double>> CVectorT;

  public:

    ScalarSweeper(complex<double> lambda, complex<double> y0)
    {
      this->_lambda  = lambda;
      this->_y0      = y0;
      this->_nf1eval = 0;
      this->_nf2eval = 0;
      this->_nf2comp = 0;
    }

    virtual ~ScalarSweeper()
    {
      cout << "Number of calls to f1eval: " << this->_nf1eval << endl;
      cout << "Number of calls to f2eval: " << this->_nf2eval << endl;
      cout << "Number of calls to f2comp: " << this->_nf2comp << endl;
    }

    void echo_error(time t)
    {
      auto& qend = pfasst::encap::as_vector<complex<double>,time>(this->get_end_state());

      CVectorT qex(qend.size());

      this->exact(qex, t);
      double max_err = abs(qend[0] - qex[0]) / abs(qex[0]);
      cout << "err: " << scientific << max_err << endl;
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

    void exact(CVectorT& q, time t)
    {
      q[0] = _y0 * exp(_lambda * t);
    }

    void exact(shared_ptr<Encapsulation> q_encap, time t)
    {
      auto& q = pfasst::encap::as_vector<complex<double>,time>(q_encap);
      exact(q, t);
    }

    void f1eval(shared_ptr<Encapsulation> f_encap, shared_ptr<Encapsulation> q_encap, time) override
    {
      auto& f = pfasst::encap::as_vector<complex<double>,time>(f_encap);
      auto& q = pfasst::encap::as_vector<complex<double>,time>(q_encap);

      // f1 = multiply with imaginary part of lambda
      f[0] = i_complex * imag(this->_lambda) * q[0];
      this->_nf1eval++;
    }

    void f2eval(shared_ptr<Encapsulation> f_encap, shared_ptr<Encapsulation> q_encap, time) override
    {
      auto& f = pfasst::encap::as_vector<complex<double>,time>(f_encap);
      auto& q = pfasst::encap::as_vector<complex<double>,time>(q_encap);

      // f2 = multiply with real part of lambda
      f[0] = real(this->_lambda) * q[0];
      this->_nf2eval++;
    }

    void f2comp(shared_ptr<Encapsulation> f_encap, shared_ptr<Encapsulation> q_encap, time, time dt,
                shared_ptr<Encapsulation> rhs_encap) override
    {
      auto& f = pfasst::encap::as_vector<complex<double>,time>(f_encap);
      auto& q = pfasst::encap::as_vector<complex<double>,time>(q_encap);
      auto& rhs = pfasst::encap::as_vector<complex<double>,time>(rhs_encap);

      // invert f2=multiply with inverse of real part of lambda
      double inv = 1 / (1 - double(dt) * real(this->_lambda));
      q[0] = inv * rhs[0];
      f[0] = real(this->_lambda) * q[0];
      this->_nf2comp++;
    }


  private:

    complex<double> _lambda, _y0;
    int _nf1eval, _nf2eval, _nf2comp;
    const complex<double> i_complex = complex<double>(0, 1);

};
#endif
