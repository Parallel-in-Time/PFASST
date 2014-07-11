/*
 * Advection/diffusion sweeper.
 */

#ifndef _ADVECTION_DIFFUSION_SWEEPER_HPP_
#define _ADVECTION_DIFFUSION_SWEEPER_HPP_

#include <complex>
#include <vector>
#include <cassert>
#include <ostream>

#include <pfasst/encap/imex_sweeper.hpp>

#include "fft.hpp"

#define PI     3.1415926535897932385
#define TWO_PI 6.2831853071795864769

using namespace std;

template<typename time = pfasst::time_precision>
class AdvectionDiffusionSweeper
  : public pfasst::encap::IMEXSweeper<time>
{
    typedef pfasst::encap::Encapsulation<time> Encapsulation;
    typedef pfasst::encap::VectorEncapsulation<double> DVectorT;

    FFT fft;

    vector<complex<double>> ddx, lap;

    double v  = 1.0;
    time   t0 = 1.0;
    double nu = 0.02;
    size_t nf1evals = 0;

  public:
    AdvectionDiffusionSweeper(size_t nvars)
    {
      ddx.resize(nvars);
      lap.resize(nvars);
      for (size_t i = 0; i < nvars; i++) {
        double kx = TWO_PI * ((i <= nvars / 2) ? int(i) : int(i) - int(nvars));
        ddx[i] = complex<double>(0.0, 1.0) * kx;
        lap[i] = (kx * kx < 1e-13) ? 0.0 : -kx * kx;
      }
    }

    ~AdvectionDiffusionSweeper()
    {
      cout << "number of f1 evals: " << nf1evals << endl;
    }

    void exact(shared_ptr<Encapsulation> q, time t)
    {
      shared_ptr<DVectorT> q_cast = dynamic_pointer_cast<DVectorT>(q);
      assert(q_cast);
      this->exact(q_cast, t);
    }

    void exact(shared_ptr<DVectorT> q, time t)
    {
      size_t n = q->size();
      double a = 1.0 / sqrt(4 * PI * nu * (t + t0));

      for (size_t i = 0; i < n; i++) {
        q->data()[i] = 0.0;
      }

      for (int ii = -2; ii < 3; ii++) {
        for (size_t i = 0; i < n; i++) {
          double x = double(i) / n - 0.5 + ii - t * v;
          q->data()[i] += a * exp(-x * x / (4 * nu * (t + t0)));
        }
      }
    }

    void echo_error(time t, bool predict = false)
    {
      shared_ptr<DVectorT> qend = dynamic_pointer_cast<DVectorT>(this->get_state(this->get_nodes().size() - 1));
      assert(qend);
      shared_ptr<DVectorT> qex = make_shared<DVectorT>(qend->size());

      exact(qex, t);

      double max = 0.0;
      for (size_t i = 0; i < qend->size(); i++) {
        double d = abs(qend->data()[i] - qex->data()[i]);
        if (d > max) { max = d; }
      }
      cout << "err: " << scientific << max
           << " (" << qend->size() << ", " << predict << ")"
           << endl;
    }

    void predict(time t, time dt, bool initial)
    {
      pfasst::encap::IMEXSweeper<time>::predict(t, dt, initial);
      echo_error(t + dt, true);
    }

    void sweep(time t, time dt)
    {
      pfasst::encap::IMEXSweeper<time>::sweep(t, dt);
      echo_error(t + dt);
    }

    void f1eval(shared_ptr<Encapsulation> f, shared_ptr<Encapsulation> q, time t)
    {
      shared_ptr<DVectorT> f_cast = dynamic_pointer_cast<DVectorT>(f);
      assert(f_cast);
      shared_ptr<DVectorT> q_cast = dynamic_pointer_cast<DVectorT>(q);
      assert(q_cast);

      this->f1eval(f_cast, q_cast, t);
    }

    void f1eval(shared_ptr<DVectorT> f, shared_ptr<DVectorT> q, time t)
    {
      double c = -v / double(q->size());

      auto* z = fft.forward(q);
      for (size_t i = 0; i < q->size(); i++) {
        z[i] *= c * ddx[i];
      }
      fft.backward(f);

      nf1evals++;
    }

    void f2eval(shared_ptr<Encapsulation> f, shared_ptr<Encapsulation> q, time t)
    {
      shared_ptr<DVectorT> f_cast = dynamic_pointer_cast<DVectorT>(f);
      assert(f_cast);
      shared_ptr<DVectorT> q_cast = dynamic_pointer_cast<DVectorT>(q);
      assert(q_cast);

      this->f2eval(f_cast, q_cast, t);
    }

    void f2eval(shared_ptr<DVectorT> f, shared_ptr<DVectorT> q, time t)
    {
      double c = nu / double(q->size());

      auto* z = fft.forward(q);
      for (size_t i = 0; i < q->size(); i++) {
        z[i] *= c * lap[i];
      }
      fft.backward(f);
    }

    void f2comp(shared_ptr<Encapsulation> f, shared_ptr<Encapsulation> q, time t, time dt,
                shared_ptr<Encapsulation> rhs)
    {
      shared_ptr<DVectorT> f_cast   = dynamic_pointer_cast<DVectorT>(f);
      assert(f_cast);
      shared_ptr<DVectorT> q_cast   = dynamic_pointer_cast<DVectorT>(q);
      assert(q_cast);
      shared_ptr<DVectorT> rhs_cast = dynamic_pointer_cast<DVectorT>(rhs);
      assert(rhs_cast);

      this->f2comp(f_cast, q_cast, t, dt, rhs_cast);
    }

    void f2comp(shared_ptr<DVectorT> f, shared_ptr<DVectorT> q, time t, time dt,
                shared_ptr<DVectorT> rhs)
    {
      auto* z = fft.forward(rhs);
      for (size_t i = 0; i < q->size(); i++) {
        z[i] /= (1.0 - nu * double(dt) * lap[i]) * double(q->size());
      }
      fft.backward(q);

      for (size_t i = 0; i < q->size(); i++) {
        f->data()[i] = (q->data()[i] - rhs->data()[i]) / double(dt);
      }
    }
};

#endif
