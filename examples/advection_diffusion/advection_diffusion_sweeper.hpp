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

#ifndef PI
#define PI 3.1415926535897932385
#endif

using namespace std;

template<typename time = pfasst::time_precision>
class AdvectionDiffusionSweeper
  : public pfasst::encap::IMEXSweeper<time>
{
    typedef pfasst::encap::Encapsulation<time> Encapsulation;

    FFT fft;

    vector<complex<double>> ddx, lap;

    double v  = 1.0;
    time   t0 = 1.0;
    double nu = 0.02;
    size_t nf1evals = 0;

    inline DVectorT& as_dvector(shared_ptr<Encapsulation> x)
    {
      shared_ptr<DVectorT> y = dynamic_pointer_cast<DVectorT>(x); assert(y);
      return *y.get();
    }

  public:
    AdvectionDiffusionSweeper(size_t nvars)
    {
      ddx.resize(nvars);
      lap.resize(nvars);
      for (size_t i = 0; i < nvars; i++) {
        double kx = 2 * PI * ((i <= nvars / 2) ? int(i) : int(i) - int(nvars));
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
      this->exact(as_dvector(q), t);
    }

    void exact(DVectorT& q, time t)
    {
      size_t n = q.size();
      double a = 1.0 / sqrt(4 * PI * nu * (t + t0));

      for (size_t i = 0; i < n; i++) {
        q[i] = 0.0;
      }

      for (int ii = -2; ii < 3; ii++) {
        for (size_t i = 0; i < n; i++) {
          double x = double(i) / n - 0.5 + ii - t * v;
          q[i] += a * exp(-x * x / (4 * nu * (t + t0)));
        }
      }
    }

    void echo_error(time t, bool predict = false)
    {
      DVectorT& qend = as_dvector(this->get_end_state());
      DVectorT  qex(qend.size());

      exact(qex, t);

      double max = 0.0;
      for (size_t i = 0; i < qend.size(); i++) {
        double d = abs(qend[i] - qex[i]);
        if (d > max) { max = d; }
      }
      cout << "err: " << scientific << max
           << " (" << qend.size() << ", " << predict << ")"
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

    void f1eval(shared_ptr<Encapsulation> _f, shared_ptr<Encapsulation> _q, time /*t*/)
    {
      DVectorT& q = as_dvector(_q);
      DVectorT& f = as_dvector(_f);

      double c = -v / double(q.size());

      auto* z = fft.forward(q);
      for (size_t i = 0; i < q.size(); i++) {
        z[i] *= c * ddx[i];
      }
      fft.backward(f);

      nf1evals++;
    }

    void f2eval(shared_ptr<Encapsulation> _f, shared_ptr<Encapsulation> _q, time /*t*/)
    {
      DVectorT& q = as_dvector(_q);
      DVectorT& f = as_dvector(_f);

      double c = nu / double(q.size());

      auto* z = fft.forward(q);
      for (size_t i = 0; i < q.size(); i++) {
        z[i] *= c * lap[i];
      }
      fft.backward(f);
    }

    void f2comp(shared_ptr<Encapsulation> _f, shared_ptr<Encapsulation> _q, time /*t*/, time dt,
                shared_ptr<Encapsulation> _rhs)
    {
      DVectorT& q = as_dvector(_q);
      DVectorT& f = as_dvector(_f);
      DVectorT& rhs = as_dvector(_rhs);

      auto* z = fft.forward(rhs);
      for (size_t i = 0; i < q.size(); i++) {
        z[i] /= (1.0 - nu * double(dt) * lap[i]) * double(q.size());
      }
      fft.backward(q);

      for (size_t i = 0; i < q.size(); i++) {
        f[i] = (q[i] - rhs[i]) / double(dt);
      }

    }

};

#endif
