/*
 * Advection/diffusion sweeper.
 */

#ifndef _ADVECTION_DIFFUSION_SWEEPER_HPP_
#define _ADVECTION_DIFFUSION_SWEEPER_HPP_

#include <complex>
#include <vector>

#include <pfasst/encap/imex_sweeper.hpp>
#include "fft.hpp"

#define PI     3.1415926535897932385
#define TWO_PI 6.2831853071795864769


using namespace std;
using pfasst::encap::Encapsulation;


template<typename ScalarT, typename timeT>
class AdvectionDiffusionSweeper : public pfasst::encap::IMEXSweeper<ScalarT, timeT>
{
    typedef pfasst::encap::VectorEncapsulation<ScalarT, timeT> DVectorT;
    FFT<ScalarT, timeT> fft;

    vector<complex<ScalarT>> ddx, lap;

    ScalarT v  = 1.0;
    timeT   t0 = 1.0;
    ScalarT nu = 0.02;
    int    nf1evals = 0;

  public:

    AdvectionDiffusionSweeper(int nvars)
    {
      ddx.resize(nvars);
      lap.resize(nvars);

      for (int i = 0; i < nvars; i++) {
        ScalarT kx = TWO_PI * (i <= nvars / 2 ? i : i - nvars);
        ddx[i] = complex<ScalarT>(0.0, 1.0) * kx;
        lap[i] = (kx * kx < 1e-13) ? 0.0 : -kx * kx;
      }
    }

    ~AdvectionDiffusionSweeper()
    {
      cout << "number of f1 evals: " << nf1evals << endl;
    }

    void exact(Encapsulation<ScalarT, timeT>* q, ScalarT t)
    {
      exact(*dynamic_cast<DVectorT*>(q), t);
    }

    void exact(DVectorT& q, ScalarT t)
    {
      int    n = q.size();
      ScalarT a = 1.0 / sqrt(4 * PI * nu * (t + t0));

      for (int i = 0; i < n; i++)
      { q[i] = 0.0; }

      for (int ii = -2; ii < 3; ii++) {
        for (int i = 0; i < n; i++) {
          ScalarT x = ScalarT(i) / n - 0.5 + ii - t * v;
          q[i] += a * exp(-x * x / (4 * nu * (t + t0)));
        }
      }
    }

    void echo_error(timeT t, bool predict = false)
    {
      auto& qend = *dynamic_cast<DVectorT*>(this->get_state(this->get_nodes().size() - 1));
      auto  qex  = DVectorT(qend.size());

      exact(qex, t);

      ScalarT max = 0.0;

      for (int i = 0; i < qend.size(); i++) {
        ScalarT d = abs(qend[i] - qex[i]);

        if (d > max)
        { max = d; }
      }

      cout << "err: " << scientific << max << " (" << qend.size() << ", " << predict << ")" << endl;
    }

    void predict(timeT t, timeT dt, bool initial)
    {
      pfasst::encap::IMEXSweeper<ScalarT, timeT>::predict(t, dt, initial);
      echo_error(t + dt, true);
    }

    void sweep(timeT t, timeT dt)
    {
      pfasst::encap::IMEXSweeper<ScalarT, timeT>::sweep(t, dt);
      echo_error(t + dt);
    }

    void f1eval(Encapsulation<ScalarT, timeT>* F, Encapsulation<ScalarT, timeT>* Q, timeT t)
    {
      auto& f = *dynamic_cast<DVectorT*>(F);
      auto& q = *dynamic_cast<DVectorT*>(Q);

      ScalarT c = -v / ScalarT(q.size());

      auto* z = fft.forward(q);

      for (int i = 0; i < q.size(); i++)
      { z[i] *= c * ddx[i]; }

      fft.backward(f);

      nf1evals++;
    }

    void f2eval(Encapsulation<ScalarT, timeT>* F, Encapsulation<ScalarT, timeT>* Q, timeT t)
    {
      auto& f = *dynamic_cast<DVectorT*>(F);
      auto& q = *dynamic_cast<DVectorT*>(Q);

      ScalarT c = nu / ScalarT(q.size());

      auto* z = fft.forward(q);

      for (int i = 0; i < q.size(); i++)
      { z[i] *= c * lap[i]; }

      fft.backward(f);
    }

    void f2comp(Encapsulation<ScalarT, timeT>* F, Encapsulation<ScalarT, timeT>* Q, timeT t, timeT dt,
                Encapsulation<ScalarT, timeT>* RHS)
    {
      auto& f   = *dynamic_cast<DVectorT*>(F);
      auto& q   = *dynamic_cast<DVectorT*>(Q);
      auto& rhs = *dynamic_cast<DVectorT*>(RHS);

      auto* z = fft.forward(rhs);

      for (int i = 0; i < q.size(); i++)
      { z[i] /= (1.0 - nu * dt * lap[i]) * ScalarT(q.size()); }

      fft.backward(q);

      for (int i = 0; i < q.size(); i++)
      { f[i] = (q[i] - rhs[i]) / dt; }
    }

};

#endif
