/*
 * Advection/diffusion sweeper.
 */

#ifndef _ADVECTION_DIFFUSION_SWEEPER_HPP_
#define _ADVECTION_DIFFUSION_SWEEPER_HPP_

#include <complex>
#include <vector>

#include <pfasst/encap/imex_sweeper.hpp>
#include "fft.hpp"

#define pi     3.1415926535897932385
#define two_pi 6.2831853071795864769


using namespace std;
using pfasst::encap::Encapsulation;


template<typename scalar, typename time>
class AdvectionDiffusionSweeper : public pfasst::encap::IMEXSweeper<scalar,time> {

  using dvector = pfasst::encap::VectorEncapsulation<scalar,time>;
  FFT<scalar,time> fft;

  vector<complex<scalar>> ddx, lap;

  scalar v  = 1.0;
  time   t0 = 1.0;
  scalar nu = 0.02;
  int    nf1evals = 0;

public:

  AdvectionDiffusionSweeper(int nvars)
  {
    ddx.resize(nvars);
    lap.resize(nvars);
    for (int i=0; i<nvars; i++) {
      scalar kx = two_pi * ( i <= nvars/2 ? i : i-nvars );
      ddx[i] = complex<scalar>(0.0, 1.0) * kx;
      lap[i] = (kx*kx < 1e-13) ? 0.0 : -kx*kx;
    }
  }

  ~AdvectionDiffusionSweeper()
  {
    cout << "number of f1 evals: " << nf1evals << endl;
  }

  void exact(Encapsulation<scalar,time>* q, scalar t)
  {
    exact(*dynamic_cast<dvector*>(q), t);
  }

  void exact(dvector& q, scalar t)
  {
    int    n = q.size();
    scalar a = 1.0/sqrt(4*pi*nu*(t+t0));

    for (int i=0; i<n; i++)
      q[i] = 0.0;

    for (int ii=-2; ii<3; ii++) {
      for (int i=0; i<n; i++) {
	scalar x = scalar(i)/n - 0.5 + ii - t * v;
	q[i] += a * exp(-x*x/(4*nu*(t+t0)));
      }
    }
  }

  void echo_error(time t, bool predict=false)
  {
    auto& qend = *dynamic_cast<dvector*>(this->get_state(this->get_nodes().size()-1));
    auto  qex  = dvector(qend.size());

    exact(qex, t);

    scalar max = 0.0;
    for (int i=0; i<qend.size(); i++) {
      scalar d = abs(qend[i]-qex[i]);
      if (d > max)
	max = d;
    }
    cout << "err: " << scientific << max << " (" << qend.size() << ", " << predict << ")" << endl;
  }

  void predict(time t, time dt, bool initial)
  {
    pfasst::encap::IMEXSweeper<scalar,time>::predict(t, dt, initial);
    echo_error(t+dt, true);
  }

  void sweep(time t, time dt)
  {
    pfasst::encap::IMEXSweeper<scalar,time>::sweep(t, dt);
    echo_error(t+dt);
  }

  void f1eval(Encapsulation<scalar,time> *F, Encapsulation<scalar,time> *Q, time t)
  {
    auto& f = *dynamic_cast<dvector*>(F);
    auto& q = *dynamic_cast<dvector*>(Q);

    cout << "f1eval " << scientific << q.norm0() << endl;

    scalar c = -v / scalar(q.size());

    auto* z = fft.forward(q);
    for (int i=0; i<q.size(); i++)
      z[i] *= c * ddx[i];
    fft.backward(f);

    nf1evals++;
  }

  void f2eval(Encapsulation<scalar,time> *F, Encapsulation<scalar,time> *Q, time t)
  {
    auto& f = *dynamic_cast<dvector*>(F);
    auto& q = *dynamic_cast<dvector*>(Q);

    scalar c = nu / scalar(q.size());

    auto* z = fft.forward(q);
    for (int i=0; i<q.size(); i++)
      z[i] *= c * lap[i];
    fft.backward(f);
  }

  void f2comp(Encapsulation<scalar,time> *F, Encapsulation<scalar,time> *Q, time t, time dt, Encapsulation<scalar,time> *RHS)
  {
    auto& f   = *dynamic_cast<dvector*>(F);
    auto& q   = *dynamic_cast<dvector*>(Q);
    auto& rhs = *dynamic_cast<dvector*>(RHS);

    auto* z = fft.forward(rhs);
    for (int i=0; i<q.size(); i++)
      z[i] /= (1.0 - nu * dt * lap[i]) * scalar(q.size());
    fft.backward(q);

    for (int i=0; i<q.size(); i++)
      f[i] = (q[i] - rhs[i]) / dt;
  }

};

#endif
