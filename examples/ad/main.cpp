/*
 * Advection/diffusion example using the encapsulated IMEX sweeper.
 */

#include <algorithm>
#include <cmath>
#include <complex>

#include <pfasst.hpp>
#include <pfasst-imex.hpp>

#include <fftw3.h>

#define pi 3.1415926535897932385
#define two_pi 6.2831853071795864769

using namespace std;

//
// config
//

const int nlevs = 1;
const int ndofs = 512;
const int nnodes = 5;
const int xrat = 2;
const int trat = 2;

const int    nsteps = 32;
const double dt = 0.01;

using namespace std;
using pfasst::Encapsulation;
using dvector = pfasst::VectorEncapsulation<double>;


//
// advection/diffusion sweeper
//

template<typename timeT>
class ADIMEX : public pfasst::IMEX<timeT> {

  fftw_plan     ffft;
  fftw_plan     ifft;
  fftw_complex* wk;

  vector<complex<double>> ddx, lap;

  double v  = 1.0;
  double t0 = 1.0;
  double nu = 0.02;

public:

  ADIMEX(vector<timeT> nodes, pfasst::VectorFactory<timeT> *factory)
  {
    this->set_nodes(nodes);
    this->set_factory(factory);

    int nvars = factory->dofs();

    wk   = fftw_alloc_complex(nvars);
    ffft = fftw_plan_dft_1d(nvars, wk, wk, FFTW_FORWARD, FFTW_ESTIMATE);
    ifft = fftw_plan_dft_1d(nvars, wk, wk, FFTW_BACKWARD, FFTW_ESTIMATE);

    ddx.resize(nvars);
    lap.resize(nvars);
    for (int i=0; i<nvars; i++) {
      double kx = two_pi * ( i <= nvars/2 ? i : i-nvars );
      ddx[i] = complex<double>(0.0, 1.0) * kx;
      lap[i] = (kx*kx < 1e-13) ? 0.0 : -kx*kx;
    }

  }

  ~ADIMEX() {
    fftw_destroy_plan(ffft);
    fftw_destroy_plan(ifft);
    fftw_free(wk);
  }

  void exact(dvector& q, double t)
  {
    int    n = q.size();
    double a = 1.0/sqrt(4*pi*nu*(t+t0));

    for (int i=0; i<n; i++)
      q[i] = 0.0;

    for (int ii=-2; ii<3; ii++) {
      for (int i=0; i<n; i++) {
	double x = double(i)/n - 0.5 + ii - t * v;
	q[i] += a * exp(-x*x/(4*nu*(t+t0)));
      }
    }
  }

  void sweep(timeT t, timeT dt)
  {
    pfasst::IMEX<timeT>::sweep(t, dt);

    auto& qend = *dynamic_cast<dvector*>(this->get_qend());
    auto  qex  = dvector(qend.size());

    exact(qex, t+dt);

    double max = 0.0;
    for (int i=0; i<qend.size(); i++) {
      double d = abs(qend[i]-qex[i]);
      if (d > max)
	max = d;
    }
    cout << "err: " << max << endl;
  }

  void f1eval(Encapsulation *F, Encapsulation *Q, timeT t)
  {
    auto& f = *dynamic_cast<dvector*>(F);
    auto& q = *dynamic_cast<dvector*>(Q);

    complex<double>* z = reinterpret_cast<complex<double>*>(wk);
    double c = -v / double(q.size());

    copy(q.begin(), q.end(), z);
    fftw_execute_dft(ffft, wk, wk);
    for (int i=0; i<q.size(); i++)
      z[i] *= c * ddx[i];
    fftw_execute_dft(ifft, wk, wk);

    for (int i=0; i<q.size(); i++)
      f[i] = real(z[i]);
  }

  void f2eval(Encapsulation *F, Encapsulation *Q, timeT t)
  {
    auto& f = *dynamic_cast<dvector*>(F);
    auto& q = *dynamic_cast<dvector*>(Q);

    complex<double>* z = reinterpret_cast<complex<double>*>(wk);
    double c = nu / double(q.size());

    copy(q.begin(), q.end(), z);
    fftw_execute_dft(ffft, wk, wk);
    for (int i=0; i<q.size(); i++)
      z[i] *= c * lap[i];
    fftw_execute_dft(ifft, wk, wk);

    for (int i=0; i<q.size(); i++)
      f[i] = real(z[i]);
  }

  void f2comp(Encapsulation *F, Encapsulation *Q, timeT t, timeT dt, Encapsulation *RHS)
  {
    auto& f   = *dynamic_cast<dvector*>(F);
    auto& q   = *dynamic_cast<dvector*>(Q);
    auto& rhs = *dynamic_cast<dvector*>(RHS);

    complex<double>* z = reinterpret_cast<complex<double>*>(wk);

    copy(rhs.begin(), rhs.end(), z);
    fftw_execute_dft(ffft, wk, wk);
    for (int i=0; i<q.size(); i++)
      z[i] /= (1.0 - nu * dt * lap[i]) * double(q.size());
    fftw_execute_dft(ifft, wk, wk);

    for (int i=0; i<q.size(); i++) {
      q[i] = real(z[i]);
      f[i] = (q[i] - rhs[i]) / dt;
    }

  }

};


//
// main
//

int main(int argc, char **argv)
{
  int ndofs  = 512;
  int nnodes = 5;

  if (nlevs == 1) {
    pfasst::SDC<double> sdc;

    auto  nodes   = pfasst::compute_nodes<double>(nnodes, "gauss-lobatto");
    auto* factory = new pfasst::VectorFactory<double>(ndofs);
    auto* sweeper = new ADIMEX<double>(nodes, factory);

    sdc.add_level(sweeper);
    sdc.set_duration(dt, nsteps, 4);
    sdc.setup();

    dvector q0(ndofs);
    sweeper->exact(q0, 0.0);
    sweeper->set_q0(&q0);

    sdc.run();
  } else {
    // ...
  }

  return 0;
}
