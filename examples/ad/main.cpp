/*
 * Advection/diffusion example using the encapsulated IMEX sweeper.
 */

#define PFASST_ENABLE_GNUPLOT

#include <algorithm>
#include <complex>
#include <cmath>
#include <map>
#include <tuple>

#include <pfasst.hpp>
#include <pfasst/encap/imex.hpp>
#include <pfasst/encap/vector.hpp>

#include <fftw3.h>

#define pi 3.1415926535897932385
#define two_pi 6.2831853071795864769

using namespace std;

//
// config
//

const int nlevs = 2;
const int xrat = 2;
const int trat = 2;

const int    nsteps = 1;
const double dt = 0.01;

typedef double scalar;

using namespace std;
using pfasst::encap::Encapsulation;
using dvector = pfasst::encap::VectorEncapsulation<scalar,double>;

//
// fft helper
//
class FFT {

  struct workspace {
    fftw_plan        ffft;
    fftw_plan        ifft;
    fftw_complex*    wk;
    complex<scalar>* z;
  };

  map<int,workspace*> workspaces;

public:

  ~FFT()
  {
    // XXX
  }

  workspace* get_workspace(int ndofs)
  {
    if (workspaces.find(ndofs) == workspaces.end()) {
      workspace* wk = new workspace;
      wk->wk = fftw_alloc_complex(ndofs);
      wk->ffft = fftw_plan_dft_1d(ndofs, wk->wk, wk->wk, FFTW_FORWARD, FFTW_ESTIMATE);
      wk->ifft = fftw_plan_dft_1d(ndofs, wk->wk, wk->wk, FFTW_BACKWARD, FFTW_ESTIMATE);
      wk->z = reinterpret_cast<complex<scalar>*>(wk->wk);
      workspaces.insert(pair<int,workspace*>(ndofs, wk));
    }

    return workspaces[ndofs];
  }

  complex<double>* forward(const dvector& x)
  {
    workspace *wk = get_workspace(x.size());
    for (unsigned int i=0; i<x.size(); i++)
      wk->z[i] = x[i];
    fftw_execute_dft(wk->ffft, wk->wk, wk->wk);
    return wk->z;
  }

  void backward(dvector &x)
  {
    workspace *wk = get_workspace(x.size());
    fftw_execute_dft(wk->ifft, wk->wk, wk->wk);
    for (unsigned int i=0; i<x.size(); i++)
      x[i] = real(wk->z[i]);
  }

} fft;


//
// advection/diffusion sweeper
//

template<typename time>
class ADIMEX : public pfasst::imex::IMEX<time> {

  vector<complex<scalar>> ddx, lap;

  scalar v  = 1.0;
  time   t0 = 1.0;
  scalar nu = 0.02;

public:

  ADIMEX(int nvars)
  {
    ddx.resize(nvars);
    lap.resize(nvars);
    for (int i=0; i<nvars; i++) {
      scalar kx = two_pi * ( i <= nvars/2 ? i : i-nvars );
      ddx[i] = complex<scalar>(0.0, 1.0) * kx;
      lap[i] = (kx*kx < 1e-13) ? 0.0 : -kx*kx;
    }

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

  void echo_error(time t)
  {
    auto& qend = *dynamic_cast<dvector*>(this->get_q(this->get_nodes().size()-1));
    auto  qex  = dvector(qend.size());

    exact(qex, t);

    scalar max = 0.0;
    for (int i=0; i<qend.size(); i++) {
      scalar d = abs(qend[i]-qex[i]);
      if (d > max)
	max = d;
    }
    cout << "err: " << scientific << max << " (" << qend.size() << ")" << endl;
  }

  void predict(time t, time dt)
  {
    // cout << "PREDICTOR" << endl;
    pfasst::imex::IMEX<time>::predict(t, dt);
    echo_error(t+dt);
  }

  void sweep(time t, time dt)
  {
    // cout << "SWEEPING" << endl;
    pfasst::imex::IMEX<time>::sweep(t, dt);
    echo_error(t+dt);
  }

  void f1eval(Encapsulation<scalar> *F, Encapsulation<scalar> *Q, time t)
  {
    auto& f = *dynamic_cast<dvector*>(F);
    auto& q = *dynamic_cast<dvector*>(Q);

    scalar c = -v / scalar(q.size());

    auto* z = fft.forward(q);
    for (int i=0; i<q.size(); i++)
      z[i] *= c * ddx[i];
    fft.backward(f);
  }

  void f2eval(Encapsulation<scalar> *F, Encapsulation<scalar> *Q, time t)
  {
    auto& f = *dynamic_cast<dvector*>(F);
    auto& q = *dynamic_cast<dvector*>(Q);

    scalar c = nu / scalar(q.size());

    auto* z = fft.forward(q);
    for (int i=0; i<q.size(); i++)
      z[i] *= c * lap[i];
    fft.backward(f);
  }

  void f2comp(Encapsulation<scalar> *F, Encapsulation<scalar> *Q, time t, time dt, Encapsulation<scalar> *RHS)
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

template<typename scalar, typename time>
class ADTRANS : public pfasst::encap::PolyInterpMixin<time> {
public:

  void interpolate(Encapsulation<scalar> *dst, const Encapsulation<scalar> *src) {
    auto& crse = *dynamic_cast<const dvector*>(src);
    auto& fine = *dynamic_cast<dvector*>(dst);

    auto* crse_z = fft.forward(crse);
    auto* fine_z = fft.get_workspace(fine.size())->z;

    for (int i=0; i<fine.size(); i++)
      fine_z[i] = 0.0;

    double c = 1.0 / crse.size();

    for (int i=0; i<crse.size()/2; i++)
      fine_z[i] = c * crse_z[i];

    for (int i=1; i<crse.size()/2; i++)
      fine_z[fine.size()-crse.size()/2+i] = c * crse_z[crse.size()/2+i];

    fft.backward(fine);
  }

  void restrict(Encapsulation<scalar> *dst, const Encapsulation<scalar> *src) {
    auto& crse = *dynamic_cast<dvector*>(dst);
    auto& fine = *dynamic_cast<const dvector*>(src);

    int xrat = fine.size() / crse.size();

    for (int i=0; i<crse.size(); i++)
      crse[i] = fine[xrat*i];
  }

};


//
// main
//

template<typename scalar, typename argsT, typename controllerT, typename buildT>
void auto_add(controllerT c, vector<pair<int,string>> nodes, vector<argsT> args, buildT build) {
  for (int l=0; l<nodes.size(); l++) {
    auto  nds = pfasst::compute_nodes<double>(get<0>(nodes[l]), get<1>(nodes[l]));
    tuple<pfasst::ISweeper*,pfasst::ITransfer*,pfasst::encap::EncapsulationFactory<scalar>*> t = build(args[l]);
  }
}

int main(int argc, char **argv)
{
  int ndofs  = 256;
  int nnodes = 5;

  if (nlevs == 1) {
    pfasst::SDC<double> sdc;

    auto  nodes   = pfasst::compute_nodes<double>(nnodes, "gauss-lobatto");
    auto* factory = new pfasst::encap::VectorFactory<scalar,double>(ndofs);
    auto* sweeper = new ADIMEX<scalar>(ndofs);

    sweeper->set_nodes(nodes);
    sweeper->set_factory(factory);

    sdc.add_level(sweeper);
    sdc.set_duration(dt, nsteps, 4);
    sdc.setup();

    dvector q0(ndofs);
    sweeper->exact(q0, 0.0);
    sweeper->set_q(&q0, 0);

    sdc.run();
  } else {
    pfasst::MLSDC<double> mlsdc;

    // vector<pair<int,string>> nodes = { { 3, "gauss-lobatto" }, { 5, "gauss-lobatto" } };
    // vector<int>              ndofs = { 256, 512 };

    auto builder = [] (int nx) {
	auto* factory  = new pfasst::encap::VectorFactory<scalar,double>(nx);
	auto* sweeper  = new ADIMEX<scalar>(nx);
	auto* transfer = new ADTRANS<scalar,double>();
	return make_tuple(sweeper, transfer, factory);
    };

    // auto_add<scalar>(mlsdc, nodes, ndofs, builder);

    for (int l=0; l<nlevs; l++) {
      auto  nodes    = pfasst::compute_nodes<double>(nnodes, "gauss-lobatto");
      auto* factory  = new pfasst::encap::VectorFactory<scalar,double>(ndofs);
      auto* sweeper  = new ADIMEX<scalar>(ndofs);
      auto* transfer = new ADTRANS<scalar,double>();

      sweeper->set_nodes(nodes);
      sweeper->set_factory(factory);

      ndofs  = ndofs / 2;
      nnodes = (nnodes-1)/2 + 1;

      mlsdc.add_level(sweeper, transfer);
    }

    mlsdc.set_duration(dt, nsteps, 4);
    mlsdc.setup();

    for (int l=0; l<nlevs; l++) {
      auto* sweeper = mlsdc.get_level<ADIMEX<scalar>>(l);
      dvector& q0 = *dynamic_cast<dvector*>(sweeper->get_q(0));
      sweeper->exact(q0, 0.0);
    }

    mlsdc.run();
  }

  return 0;
}
