/*
 * Scalar test equation example using an encapsulated IMEX sweeper.
 * Solves Dahlquist's test equation
 * y' = a y + i b y
 * treating the real part implicitly and the imaginary part explicitly.
 */
 
#include<complex>

#include<pfasst.hpp>
#include<pfasst/sdc.hpp>
#include<pfasst/encap/vector.hpp>

#include "scalar_sweeper.hpp"

double run_scalar_sdc()
{
  pfasst::SDC<> sdc;

  const size_t nsteps = 2;
  const double dt     = 1.0;
  const size_t nnodes = 4;
  const size_t niters = 6;
  const complex<double> lambda = complex<double>(-1.0, 1.0);
  const complex<double> y0     = complex<double>(1.0, 0.0);

  // so far, only lobatto nodes appear to be working
  auto nodes   = pfasst::compute_nodes(nnodes, "gauss-lobatto");
  auto factory = make_shared<pfasst::encap::VectorFactory<complex<double>>>(1);
  auto sweeper = make_shared<ScalarSweeper<>>(lambda, y0);

  sweeper->set_nodes(nodes);
  sweeper->set_factory(factory);

  sdc.add_level(sweeper);
  sdc.set_duration(0.0, dt * nsteps, dt, niters);
  sdc.setup();

  auto q0 = sweeper->get_state(0);
  sweeper->exact(q0, 0.0);

  sdc.run();

  return sweeper->get_errors();
}

#ifndef PFASST_UNIT_TESTING
int main(int /*argc*/, char** /*argv*/)
{
  run_scalar_sdc();
}
#endif
