#include<complex>

#include<pfasst.hpp>
#include<pfasst/sdc.hpp>
#include<pfasst/encap/vector.hpp>

#include "scalar_sweeper.hpp"

int main(int, char**)
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

}
