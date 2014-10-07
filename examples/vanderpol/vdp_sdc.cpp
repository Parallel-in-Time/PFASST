/**
 * Solve the van der Pol oscillator
 *
 */

#include<pfasst.hpp>
#include<pfasst/sdc.hpp>
#include<pfasst/encap/vector.hpp>

#include "vdp_sweeper.hpp"

double run_vdp_sdc(const size_t nsteps, const double dt, const size_t nnodes,
                      const size_t niters, const double nu,
                      const pfasst::QuadratureType nodetype)
{
  pfasst::SDC<> sdc;

  auto nodes = pfasst::compute_nodes(nnodes, nodetype);

  // van der Pol oscillator (as first order system) has two components
  auto factory = make_shared<pfasst::encap::VectorFactory<double>>(2);
  auto sweeper = make_shared<VdpSweeper<>>(nu);

  sweeper->set_nodes(nodes);
  sweeper->set_factory(factory);

  sdc.add_level(sweeper);

  // Final time Tend = dt*nsteps
  sdc.set_duration(0.0, dt*nsteps, dt, niters);
  sdc.setup();

  auto q0 = sweeper->get_state(0);
  sweeper->exact(q0, 0.0);

  sdc.run();

  return sweeper->get_errors();
}

#ifndef PFASST_UNIT_TESTING
/**
 * Main routine running the scalar example with a preset parameters
 */
int main(int /*argc*/, char** /*argv*/)
{
  const size_t nsteps = 2;
  const double dt     = 1.0;
  const size_t nnodes = 4;
  const size_t niters = 6;
  const double nu     = 0.1;
  const pfasst::QuadratureType nodetype = pfasst::QuadratureType::GaussLobatto;

  run_vdp_sdc(nsteps, dt, nnodes, niters, nu, nodetype);
}
#endif
