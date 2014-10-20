/**
 * Run SDC for the van der Pol oscillator.
 *
 */

#include<pfasst/sdc.hpp>
#include<pfasst/encap/vector.hpp>

#include "vdp_sweeper.hpp"

double run_vdp_sdc(const size_t nsteps, const double dt, const size_t nnodes,
                   const size_t niters, const double nu, const double x0,
                   const double y0, const pfasst::QuadratureType nodetype)
{
  pfasst::SDC<> sdc;
  
  auto nodes = pfasst::compute_nodes(nnodes, nodetype);
  
  // van der Pol oscillator (as first order system) has two components
  auto factory = make_shared<pfasst::encap::VectorFactory<double>>(2);
  // input is parameter nu and initial values for position and velocity
  auto sweeper = make_shared<VdpSweeper<>>(nu, x0, y0);
  
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
  const double x0 = 2.0, y0 = 0.0;
  const size_t nsteps = 1000;
  const double dt     = 50.0/( (double) nsteps);
  const size_t nnodes = 3;
  const size_t niters = 1;
  const double nu     = 5.0;
  const pfasst::QuadratureType nodetype = pfasst::QuadratureType::GaussLegendre;
  
  std::cout << "Used timestep: " << dt << std::endl;
  
  run_vdp_sdc(nsteps, dt, nnodes, niters, nu, x0, y0, nodetype);
}
#endif
