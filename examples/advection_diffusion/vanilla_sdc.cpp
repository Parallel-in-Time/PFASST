/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a vanilla SDC sweeper.
 */

#include <cstdlib>
#include <memory>

#include <fftw3.h>

#include <pfasst.hpp>
#include "pfasst/config.hpp"
#include <pfasst/sdc.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"

error_map run_vanilla_sdc()
{
  pfasst::SDC<> sdc;

  const size_t nsteps = pfasst::config::get_value<size_t>("num_steps", 4);
  const double dt     = pfasst::config::get_value<double>("delta_step", 0.01);
  const size_t nnodes = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t ndofs  = pfasst::config::get_value<size_t>("spatial_dofs", 64);
  const size_t niters = pfasst::config::get_value<size_t>("num_iter", 4);
  const pfasst::QuadratureType quad_type = \
    pfasst::config::get_value<pfasst::QuadratureType>("nodes_type", pfasst::QuadratureType::GaussLobatto);

  auto nodes   = pfasst::compute_nodes(nnodes, quad_type);
  auto factory = make_shared<pfasst::encap::VectorFactory<double>>(ndofs);
  auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);

  sweeper->set_nodes(nodes);
  sweeper->set_factory(factory);

  sdc.add_level(sweeper);
  sdc.set_duration(0.0, nsteps*dt, dt, niters);
  sdc.setup();

  auto q0 = sweeper->get_state(0);
  sweeper->exact(q0, 0.0);

  sdc.run();

  fftw_cleanup();

  return sweeper->get_errors();
}


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  // First we want to enable command line options for the Advection-Diffusion Sweeper ...
  AdvectionDiffusionSweeper<>::enable_config_options();
  // ... then we initialize all options other default options and parse given parameters ...
  pfasst::init(argc, argv);

  run_vanilla_sdc();
}
#endif
