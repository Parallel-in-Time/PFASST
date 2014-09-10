/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a vanilla SDC sweeper.
 */

#include <cstdlib>
#include <memory>

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/sdc.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"

error_map run_vanilla_sdc()
{
  pfasst::SDC<> sdc;

  const size_t nsteps = 4;
  const double dt     = 0.01;
  const size_t nnodes = 5;
  const size_t ndofs  = 64;
  const size_t niters = 4;

   // auto quad    = pfasst::quadrature::quadrature_factory(nnodes, pfasst::quadrature::QuadratureType::GaussLobatto);
  auto quad    = pfasst::quadrature::quadrature_factory(nnodes-2, pfasst::quadrature::QuadratureType::GaussLegendre);
  auto factory = make_shared<pfasst::encap::VectorFactory<double>>(ndofs);
  auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);

  sweeper->set_quadrature(quad);
  sweeper->set_factory(factory);

  sdc.add_level(sweeper);
  sdc.set_duration(0.0, nsteps*dt, dt, niters);
  sdc.setup();

  auto q0 = sweeper->get_start_state();
  sweeper->exact(q0, 0.0);

  sdc.run();

  fftw_cleanup();

  return sweeper->get_errors();
}


#ifndef PFASST_UNIT_TESTING
int main(int /*argc*/, char** /*argv*/)
{
  run_vanilla_sdc();
}
#endif
