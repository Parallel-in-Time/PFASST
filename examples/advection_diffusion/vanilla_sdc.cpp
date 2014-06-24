/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a vanilla SDC sweeper.
 */

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/sdc.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"

int main(int argc, char **argv)
{
  pfasst::SDC<double> sdc;

  const int    nsteps = 4;
  const double dt     = 0.01;
  const int    nnodes = 5;
  const int    ndofs  = 64;
  const int    niters = 4;

  auto  nodes   = pfasst::compute_nodes<double>(nnodes, "gauss-lobatto");
  auto* factory = new pfasst::encap::VectorFactory<double,double>(ndofs);
  auto* sweeper = new AdvectionDiffusionSweeper<double>(ndofs);

  sweeper->set_nodes(nodes);
  sweeper->set_factory(factory);

  sdc.add_level(sweeper);
  sdc.set_duration(dt, nsteps, niters);
  sdc.setup();

  auto* q0 = sweeper->get_state(0);
  sweeper->exact(q0, 0.0);

  sdc.run();

  fftw_cleanup();
}
