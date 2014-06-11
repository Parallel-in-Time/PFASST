/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a multi-level SDC sweeper.
 */

#include <pfasst.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"

int main(int argc, char **argv)
{
  pfasst::MLSDC<double> mlsdc;

  const int    nlevs  = 2;
  const int    nsteps = 4;
  const double dt     = 0.01;
  const int    niters = 4;
  const int    xrat   = 2;
  const int    trat   = 2;

  int    nnodes = 5;
  int    ndofs  = 128;

  for (int l=0; l<nlevs; l++) {
    auto  nodes    = pfasst::compute_nodes<double>(nnodes, "gauss-lobatto");
    auto* factory  = new pfasst::encap::VectorFactory<double,double>(ndofs);
    auto* sweeper  = new AdvectionDiffusionSweeper<double,double>(ndofs);
    auto* transfer = new SpectralTransfer1D<double,double>();

    sweeper->set_nodes(nodes);
    sweeper->set_factory(factory);

    ndofs  = ndofs / xrat;
    nnodes = (nnodes-1) / trat + 1;

    mlsdc.add_level(sweeper, transfer);
  }

  mlsdc.set_duration(dt, nsteps, niters);
  mlsdc.setup();

  for (int l=0; l<nlevs; l++) {
    auto* sweeper = mlsdc.get_level<AdvectionDiffusionSweeper<double,double>>(l);
    auto* q0 = sweeper->get_state(0);
    sweeper->exact(q0, 0.0);
  }

  mlsdc.run();

  return 0;
}
