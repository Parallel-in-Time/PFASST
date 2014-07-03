/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a (serial) multi-level SDC sweeper.
 */

#include <memory>

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/mlsdc.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"

using namespace pfasst;
using namespace pfasst::encap;

int main(int argc, char** argv)
{
  MLSDC<> mlsdc;

  const size_t nlevs  = 2;
  const size_t nsteps = 4;
  const double dt     = 0.01;
  const size_t niters = 4;
  const int    xrat   = 2;
  const int    trat   = 2;

  size_t nnodes = 5;
  size_t ndofs  = 128;

  /*
   * build space/time discretisation levels and add them to mlsdc
   * controller.  this loop adds the finest level first, and
   * subsequently refines in time (accoring to 'trat') and space
   * (according to 'xrat').
   */
  for (size_t l = 0; l < nlevs; l++) {
    auto  nodes    = compute_nodes<double>(nnodes, "gauss-lobatto");
    shared_ptr<VectorFactory<double>> factory  = make_shared<VectorFactory<double>>(ndofs);
    AdvectionDiffusionSweeper<>*      sweeper  = new AdvectionDiffusionSweeper<>(ndofs);
    SpectralTransfer1D<>*             transfer = new SpectralTransfer1D<>();

    sweeper->set_nodes(nodes);
    sweeper->set_factory(factory);

    mlsdc.add_level(sweeper, transfer);

    ndofs  = ndofs / xrat;
    nnodes = (nnodes - 1) / trat + 1;
  }

  /*
   * set up the mlsdc controller (which in turn calls 'setup' on the
   * sweepers that were added above).  this stage typically
   * preallocates various buffers that the sweepers need.
   */
  mlsdc.setup();

  /*
   * set initial conditions on each level
   */
  auto* sweeper = mlsdc.get_level<AdvectionDiffusionSweeper<>>(mlsdc.nlevels() - 1);
  auto q0 = sweeper->get_state(0);
  sweeper->exact(q0, 0.0);

  /*
   * run mlsdc!
   */
  mlsdc.set_duration(dt, nsteps, niters);
  mlsdc.run();

  fftw_cleanup();
}
