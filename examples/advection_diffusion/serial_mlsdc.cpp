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

error_map run_serial_mlsdc()
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
    auto quad     = pfasst::quadrature::quadrature_factory(nnodes, pfasst::quadrature::QuadratureType::GaussLobatto);
    auto factory  = make_shared<VectorFactory<double>>(ndofs);
    auto sweeper  = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
    auto transfer = make_shared<SpectralTransfer1D<>>();

    sweeper->set_quadrature(quad);
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
  auto sweeper = mlsdc.get_finest<AdvectionDiffusionSweeper<>>();
  auto q0 = sweeper->get_start_state();
  sweeper->exact(q0, 0.0);

  /*
   * run mlsdc!
   */
  mlsdc.set_duration(0.0, nsteps*dt, dt, niters);
  mlsdc.run();

  fftw_cleanup();

  return sweeper->get_errors();
}

#ifndef PFASST_UNIT_TESTING
int main(int /*argc*/, char** /*argv*/)
{
  run_serial_mlsdc();
}
#endif
