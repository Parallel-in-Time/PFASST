/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a vanilla SDC sweeper.
 */

#include <memory>

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/sdc.hpp>

using pfasst::SDC;
using pfasst::sdc_factory;

#include "advection_diffusion_sweeper.hpp"

int main(int /*argc*/, char** /*argv*/)
{
  typedef double time_precision;

  const size_t         nsteps        = 4;
  const time_precision dt            = 0.01;
  const size_t         nnodes        = 5;
  const size_t         ndofs         = 64;
  const size_t         niters        = 4;
  const time_precision initial_value = 0.0;

  auto sdc = sdc_factory<AdvectionDiffusionSweeper<time_precision>,
                         pfasst::encap::VectorFactory<time_precision>,
                         time_precision>(niters, nsteps, dt, ndofs, nnodes, "gauss-lobatto", 
                                         initial_value);

  sdc.run();

  fftw_cleanup();
}
