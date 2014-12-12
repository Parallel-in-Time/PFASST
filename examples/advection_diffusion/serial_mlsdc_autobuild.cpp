/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a (serial) multi-level SDC sweeper.  It is
 * functionally exactly the same as ex2.cpp, but uses the 'auto
 * builder' to shorten the build and setup stages of the MLSDC
 * controller.
 */

#include <cassert>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
using namespace std;

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/mlsdc.hpp>
#include <pfasst/encap/automagic.hpp>
#include <pfasst/encap/vector.hpp>
using namespace pfasst;
using namespace pfasst::encap;

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"
using pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper;
using pfasst::examples::advection_diffusion::SpectralTransfer1D;


int main(int /*argc*/, char** /*argv*/)
{
  MLSDC<> mlsdc;

  const size_t nsteps = 4;
  const double dt     = 0.01;
  const size_t niters = 4;

  vector<pair<size_t, pfasst::quadrature::QuadratureType>> nodes = {
    { 3, pfasst::quadrature::QuadratureType::GaussLobatto },
    { 5, pfasst::quadrature::QuadratureType::GaussLobatto }
  };

  vector<size_t> ndofs = { 64, 128 };

  /*
   * the 'build' function is called once for each level, and returns a
   * tuple containing a sweeper, encapsulation factory, and transfer
   * routines.  in this case our builder is a lambda function that
   * captures the 'ndofs' variable from above.
   */
  auto build_level = [ndofs](size_t level) {
    auto factory  = make_shared<VectorFactory<double>>(ndofs[level]);
    auto sweeper  = make_shared<AdvectionDiffusionSweeper<>>(ndofs[level]);
    auto transfer = make_shared<SpectralTransfer1D<>>();

    return AutoBuildTuple<>(sweeper, transfer, factory);
  };

  /*
   * the 'initial' function is called once for each level to set the
   * intial conditions.
   */
  auto initial = [](shared_ptr<EncapSweeper<>> sweeper, shared_ptr<Encapsulation<>> q0) {
    auto ad = dynamic_pointer_cast<AdvectionDiffusionSweeper<>>(sweeper);
    assert(ad);
    ad->exact(q0, 0.0);
  };

  auto_build(mlsdc, nodes, build_level);
  auto_setup(mlsdc, initial);
  mlsdc.set_duration(0.0, nsteps*dt, dt, niters);
  mlsdc.run();

  fftw_cleanup();
}
