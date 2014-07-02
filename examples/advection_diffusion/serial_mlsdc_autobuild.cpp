/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a (serial) multi-level SDC sweeper.  It is
 * functionally exactly the same as ex2.cpp, but uses the 'auto
 * builder' to shorten the build and setup stages of the MLSDC
 * controller.
 */

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/mlsdc.hpp>
#include <pfasst/encap/automagic.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"

using namespace std;
using namespace pfasst;
using namespace pfasst::encap;

int main(int argc, char **argv)
{
  MLSDC<> mlsdc;

  const size_t nsteps = 4;
  const double dt     = 0.01;
  const size_t niters = 4;

  vector<pair<size_t,string>> nodes = {
    { 3, "gauss-lobatto" },
    { 5, "gauss-lobatto" }
  };

  vector<size_t> ndofs = { 64, 128 };

  /*
   * the 'build' function is called once for each level, and returns a
   * tuple containing a sweeper, encapsulation factory, and transfer
   * routines.  in this case our builder is a lambda function that
   * captures the 'ndofs' variable from above.
   */
  auto build_level = [ndofs] (size_t level) {
    auto* factory  = new VectorFactory<double>(ndofs[level]);
    auto* sweeper  = new AdvectionDiffusionSweeper<>(ndofs[level]);
    auto* transfer = new SpectralTransfer1D<>();

    return auto_build_tuple<>(sweeper,transfer,factory);
  };

  /*
   * the 'initial' function is called once for each level to set the
   * intial conditions.
   */
  auto initial = [] (EncapSweeper<> *sweeper,
		     Encapsulation<> *q0) {
    auto* ad = dynamic_cast<AdvectionDiffusionSweeper<>*>(sweeper);
    ad->exact(q0, 0.0);
  };

  auto_build(mlsdc, nodes, build_level);
  auto_setup(mlsdc, initial);
  mlsdc.set_duration(dt, nsteps, niters);
  mlsdc.run();

  fftw_cleanup();
}
