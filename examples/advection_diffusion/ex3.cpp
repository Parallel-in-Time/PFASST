/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a (serial) multi-level SDC sweeper.  It is
 * functionally exactly the same as ex2.cpp, but uses the 'auto
 * builder' to shorten the build and setup stages of the MLSDC
 * controller.
 */

#include <tuple>

#include <pfasst.hpp>
#include <pfasst/encap/automagic.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"

using namespace std;
using namespace pfasst;
using namespace pfasst::encap;

int main(int argc, char **argv)
{
  MLSDC<double> mlsdc;

  const int    nsteps = 4;
  const double dt     = 0.01;
  const int    niters = 4;

  vector<pair<int,string>> nodes = {
    { 3, "gauss-lobatto" },
    { 5, "gauss-lobatto" }
  };

  vector<int> ndofs = { 64, 128 };

  /*
   * the 'build' function is called once for each level, and returns a
   * tuple containing a sweeper, encapsulation factory, and transfer
   * routines.  in this case our builder is a lambda function that
   * captures the 'ndofs' variable from above.
   */
  auto build_level = [ndofs] (unsigned int level) {
    auto* factory  = new VectorFactory<double,double>(ndofs[level]);
    auto* sweeper  = new AdvectionDiffusionSweeper<double,double>(ndofs[level]);
    auto* transfer = new SpectralTransfer1D<double,double>();

    return auto_build_tuple<double,double>(sweeper,transfer,factory);
  };

  /*
   * the 'initial' function is called once for each level to set the
   * intial conditions.
   */
  auto initial = [] (EncapSweeper<double,double> *sweeper, Encapsulation<double,double> *q0) {
    auto* ad = dynamic_cast<AdvectionDiffusionSweeper<double,double>*>(sweeper);
    ad->exact(q0, 0.0);
  };

  auto_build<double,double>(mlsdc, nodes, build_level);
  auto_setup<double,double>(mlsdc, initial);
  mlsdc.set_duration(dt, nsteps, niters);
  mlsdc.run();

  return 0;
}
