/**
 * Advection-Diffusion MLSDC with _auto builder_.
 *
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/serial_mlsdc_autobuild.cpp
 * @since v0.1.0
 */
#include <cassert>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
using namespace std;

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>
#include <pfasst/controller/mlsdc.hpp>
#include <pfasst/encap/automagic.hpp>
#include <pfasst/encap/vector.hpp>
using namespace pfasst;
using namespace pfasst::encap;

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"
using pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper;
using pfasst::examples::advection_diffusion::SpectralTransfer1D;


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Advection/diffusion example using an encapsulated IMEX sweeper.
       *
       * This example uses a (serial) multi-level SDC sweeper.
       * It is functionally exactly the same as serial_mlsdc.cpp, but uses the _auto builder_ to
       * shorten the build and setup stages of the MLSDC controller.
       *
       * @ingroup AdvectionDiffusion
       */
      tuple<error_map, residual_map> run_serial_mlsdc_autobuild()
      {
        MLSDC<> mlsdc;

        const size_t nsteps = config::get_value<size_t>("num_steps", 4);
        const double dt     = config::get_value<double>("delta_step", 0.01);
        const size_t niters = config::get_value<size_t>("num_iter", 4);

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

        tuple<error_map, residual_map> rinfo;
        get<0>(rinfo) = mlsdc.get_finest<AdvectionDiffusionSweeper<>>()->get_errors();
        for (auto l = mlsdc.coarsest(); l <= mlsdc.finest(); ++l) {
          get<1>(rinfo).insert(pair<size_t, error_map>(l.level, l.current<AdvectionDiffusionSweeper<>>()->get_residuals()));
        }
        return rinfo;
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::init(argc, argv,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_opts,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_logs);
  pfasst::examples::advection_diffusion::run_serial_mlsdc_autobuild();
}
#endif
