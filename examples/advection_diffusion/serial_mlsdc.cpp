/**
 * Advection-Diffusion with serial MLSDC.
 *
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/serial_mlsdc.cpp
 * @since v0.1.0
 */
#include <memory>
using namespace std;

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>
#include <pfasst/controller/mlsdc.hpp>
#include <pfasst/encap/vector.hpp>
using namespace pfasst::encap;

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"


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
       *
       * @ingroup AdvectionDiffusion
       */
      tuple<error_map, residual_map> run_serial_mlsdc(size_t nlevs,
                                                      size_t nsteps_in=4,
                                                      double step_size_in=0.01,
                                                      size_t num_iter_in=8,
                                                      size_t nnodes_in=5,
                                                      size_t ndofs_in=128)
      {
        MLSDC<> mlsdc;

        const size_t nsteps = config::get_value<size_t>("num_steps", nsteps_in);
        const double dt     = config::get_value<double>("step_size", step_size_in);
        const size_t niters = config::get_value<size_t>("num_iter", num_iter_in);
        const int    xrat   = 2;
        const int    trat   = 2;

        size_t nnodes = config::get_value<size_t>("num_nodes", nnodes_in);
        size_t ndofs  = config::get_value<size_t>("spatial_dofs", ndofs_in);

        const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 0.0);
        const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);

        /*
         * build space/time discretisation levels and add them to mlsdc
         * controller.  this loop adds the finest level first, and
         * subsequently refines in time (accoring to 'trat') and space
         * (according to 'xrat').
         */
        for (size_t l = 0; l < nlevs; l++) {
          auto quad     = quadrature::quadrature_factory(nnodes, quadrature::QuadratureType::GaussLobatto);
          auto factory  = make_shared<VectorFactory<double>>(ndofs);
          auto sweeper  = make_shared<AdvectionDiffusionSweeper<>>(ndofs);
          auto transfer = make_shared<SpectralTransfer1D<>>();

          ML_LOG(INFO, "expected quadrature error: " << quad->expected_error() << " (" << nnodes << ")");

          sweeper->set_quadrature(quad);
          sweeper->set_factory(factory);
          sweeper->set_residual_tolerances(abs_res_tol, rel_res_tol);

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
        //  sweeper->set_residual_tolerances(1e-5, 0.0);

        /*
         * run mlsdc!
         */
        mlsdc.set_duration(0.0, nsteps*dt, dt, niters);
        mlsdc.set_options();
        mlsdc.run();

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
  pfasst::examples::advection_diffusion::run_serial_mlsdc(3);
  fftw_cleanup();
}
#endif
