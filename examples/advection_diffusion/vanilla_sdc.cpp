/*
 * Advection/diffusion example using an encapsulated IMEX sweeper.
 *
 * This example uses a vanilla SDC sweeper.
 */

#include <cstdlib>
#include <memory>

#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/config.hpp>
#include <pfasst/sdc.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"

namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      error_map run_vanilla_sdc(double abs_residual_tol)
      {
        SDC<> sdc;

        const size_t nsteps = config::get_value<size_t>("num_steps", 4);
        const double dt     = config::get_value<double>("delta_step", 0.01);
        const size_t nnodes = config::get_value<size_t>("num_nodes", 3);
        const size_t ndofs  = config::get_value<size_t>("spatial_dofs", 64);
        const size_t niters = config::get_value<size_t>("num_iter", 4);
        const quadrature::QuadratureType quad_type = \
          config::get_value<quadrature::QuadratureType>("nodes_type", quadrature::QuadratureType::GaussLegendre);

        auto quad    = quadrature::quadrature_factory(nnodes, quad_type);
        auto factory = make_shared<encap::VectorFactory<double>>(ndofs);
        auto sweeper = make_shared<AdvectionDiffusionSweeper<>>(ndofs);

        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);
        sweeper->set_residual_tolerances(abs_residual_tol, 0.0);

        sdc.add_level(sweeper);
        sdc.set_duration(0.0, nsteps*dt, dt, niters);
        sdc.setup();

        auto q0 = sweeper->get_start_state();
        sweeper->exact(q0, 0.0);

        sdc.run();

        fftw_cleanup();

        return sweeper->get_errors();
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::enable_config_options();
  pfasst::init(argc, argv);

  pfasst::examples::advection_diffusion::run_vanilla_sdc(0.0);
}
#endif
