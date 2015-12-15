/**
 * Advection-Diffusion with MPI-enabled PFASST.
 *
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/mpi_pfasst.cpp
 * @since v0.2.0
 */

#include <cassert>
#include <memory>
#include <vector>
#include <utility>
using namespace std;

#include <mpi.h>
#include <fftw3.h>

#include <pfasst.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/pfasst.hpp>
#include <pfasst/mpi_communicator.hpp>
#include <pfasst/encap/vector.hpp>

#include "advection_diffusion_sweeper.hpp"
#include "spectral_transfer_1d.hpp"

using namespace pfasst::encap;
using namespace pfasst::mpi;

namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Advection/diffusion example using an encapsulated IMEX sweeper.
       *
       * This example uses MPI PFASST.
       *
       * @ingroup AdvectionDiffusion
       */
      error_map run_mpi_pfasst(const double abs_res_tol, const double rel_res_tol,
                               const size_t niters, const size_t nsteps, const double dt,
                               const size_t ndofs_f, const size_t ndofs_c,
                               const size_t nnodes_f, const size_t nnodes_c)
      {
        ML_CLOG(INFO, "Advec", "abs_res_tol: " << abs_res_tol << ", "
                               << "rel_res_tol: " << rel_res_tol << ", "
                               << "niter: " << niters << ", "
                               << "nsteps: " << nsteps << ", "
                               << "dt: " << dt << ", "
                               << "ndofs (f-c): " << ndofs_f << "-" << ndofs_c << ", "
                               << "nnodes (f-c): " << nnodes_f << "-" << nnodes_c);

        MPICommunicator comm(MPI_COMM_WORLD);
        PFASST<> pf;

        auto quad_c     = quadrature::quadrature_factory(nnodes_c, quadrature::QuadratureType::GaussLobatto);
        auto factory_c  = make_shared<VectorFactory<double>>(ndofs_c);
        auto sweeper_c  = make_shared<AdvectionDiffusionSweeper<>>(ndofs_c);
        auto transfer_c = make_shared<SpectralTransfer1D<>>();

        sweeper_c->set_quadrature(quad_c);
        sweeper_c->set_factory(factory_c);
        sweeper_c->set_residual_tolerances(abs_res_tol, rel_res_tol);

        auto quad_f     = quadrature::quadrature_factory(nnodes_f, quadrature::QuadratureType::GaussLobatto);
        auto factory_f  = make_shared<VectorFactory<double>>(ndofs_f);
        auto sweeper_f  = make_shared<AdvectionDiffusionSweeper<>>(ndofs_f);
        auto transfer_f = make_shared<SpectralTransfer1D<>>();

        sweeper_f->set_quadrature(quad_f);
        sweeper_f->set_factory(factory_f);
        sweeper_f->set_residual_tolerances(abs_res_tol, rel_res_tol);

        ML_LOG(INFO, "expected quadrature error: " << quad_c->expected_error() << " (" << nnodes_c << ")");
        ML_LOG(INFO, "expected quadrature error: " << quad_f->expected_error() << " (" << nnodes_f << ")");

        pf.add_level(sweeper_f, transfer_f);
        pf.add_level(sweeper_c, transfer_c);
        pf.setup();

        auto q0 = sweeper_f->get_start_state();
        sweeper_f->exact(q0, 0.0);

        pf.set_comm(&comm);
        pf.set_duration(0.0, nsteps * dt, dt, niters);
        pf.set_nsweeps({2, 1});
        pf.get_finest<AdvectionDiffusionSweeper<>>()->set_residual_tolerances(abs_res_tol, rel_res_tol);
        pf.set_options();
        pf.run();

        auto fine = pf.get_finest<AdvectionDiffusionSweeper<>>();
        return fine->get_errors();
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst


#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  pfasst::init(argc, argv,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_opts,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_logs);

  const double tend        = pfasst::config::get_value<double>("tend", 0.04);
  const double dt          = pfasst::config::get_value<double>("dt", 0.01);
  const size_t nnodes_f    = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t ndofs_f     = pfasst::config::get_value<size_t>("spatial_dofs", 128);
  const size_t niters      = pfasst::config::get_value<size_t>("num_iter", 4);
  const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 0.0);
  const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);

  const size_t nsteps = tend / dt;
  const size_t nnodes_c = (nnodes_f + 1) / 2;
  const size_t ndofs_c = ndofs_f / 2;

  pfasst::examples::advection_diffusion::run_mpi_pfasst(abs_res_tol, rel_res_tol,
                                                        niters, nsteps, dt,
                                                        ndofs_f, ndofs_c, nnodes_f, nnodes_c);
  fftw_cleanup();
  MPI_Finalize();
}
#endif
