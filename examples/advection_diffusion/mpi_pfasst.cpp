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
#include <pfasst/controller/pfasst.hpp>
#include <pfasst/mpi_communicator.hpp>
#include <pfasst/encap/automagic.hpp>
#include <pfasst/encap/mpi_vector.hpp>

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
      error_map run_mpi_pfasst(double abs_residual_tol, size_t niters=4)
      {
        const size_t nsteps = 4;
        const double dt     = 0.01;

        vector<pair<size_t, quadrature::QuadratureType>> nodes = {
          { 3, quadrature::QuadratureType::GaussLobatto },
          { 5, quadrature::QuadratureType::GaussLobatto }
        };

        vector<size_t> ndofs = { 64, 128 };

        auto build_level = [ndofs](size_t level) {
          auto factory  = make_shared<MPIVectorFactory<double>>(ndofs[level]);
          auto sweeper  = make_shared<AdvectionDiffusionSweeper<>>(ndofs[level]);
          auto transfer = make_shared<SpectralTransfer1D<>>();

          return AutoBuildTuple<>(sweeper, transfer, factory);
        };

        auto initial = [](shared_ptr<EncapSweeper<>> sweeper, shared_ptr<Encapsulation<>> q0) {
          auto ad = dynamic_pointer_cast<AdvectionDiffusionSweeper<>>(sweeper);
          assert(ad);
          ad->exact(q0, 0.0);
        };

        MPICommunicator comm(MPI_COMM_WORLD);
        PFASST<> pf;

        auto_build(pf, nodes, build_level);
        auto_setup(pf, initial);

        pf.set_comm(&comm);
        pf.set_duration(0.0, nsteps * dt, dt, niters);
        pf.set_nsweeps({2, 1});
        pf.get_finest<AdvectionDiffusionSweeper<>>()->set_residual_tolerances(abs_residual_tol, 0.0);
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
  pfasst::examples::advection_diffusion::run_mpi_pfasst(0.0);
  fftw_cleanup();
  MPI_Finalize();
}
#endif
