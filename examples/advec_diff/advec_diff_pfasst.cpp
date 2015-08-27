#include <memory>
using namespace std;

#include <mpi.h>

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/vector.hpp>
#include <pfasst/comm/mpi_p2p.hpp>
#include <pfasst/controller/two_level_pfasst.hpp>
#include <pfasst/transfer/spectral_1d.hpp>

#include "advec_diff_sweeper.hpp"

using pfasst::encap::VectorEncapsulation;
using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::Spectral1DTransfer;
using pfasst::TwoLevelPfasst;
typedef pfasst::comm::MpiP2P CommType;

using pfasst::examples::advec_diff::AdvecDiff;

typedef VectorEncapsulation<double, double>                           EncapType;
typedef AdvecDiff<pfasst::sweeper_traits<typename EncapType::traits>> SweeperType;
typedef pfasst::transfer_traits<SweeperType, SweeperType, 2>          TransferTraits;
typedef Spectral1DTransfer<TransferTraits>                            TransferType;


namespace pfasst
{
  namespace examples
  {
    namespace advec_diff
    {
      void run_pfasst(const size_t& ndofs, const size_t& nnodes, const QuadratureType& quad_type,
                      const double& t_0, const double& dt, const double& t_end, const size_t& niter)
      {
        TwoLevelPfasst<TransferType, CommType> pfasst;
        pfasst.communicator() = make_shared<CommType>(MPI_COMM_WORLD);

        auto coarse = make_shared<SweeperType>(ndofs);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = make_shared<SweeperType>(ndofs * 2);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        auto transfer = make_shared<TransferType>();

        pfasst.add_sweeper(coarse, true);
        pfasst.add_sweeper(fine, false);
        pfasst.add_transfer(transfer);
        pfasst.set_options();

        pfasst.status()->time() = t_0;
        pfasst.status()->dt() = dt;
        pfasst.status()->t_end() = t_end;
        pfasst.status()->max_iterations() = niter;

        pfasst.setup();

        coarse->initial_state() = coarse->exact(pfasst.get_status()->get_time());
        fine->initial_state() = fine->exact(pfasst.get_status()->get_time());

        pfasst.run();
        pfasst.post_run();
      }
    }  // ::pfasst::examples::heat1d
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  MPI_Init(&argc, &argv);

  pfasst::init(argc, argv, SweeperType::init_opts);
  pfasst::Status<double>::create_mpi_datatype();

  const size_t ndofs = get_value<size_t>("num_dofs", 4);
  const size_t nnodes = get_value<size_t>("num_nodes", 3);
  const QuadratureType quad_type = QuadratureType::GaussRadau;
  const double t_0 = 0.0;
  const double dt = get_value<double>("dt", 0.01);
  double t_end = get_value<double>("tend", -1);
  size_t nsteps = get_value<size_t>("num_steps", 0);
  if (t_end == -1 && nsteps == 0) {
    CLOG(ERROR, "USER") << "Either t_end or num_steps must be specified.";
    throw runtime_error("either t_end or num_steps must be specified");
  } else if (t_end != -1 && nsteps != 0) {
    if (!pfasst::almost_equal(t_0 + nsteps * dt, t_end)) {
      CLOG(ERROR, "USER") << "t_0 + nsteps * dt != t_end ("
      << t_0 << " + " << nsteps << " * " << dt << " = " << (t_0 + nsteps * dt)
      << " != " << t_end << ")";
      throw runtime_error("t_0 + nsteps * dt != t_end");
    }
  } else if (nsteps != 0) {
    t_end = t_0 + dt * nsteps;
  }
  const size_t niter = get_value<size_t>("num_iters", 5);

  pfasst::examples::advec_diff::run_pfasst(ndofs, nnodes, quad_type, t_0, dt, t_end, niter);

  pfasst::Status<double>::free_mpi_datatype();

  MPI_Finalize();

  return 0;
}
