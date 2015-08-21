#include <memory>
using namespace std;

#include <pfasst.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/vector.hpp>
#include <pfasst/controller/two_level_mlsdc.hpp>
#include <pfasst/transfer/spectral_1d.hpp>

#include "heat1d_sweeper.hpp"

using pfasst::encap::VectorEncapsulation;
using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::Spectral1DTransfer;
using pfasst::TwoLevelMLSDC;

using pfasst::examples::heat1d::Heat1D;

typedef VectorEncapsulation<double, double>                        EncapType;
typedef Heat1D<pfasst::sweeper_traits<typename EncapType::traits>> SweeperType;
typedef pfasst::transfer_traits<SweeperType, SweeperType, 2>       TransferTraits;
typedef Spectral1DTransfer<TransferTraits>                         TransferType;


namespace pfasst
{
  namespace examples
  {
    namespace heat1d
    {
      void run_mlsdc(const size_t& ndofs, const size_t& coarse_factor, const size_t& nnodes,
                     const QuadratureType& quad_type, const double& t_0, const double& dt,
                     const double& t_end, const size_t& niter)
      {
        TwoLevelMLSDC<TransferType> mlsdc;

        auto coarse = make_shared<SweeperType>(ndofs / coarse_factor);
        coarse->quadrature() = quadrature_factory<double>(nnodes, quad_type);
        auto fine = make_shared<SweeperType>(ndofs);
        fine->quadrature() = quadrature_factory<double>(nnodes, quad_type);

        auto transfer = make_shared<TransferType>();

        mlsdc.add_sweeper(coarse, true);
        mlsdc.add_sweeper(fine, false);
        mlsdc.add_transfer(transfer);
        mlsdc.set_options();

        mlsdc.status()->time() = t_0;
        mlsdc.status()->dt() = dt;
        mlsdc.status()->t_end() = t_end;
        mlsdc.status()->max_iterations() = niter;

        mlsdc.setup();

        coarse->initial_state() = coarse->exact(mlsdc.get_status()->get_time());
        fine->initial_state() = fine->exact(mlsdc.get_status()->get_time());

        mlsdc.run();
        mlsdc.post_run();
      }
    }  // ::pfasst::examples::advec_diff
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  using pfasst::config::get_value;
  using pfasst::quadrature::QuadratureType;

  pfasst::init(argc, argv, SweeperType::init_opts);

  const size_t ndofs = get_value<size_t>("num_dofs", 8);
  const size_t coarse_factor = get_value<size_t>("coarse_factor", 2);
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

  pfasst::examples::heat1d::run_mlsdc(ndofs, coarse_factor, nnodes, quad_type, t_0, dt, t_end, niter);
}
