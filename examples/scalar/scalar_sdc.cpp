/**
 * Scalar test equation example using an encapsulated IMEX sweeper.
 *
 * @ingroup ScalarFiles
 * @file examples/scalar/scalar_sdc.cpp
 * @since v0.1.0
 */

#include <complex>
using namespace std;

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/encap/vector.hpp>

#include "scalar_sweeper.hpp"

namespace pfasst
{
  namespace examples
  {
    namespace scalar
    {
      /**
       * Scalar test equation example using an encapsulated IMEX sweeper.
       *
       * Solves Dahlquist's test equation with single level SDC
       *
       * \\( u' = \\lambda*u \\quad\\text{ , } u(0) = u_0 \\)
       *
       * with complex lambda, treating the real part implicitly and the imaginary part
       * explicitly.
       *
       * @ingroup Scalar
       */
      double run_scalar_sdc(const size_t nsteps, const double dt, const size_t nnodes,
                            const size_t niters, const complex<double> lambda,
                            const quadrature::QuadratureType nodetype)
      {
        SDC<> sdc;

        // For test equation, set initial value to \\( 1+0i \\)
        const complex<double> y0 = complex<double>(1.0, 0.0);

        auto quad = quadrature::quadrature_factory(nnodes, nodetype);

        // This is a scalar example, so we use the encap::VectorFactory with fixed length of 1 and
        //  complex type.
        auto factory = make_shared<encap::VectorFactory<complex<double>>>(1);
        auto sweeper = make_shared<ScalarSweeper<>>(lambda, y0);

        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);

        sdc.add_level(sweeper);

        // Final time Tend = dt*nsteps
        sdc.set_duration(0.0, dt*nsteps, dt, niters);
        sdc.setup();

        auto q0 = sweeper->get_start_state();
        sweeper->exact(q0, 0.0);

        sdc.run();

        return sweeper->get_errors();
      }
    }  // ::pfasst::examples::scalar
  }  // ::pfasst::examples
}  // ::pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  const size_t nsteps = 2;
  const double dt     = 1.0;
  const size_t nnodes = 4;
  const size_t niters = 6;
  const complex<double> lambda = complex<double>(-1.0, 1.0);
  const pfasst::quadrature::QuadratureType nodetype = pfasst::quadrature::QuadratureType::GaussLobatto;

  pfasst::init(argc, argv);

  pfasst::examples::scalar::run_scalar_sdc(nsteps, dt, nnodes, niters, lambda, nodetype);
}
#endif
