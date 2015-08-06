#include <memory>
using namespace std;

#include <pfasst/logging.hpp>
#include <pfasst/quadrature.hpp>
#include <pfasst/encap/vector.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/transfer/spectral_1d.hpp>

#include "heat1d_sweeper.hpp"

using pfasst::encap::VectorEncapsulation;
using pfasst::quadrature::quadrature_factory;
using pfasst::quadrature::QuadratureType;
using pfasst::Spectral1DTransfer;
using pfasst::SDC;

using pfasst::examples::heat1d::Heat1D;

typedef VectorEncapsulation<double, double>                                      EncapType;
typedef Heat1D<pfasst::sweeper_traits<typename EncapType::traits>>               SweeperType;
typedef Spectral1DTransfer<pfasst::transfer_traits<SweeperType, SweeperType, 1>> TransferType;


namespace pfasst
{
  namespace examples
  {
    namespace heat1d
    {
      void run()
      {
        SDC<TransferType> sdc;

        auto sweeper = make_shared<SweeperType>(8);
        sweeper->quadrature() = quadrature_factory<double>(3, QuadratureType::GaussRadau);

        sdc.add_sweeper(sweeper);

        sdc.status()->dt() = 0.01;
        sdc.status()->t_end() = 0.01;
        sdc.status()->max_iterations() = 2;

        sdc.setup();

        sweeper->initial_state()->zero();
//         fill(sweeper->initial_state()->data().begin(), sweeper->initial_state()->data().end(), 0.5);
//         sweeper->initial_state()->data()[2] = 0.5;
//         sweeper->initial_state()->data()[3] = 1.0;
//         sweeper->initial_state()->data()[4] = 0.5;
        sweeper->initial_state() = sweeper->exact(0.0);

        sdc.run();
      }
    }  // ::pfasst::examples::advec_diff
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  pfasst::log::start_log(argc, argv);

  pfasst::examples::heat1d::run();
}
