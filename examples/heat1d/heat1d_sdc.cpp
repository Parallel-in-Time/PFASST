#include <memory>
using namespace std;

#include <pfasst.hpp>
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

        auto sweeper = make_shared<SweeperType>(32);
        sweeper->quadrature() = quadrature_factory<double>(3, QuadratureType::GaussRadau);

        sdc.add_sweeper(sweeper);

        sdc.status()->time() = 0.0;
        sdc.status()->dt() = 0.05;
        sdc.status()->t_end() = 0.05;
        sdc.status()->max_iterations() = 2;

        sdc.setup();

        sweeper->initial_state() = sweeper->exact(sdc.get_status()->get_time());

        sdc.run();
      }
    }  // ::pfasst::examples::advec_diff
  } // ::pfasst::examples
}  // ::pfasst


int main(int argc, char** argv)
{
  pfasst::init(argc, argv);

  pfasst::examples::heat1d::run();
}
