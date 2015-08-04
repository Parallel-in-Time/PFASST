#include <memory>
using namespace std;

#include <pfasst/logging.hpp>
#include <pfasst/encap/vector.hpp>
using pfasst::encap::VectorEncapsulation;

#include "advec_diff_sweeper.hpp"
using pfasst::examples::AdvecDiff;

typedef VectorEncapsulation<double, double>                           EncapType;
typedef AdvecDiff<pfasst::sweeper_traits<typename EncapType::traits>> SweeperType;


int main(int argc, char** argv)
{
  pfasst::log::start_log(argc, argv);

  auto sweeper = make_shared<SweeperType>(16);
}
