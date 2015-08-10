#include "fixtures/test_helpers.hpp"

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
typedef pfasst::vector_encap_traits<double, double> VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait> VectorEncapsulation;

#include "examples/heat1d/heat1d_sweeper.hpp"
template<class... Ts>
using Sweeper = pfasst::examples::heat1d::Heat1D<Ts...>;


class Setup
  : public ::testing::Test
{
  protected:
    typedef          Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> sweeper_type;
    typedef typename sweeper_type::encap_type                          encap_type;

    sweeper_type sweeper;

    shared_ptr<pfasst::Status<double>> status = make_shared<pfasst::Status<double>>();

    virtual void SetUp()
    {}
};



TEST_MAIN()
