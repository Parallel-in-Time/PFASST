#include "fixtures/test_helpers.hpp"

#include <pfasst/transfer/spectral_1d.hpp>
using pfasst::Spectral1DTransfer;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
typedef pfasst::vector_encap_traits<double, double> VectorEncapTrait;

#include "pfasst/sweeper/interface.hpp"
typedef pfasst::Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> Sweeper;

typedef ::testing::Types<Spectral1DTransfer<pfasst::transfer_traits<Sweeper, Sweeper>>> Spectral1DTransferTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Spectral1DTransfer, Concepts, Spectral1DTransferTypes);


class Interpolation
  : public ::testing::Test
{
  protected:
    typedef Spectral1DTransfer<pfasst::transfer_traits<Sweeper, Sweeper>> transfer_type;
    typedef pfasst::encap::VectorEncapsulation<double, double>            encap_type;

    transfer_type          transfer;
    shared_ptr<Sweeper>    coarse_sweeper;
    shared_ptr<Sweeper>    fine_sweeper;
    shared_ptr<encap_type> coarse_encap;
    shared_ptr<encap_type> fine_encap;

    virtual void SetUp()
    {
      this->coarse_encap = make_shared<encap_type>(vector<double>(3, 1.0));
      this->fine_encap = make_shared<encap_type>(6);
    }
};

TEST_F(Interpolation, interpolate_constant)
{
  transfer.interpolate_data(this->coarse_encap, this->fine_encap);
  EXPECT_THAT(this->fine_encap->data(), Pointwise(Eq(), vector<double>(6, 1.0)));
}


TEST_MAIN()
