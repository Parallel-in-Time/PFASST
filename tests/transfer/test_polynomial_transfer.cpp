#include "fixtures/test_helpers.hpp"

#include <pfasst/transfer/polynomial.hpp>
using pfasst::PolynomialTransfer;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
typedef pfasst::vector_encap_traits<double, double> VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait> VectorEncapsulation;

// #include "sweeper/mocks.hpp"
// typedef SweeperMock<double, VectorEncapsulation> Sweeper;
#include "pfasst/sweeper/interface.hpp"
typedef pfasst::Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> Sweeper;

typedef ::testing::Types<PolynomialTransfer<pfasst::transfer_traits<Sweeper, Sweeper>>> PolynomialTransferTypes;
INSTANTIATE_TYPED_TEST_CASE_P(PolynomialTransfer, Concepts, PolynomialTransferTypes);


TEST_MAIN()
