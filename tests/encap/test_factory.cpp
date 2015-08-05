#include "fixtures/test_helpers.hpp"

#include <vector>
using namespace std;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/encapsulation.hpp>
#include <pfasst/encap/vector.hpp>
typedef pfasst::vector_encap_traits<double, double> VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait> VectorEncapsulation;
typedef pfasst::encap::EncapsulationFactory<VectorEncapTrait> VectorEncapsulationFactory;

typedef ::testing::Types<VectorEncapsulationFactory> FactoryTypes;
INSTANTIATE_TYPED_TEST_CASE_P(VectorEncapFactory, Concepts, FactoryTypes);


TEST(VectorEncapFactory, set_size_after_initialization)
{
  VectorEncapsulationFactory factory;
  EXPECT_THAT(factory.size(), Eq(0));

  factory.set_size(3);
  EXPECT_THAT(factory.size(), Eq(3));
}

TEST(VectorEncapFactory, takes_fixed_size)
{
  VectorEncapsulationFactory factory(3);
  EXPECT_THAT(factory.size(), Eq(3));
}

TEST(VectorEncapFactory, produces_encapsulated_vectors)
{
  VectorEncapsulationFactory factory(3);
  auto encap = factory.create();
  EXPECT_THAT(encap, NotNull());
  EXPECT_THAT(encap->get_data(), SizeIs(3));
}

TEST_MAIN()
