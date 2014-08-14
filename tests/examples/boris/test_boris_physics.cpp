#include <memory>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../examples/boris/physics.hpp"

#include "mocks.hpp"

using namespace ::testing;

typedef MockParticle<double, double> MockParticleT;

typedef ElectricField<double, double, MockParticle> EFieldT;
typedef MagneticField<double, double, MockParticle> BFieldT;
typedef MockEOperator<double, double> EnergyOperatorT;


TEST(ElectricFieldTest, Instantiation)
{
  EFieldT default_ctor;
  EXPECT_THAT(default_ctor.omega_z, DoubleEq(1.0));

  EFieldT special_ctor = EFieldT(0.5);
  EXPECT_THAT(special_ctor.omega_z, DoubleEq(0.5));
}

TEST(ElectricFieldTest, OmegaZ)
{
  EFieldT default_ctor;
  default_ctor.omega_z = 0.0;
  EXPECT_THAT(default_ctor.omega_z, DoubleEq(0.0));
}

TEST(ElectricFieldTest, Evaluation)
{
  EFieldT default_ctor;
  shared_ptr<const MockParticleT> p = make_shared<const MockParticleT>();
  EXPECT_THROW(default_ctor.evaluate(p, 0.0), NotImplementedYet);
}


TEST(MagneticFieldTest, Instantiation)
{
  BFieldT default_ctor;
  EXPECT_THAT(default_ctor.omega_c, DoubleEq(1.0));

  BFieldT special_ctor = BFieldT(0.5);
  EXPECT_THAT(special_ctor.omega_c, DoubleEq(0.5));
}

TEST(MagneticFieldTest, OmegaC)
{
  BFieldT default_ctor;
  default_ctor.omega_c = 0.0;
  EXPECT_THAT(default_ctor.omega_c, DoubleEq(0.0));
}

TEST(MagneticFieldTest, Evaluation)
{
  BFieldT default_ctor;
  shared_ptr<const MockParticleT> p = make_shared<const MockParticleT>();
  EXPECT_THROW(default_ctor.evaluate(p, 0.0), NotImplementedYet);
}


TEST(EnergyOperatorTest, Evaluation)
{
  EnergyOperatorT e_operator;

  shared_ptr<const typename EnergyOperatorT::e_field_type> EField = \
    make_shared<const typename EnergyOperatorT::e_field_type>();
  shared_ptr<const typename EnergyOperatorT::b_field_type> BField = \
    make_shared<const typename EnergyOperatorT::b_field_type>();
  shared_ptr<const typename EnergyOperatorT::particle_type> particle = \
    make_shared<const typename EnergyOperatorT::particle_type>();

  ON_CALL(e_operator, evaluate(_, _, _, _))
    .WillByDefault(Return(1.0));
  EXPECT_CALL(e_operator, evaluate(_, _, _, _))
    .Times(1);
  EXPECT_THAT(e_operator.evaluate(particle, 0.0, EField, BField), DoubleEq(1.0));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}