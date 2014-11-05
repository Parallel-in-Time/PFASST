#include <memory>
#include <vector>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include "../examples/boris/physics.hpp"

using namespace pfasst;
using namespace pfasst::examples::boris;

#include "mocks.hpp"

typedef MockParticle<double, double> MockParticleT;

typedef ElectricField<double, double, MockParticle> EFieldT;
typedef MagneticField<double, double, MockParticle> BFieldT;
typedef MockEOperator<double, double> EnergyOperatorT;


TEST(ElectricFieldTest, Instantiation)
{
  EFieldT default_ctor;
  EXPECT_THAT(default_ctor.omega_e, DoubleEq(1.0));

  EFieldT special_ctor = EFieldT(0.5);
  EXPECT_THAT(special_ctor.omega_e, DoubleEq(0.5));
}

TEST(ElectricFieldTest, OmegaZ)
{
  EFieldT default_ctor;
  default_ctor.omega_e = 0.0;
  EXPECT_THAT(default_ctor.omega_e, DoubleEq(0.0));
}

TEST(ElectricFieldTest, Evaluation)
{
  EFieldT default_ctor;
  vector<shared_ptr<MockParticleT>> particles { make_shared<MockParticleT>(),
                                                make_shared<MockParticleT>(),
                                                make_shared<MockParticleT>() };
  EXPECT_THROW(default_ctor.evaluate(particles, 0, 0.0), NotImplementedYet);
}


TEST(MagneticFieldTest, Instantiation)
{
  BFieldT default_ctor;
  EXPECT_THAT(default_ctor.omega_b, DoubleEq(1.0));

  BFieldT special_ctor = BFieldT(0.5);
  EXPECT_THAT(special_ctor.omega_b, DoubleEq(0.5));
}

TEST(MagneticFieldTest, OmegaC)
{
  BFieldT default_ctor;
  default_ctor.omega_b = 0.0;
  EXPECT_THAT(default_ctor.omega_b, DoubleEq(0.0));
}

TEST(MagneticFieldTest, Evaluation)
{
  BFieldT default_ctor;
  vector<shared_ptr<MockParticleT>> particles { make_shared<MockParticleT>(),
                                                make_shared<MockParticleT>(),
                                                make_shared<MockParticleT>() };
  EXPECT_THROW(default_ctor.evaluate(particles, 0, 0.0), NotImplementedYet);
}


TEST(EnergyOperatorTest, Evaluation)
{
  EnergyOperatorT e_operator;

  vector<shared_ptr<typename EnergyOperatorT::particle_type>> particles \
    { make_shared<typename EnergyOperatorT::particle_type>(),
      make_shared<typename EnergyOperatorT::particle_type>(),
      make_shared<typename EnergyOperatorT::particle_type>() };

  ON_CALL(e_operator, evaluate(_, _))
    .WillByDefault(Return(1.0));
  EXPECT_CALL(e_operator, evaluate(_, _))
    .Times(1);
  EXPECT_THAT(e_operator.evaluate(particles, 0.0), DoubleEq(1.0));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
