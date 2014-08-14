#include <memory>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../examples/boris/boris_sweeper.hpp"
#include "mocks.hpp"

using namespace ::testing;

typedef Position3DEncapsulation<double, double> Position3D;
typedef Velocity3DEncapsulation<double, double> Velocity3D;
typedef Acceleration3DEncapsulation<double, double> Acceleration3D;
typedef Particle3DEncapsulation<double, double> Particle3D;



int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
