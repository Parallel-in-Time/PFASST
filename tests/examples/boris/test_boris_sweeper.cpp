#include <memory>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include "../examples/boris/boris_sweeper.hpp"
using namespace pfasst::examples::boris;

#include "mocks.hpp"


typedef Position3DEncapsulation<double, double> Position3D;
typedef Velocity3DEncapsulation<double, double> Velocity3D;
typedef Acceleration3DEncapsulation<double, double> Acceleration3D;
typedef Particle3DEncapsulation<double, double> Particle3D;



int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
