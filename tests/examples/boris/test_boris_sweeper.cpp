#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include "pfasst/logging.hpp"

#define PFASST_UNIT_TESTING
#include "../examples/boris/boris_sweeper.hpp"
#include "../examples/boris/boris_sdc.cpp"
#undef PFASST_UNIT_TESTING
using namespace pfasst::examples::boris;


TEST(EnergyDriftAndResidual, SingleStep)
{
  const size_t num_iter = 9;
  auto errors_map = run_boris_sdc<double>(1, 0.015625, 5, 1, num_iter+1, 0.0, 0.0);
  ASSERT_THAT(errors_map, SizeIs(num_iter+1));

  auto final_error = errors_map.rbegin()->second;

  EXPECT_THAT(final_error.e_drift, DoubleNear(0.0, 2e-12));
  EXPECT_THAT(final_error.res, DoubleNear(0.0, 1.5e-14));
}

TEST(EnergyDriftAndResidual, MultiStep)
{
  const size_t num_iter = 9;
  const size_t num_steps = 10;
  auto errors_map = run_boris_sdc<double>(num_steps, 0.015625, 5, 1, num_iter+1, 0.0, 0.0);
  ASSERT_THAT(errors_map, SizeIs((num_iter+1) * num_steps));

  auto final_error = errors_map.rbegin()->second;

  EXPECT_THAT(final_error.e_drift, DoubleNear(0.0, 1.1e-11));
  EXPECT_THAT(final_error.res, DoubleNear(0.0, 1.5e-14));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  pfasst::log::start_log(argc, argv);
  return RUN_ALL_TESTS();
}
