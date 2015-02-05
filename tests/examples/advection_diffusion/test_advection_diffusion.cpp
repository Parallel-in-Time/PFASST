/*
 * Tests for the advection-diffusion examples.
 */

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#define PFASST_UNIT_TESTING
#include "../examples/advection_diffusion/vanilla_sdc.cpp"
#include "../examples/advection_diffusion/serial_mlsdc.cpp"
#undef PFASST_UNIT_TESTING
using namespace pfasst::examples::advection_diffusion;


MATCHER(DoubleNear, "")
{
  return abs(get<0>(arg) - get<1>(arg)) < 1e-15;
}

MATCHER(DoubleLess, "")
{
  return get<0>(arg) < get<1>(arg);
}

TEST(ErrorTest, VanillaSDC)
{
  typedef error_map::value_type vtype;

  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto errors = run_vanilla_sdc(0.0);
  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
                                             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> tol = { 7e-9, 7e-9, 7e-9, 7e-9 };
  vector<double> err;
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      err.push_back(get_error(x));
    }
  }

  EXPECT_THAT(err, testing::Pointwise(DoubleLess(), tol));
  ASSERT_EQ(max_iter, (size_t) 3);
}

TEST(AdaptiveErrorTest, VanillaSDC)
{
  typedef error_map::value_type vtype;

  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto errors = run_vanilla_sdc(1.e-6);
  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
                                             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> tol = { 5e-8, 5e-8, 5e-8, 5e-8 };
  vector<double> err;
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      err.push_back(get_error(x));
    }
  }

  EXPECT_THAT(err, testing::Pointwise(DoubleLess(), tol));
  ASSERT_EQ(max_iter, (size_t) 2);
}

TEST(RelativeAdaptiveErrorTest, VanillaSDC)
{
  typedef error_map::value_type vtype;

  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto errors = run_vanilla_sdc(0.0, 1.e-6);
  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
                                             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> tol = { 5e-8, 5e-8, 5e-8, 5e-8 };
  vector<double> err;
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      err.push_back(get_error(x));
    }
  }

  EXPECT_THAT(err, testing::Pointwise(DoubleLess(), tol));
  ASSERT_EQ(max_iter, (size_t) 2);
}

TEST(ErrorTest, SerialMLSDC)
{
  typedef error_map::value_type vtype;

  auto errors_and_residuals = run_serial_mlsdc(2);
  auto errors = get<0>(errors_and_residuals);
  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> tol = { 8e-10, 8e-10, 8e-10, 8e-10 };
  vector<double> err;
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      err.push_back(get_error(x));
    }
  }

  EXPECT_THAT(err, testing::Pointwise(DoubleLess(), tol));
}

TEST(FASTest, SerialMLSDC)
{
  typedef error_map::key_type ktype;

  auto errors_and_residuals = run_serial_mlsdc(3);
  auto residuals = get<1>(errors_and_residuals);

  ASSERT_NEAR(residuals[2][ktype(3, 0)], 0.000667207, 1.e-8);
  ASSERT_NEAR(residuals[0][ktype(3, 0)], 6.23966e-07, 1.e-12);
  ASSERT_NEAR(residuals[1][ktype(3, 0)], 1.27783e-08, 1.e-12);
  ASSERT_NEAR(residuals[2][ktype(3, 1)], 6.60607e-07, 1.e-12);
  ASSERT_NEAR(residuals[0][ktype(3, 1)], 5.19702e-10, 1.e-14);
  ASSERT_NEAR(residuals[1][ktype(3, 1)], 2.59963e-10, 1.e-12);
  ASSERT_NEAR(residuals[2][ktype(3, 2)], 8.89424e-09, 1.e-12);
  ASSERT_NEAR(residuals[0][ktype(3, 2)], 8.28716e-11, 1.e-14);
  ASSERT_NEAR(residuals[1][ktype(3, 2)], 4.54949e-11, 1.e-14);
  ASSERT_NEAR(residuals[2][ktype(3, 3)], 1.04101e-10, 1.e-12);
  ASSERT_NEAR(residuals[0][ktype(3, 3)], 8.35953e-11, 1.e-15);
  ASSERT_NEAR(residuals[1][ktype(3, 3)], 4.20877e-11, 1.e-15);
  ASSERT_NEAR(residuals[2][ktype(3, 4)], 2.18056e-12, 1.e-15);
  ASSERT_NEAR(residuals[0][ktype(3, 4)], 8.34365e-11, 1.e-15);
  ASSERT_NEAR(residuals[1][ktype(3, 4)], 4.19699e-11, 1.e-15);
  ASSERT_NEAR(residuals[2][ktype(3, 5)], 7.18701e-13, 1.e-15);
  ASSERT_NEAR(residuals[0][ktype(3, 5)], 8.34336e-11, 1.e-15);
  ASSERT_NEAR(residuals[1][ktype(3, 5)], 4.19691e-11, 1.e-15);
  ASSERT_NEAR(residuals[2][ktype(3, 6)], 7.07797e-13, 1.e-15);
  ASSERT_NEAR(residuals[0][ktype(3, 6)], 8.34340e-11, 1.e-15);
  ASSERT_NEAR(residuals[1][ktype(3, 6)], 4.19693e-11, 1.e-15);
  ASSERT_NEAR(residuals[2][ktype(3, 7)], 7.07356e-13, 1.e-15);
  ASSERT_NEAR(residuals[0][ktype(3, 7)], 8.34338e-11, 1.e-15);
  ASSERT_NEAR(residuals[1][ktype(3, 7)], 4.19698e-11, 1.e-15);
  ASSERT_NEAR(residuals[2][ktype(3, 8)], 7.07458e-13, 1.e-15);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  pfasst::log::start_log(argc, argv);
  pfasst::log::add_custom_logger("Advec");
  return RUN_ALL_TESTS();
}
