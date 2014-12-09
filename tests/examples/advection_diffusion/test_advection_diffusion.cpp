/*
 * Tests for the advection-diffusion examples.
 */

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#define PFASST_UNIT_TESTING
#include "../examples/advection_diffusion/vanilla_sdc.cpp"
#include "../examples/advection_diffusion/serial_mlsdc.cpp"
#undef PFASST_UNIT_TESTING

using namespace std;


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

  {
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
    ASSERT_EQ(max_iter, (size_t)3);
  }

  {
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
    ASSERT_EQ(max_iter, (size_t)2);
  }
}

TEST(ErrorTest, SerialMLSDC)
{
  typedef error_map::value_type vtype;

  auto errors = run_serial_mlsdc();
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

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
