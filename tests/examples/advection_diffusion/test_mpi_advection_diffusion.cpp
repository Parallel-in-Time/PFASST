/*
 * Tests for the advection-diffusion examples.
 */

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <mpi.h>

#define PFASST_UNIT_TESTING
#include "../examples/advection_diffusion/mpi_pfasst.cpp"
#undef PFASST_UNIT_TESTING


using namespace std;

// MATCHER(DoubleNear, "")
// {
//   return abs(get<0>(arg) - get<1>(arg)) < 1e-15;
// }

// MATCHER(DoubleLess, "")
// {
//   return get<0>(arg) < get<1>(arg);
// }

TEST(ErrorTest, MPIPFASST)
{
  typedef error_map::value_type vtype;

  auto errors = run_mpi_pfasst();
  auto get_step  = [](const vtype x) { return get<0>(get<0>(x)); };
  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> ex = { 1.224103e-10, 5.145808e-10, 3.389905e-9, 1.920198e-7 };
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      EXPECT_NEAR(get_error(x), ex[get_step(x)], 1e-12);
    }
  }
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  auto result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
