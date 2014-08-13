/*
 * Tests for the advection-diffusion examples.
 */

#include <algorithm>
#include <iostream>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace ::testing;

#include <mpi.h>

#define PFASST_UNIT_TESTING
#include "../examples/advection_diffusion/mpi_pfasst.cpp"
#undef PFASST_UNIT_TESTING


using namespace std;

TEST(ErrorTest, MPIPFASST)
{
  typedef error_map::value_type vtype;

  auto errors = run_mpi_pfasst();
  auto get_step  = [](const vtype x) { return get<0>(get<0>(x)); };
  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> ex = { 1.1168e-13, 4.8849e-13, 5.3268e-13, 2.3059e-12 };
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      EXPECT_NEAR(get_error(x), ex[get_step(x)], 1e-14);
    }
  }
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  int rank = 0;
  int result = 1;  // GTest return value 1 (failure), 0 (success)
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    result = RUN_ALL_TESTS();
  }
  MPI_Finalize();
  return result;
}
