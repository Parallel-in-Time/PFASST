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
using namespace ::testing;

#include <mpi.h>

#define PFASST_UNIT_TESTING
#include "../examples/advection_diffusion/mpi_pfasst.cpp"
#undef PFASST_UNIT_TESTING
using namespace pfasst::examples::advection_diffusion;


TEST(ErrorTest, MPIPFASST)
{
  typedef error_map::value_type vtype;

  auto errors = run_mpi_pfasst(0.0, 0.0, 4, 4, 0.01, 128, 64, 5, 3);
  auto get_step  = [](const vtype x) { return get<0>(get<0>(x)); };
  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> ub = { 1.e-12, 1.e-12, 2.5e-12, 5.e-12 };
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      EXPECT_LE(get_error(x), ub[get_step(x)]);
    }
  }
}

TEST(AdaptiveErrorTest, MPIPFASST)
{
  typedef error_map::value_type vtype;

  auto errors = run_mpi_pfasst(1.e-8, 0.0, 12, 4, 0.01, 128, 64, 5, 3);
  auto get_step  = [](const vtype x) { return get<0>(get<0>(x)); };
  auto get_iter  = [](const vtype x) { return get<1>(get<0>(x)); };
  auto get_error = [](const vtype x) { return get<1>(x); };

  auto max_iter = get_iter(*std::max_element(errors.begin(), errors.end(),
             [get_iter](const vtype p1, const vtype p2) { return get_iter(p1) < get_iter(p2); }));

  vector<double> ub = { 5e-8, 5e-8, 5e-8, 5e-8 };
  for (auto& x: errors) {
    if (get_iter(x) == max_iter) {
      EXPECT_LE(get_error(x), ub[get_step(x)]);
    }
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  vector<size_t> ic = { 1, 1, 2, 2 };
  ASSERT_EQ(max_iter, (size_t) ic[rank]);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  pfasst::init(argc, argv,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_opts,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_logs);
  int result = 1, max_result;  // GTest return value 1 (failure), 0 (success)
  // cppcheck-suppress redundantAssignment
  result = RUN_ALL_TESTS();
  MPI_Allreduce(&result, &max_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return max_result;
}
