/*
 * Tests for the scalar example solving the test equation
 */ 

#include <cmath>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#define PFASST_UNIT_TESTING
#include "../examples/scalar/scalar_sdc.cpp"
#undef PFASST_UNIT_TESTING


/*
 * For Lobatto nodes, the resulting method should of order 2*M-2 with M=number of nodes
 * The test below verifies that the code approximately (up to a safety factor) reproduces
 * the theoretically expected rate of convergence 
 */
TEST(ConvergenceTest, ScalarSDC)
{
  const complex<double> lambda = complex<double>(-1.0, 1.0);
  const double Tend = 4.0;
  const vector<size_t> nsteps = { 2, 5, 10, 15, 20 };
  size_t nsteps_l = nsteps.size();

  vector<double> err(nsteps.size());
  vector<double> convrate(nsteps.size()-1);

  double dt;

  // Test converge up to 6 nodes: For more nodes, errors are of the order of machine
  // precision and convergence can no longer be monitored
  for ( size_t nnodes = 2; nnodes<=6; ++nnodes)
  {
    // Expect convergence rate of 2*nodes-2 from collocation formula,
    // doing an identical number of iteration should suffice to reach this as each
    // iteration should increase order by one
    size_t niters = 2*nnodes-2;

    for ( size_t i = 0; i<=nsteps_l-1; ++i)
    {
      dt = Tend/double(nsteps[i]);
      err[i] = run_scalar_sdc(nsteps[i], dt, nnodes, niters, lambda);
    }

    for ( size_t i = 0; i<=nsteps_l-2; ++i)
    {
      convrate[i] = log10(err[i+1]/err[i])/log10(double(nsteps[i])/double(nsteps[i+1]));

      EXPECT_THAT(convrate[i],
                  testing::DoubleNear(double(2*nnodes-2), 0.99)) << "Convergence rate at node "
                                                                 << i
                                                                 << " not within expected range";
    }
  }
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
