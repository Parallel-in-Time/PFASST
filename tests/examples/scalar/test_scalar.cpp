/*
 * Tests for the scalar example solving the test equation
 */ 

#include <cmath>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#define PFASST_UNIT_TESTING
#include "../examples/scalar/scalar_sdc.cpp"
#undef PFASST_UNIT_TESTING

MATCHER (DoubleMore, "")
{
  return get<0>(arg) > get<1>(arg);
}

/*
 * For Lobatto nodes, the resulting method should of order 2*M-2 with M=number of nodes
 * The test below verifies that the code approximately (up to a safety factor) reproduces
 * the theoretically expected rate of convergence 
 */
TEST(ConvergenceTest, ScalarSDC)
{
    
  const complex<double> lambda = complex<double>(-1.0, 1.0);
  const double Tend = 4.0;
  const vector<int> nsteps = { 2, 5, 10, 15, 20 };
  int nsteps_l = nsteps.size();

  vector<double> err(nsteps.size());
  vector<double> convrate(nsteps.size()-1);
  vector<double> expected_cr(nsteps.size()-1); 
    
  double dt;
    
  // Test converge up to 6 nodes: For more nodes, errors are of the order of machine 
  // precision and convergence can no longer be monitored  
  for ( int nnodes = 2; nnodes<=6; ++nnodes)
  {
    // Expect convergence rate of 2*nodes-2 from collocation formula,
    // doing an identical number of iteration should suffice to reach this as each
    // iteration should increase order by one
    int niters = 2*nnodes-2;
  
    for ( int i = 0; i<=nsteps_l-1; ++i)
    {
      dt = Tend/double(nsteps[i]);
      err[i] = run_scalar_sdc(nsteps[i], dt, nnodes, niters, lambda);
    }
  
    for (int i = 0; i<=nsteps_l-2; ++i)
    {
      convrate[i] = log10(err[i+1]/err[i])/log10(double(nsteps[i])/double(nsteps[i+1]));
      // The expected convergence rate for Lobatto nodes is 2*nnodes-2, but because
      // it will typically be not matched exactly, put in a security factor
      expected_cr[i] = 0.9*double(2*nnodes-2);
    }
  
    // NOTE: There is probably a much more elegant way to test this, because
    // expected_cr contains the same value in all entries. But I could not so far figure
    // out how to build a more clever MATCHER here so far....
    EXPECT_THAT(convrate, testing::Pointwise(DoubleMore(), expected_cr ));
  
  }
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}