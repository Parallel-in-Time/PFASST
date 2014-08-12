/*
 * Tests for the scalar example solving the test equation
 */
#include <cmath>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace ::testing;

#define PFASST_UNIT_TESTING
#include "../examples/scalar/scalar_sdc.cpp"
#undef PFASST_UNIT_TESTING

/*
 * parameterized test fixture with number of nodes as parameter
 */
class ConvergenceTest
  : public ::testing::TestWithParam<size_t>
{
  protected:
    size_t nnodes; // parameter

    complex<double> lambda = complex<double>(-1.0, 1.0);
    double Tend = 4.0;
    vector<size_t> nsteps = { 2, 5, 10, 15, 20 };
    size_t nsteps_l;
    vector<double> err;
    vector<double> convrate_lobatto;
    vector<double> convrate_legendre;
    double dt;
    size_t niters;
    pfasst::QuadratureType nodetype = pfasst::QuadratureType::GaussLobatto;

  public:
    virtual void SetUp()
    {
      this->nnodes = this->GetParam();
      this->nsteps_l = this->nsteps.size();
      this->err.resize(this->nsteps.size());
      this->convrate_lobatto.resize(this->nsteps.size());
    
      // Run first for Lobatto nodes which is the default nodetype

      // Expect convergence rate of 2*nodes-2 from collocation formula,
      // doing an identical number of iteration should suffice to reach this as each
      // iteration should increase order by one
      this->niters = 2 * nnodes - 2;

      // run to compute errors
      for (size_t i = 0; i <= nsteps_l - 1; ++i) {
        dt = Tend / double(nsteps[i]);
        err[i] = run_scalar_sdc(nsteps[i], dt, nnodes, niters, lambda, nodetype);
      }

      // compute convergence rates
      for (size_t i = 0; i <= nsteps_l - 2; ++i) {
        convrate_lobatto[i] = log10(err[i+1] / err[i]) / log10(double(nsteps[i]) / double(nsteps[i + 1]));
      }
       
      /*    
       * Compute convergence rates for Legendre nodes. 
       * Note that this is not a particularly nice way to do the test, because both types 
       * of nodes are tested in the same test. Having separate tests would be nicer...
       */ 
      nodetype = pfasst::QuadratureType::GaussLegendre;
      this->convrate_legendre.resize(this->nsteps.size());
      
      // refine parameter
      complex<double> lambda = complex<double>(-1.0, 2.0);      
      this->Tend = 6.0;
      this->nsteps = { 2, 4, 6, 8, 10 };
      this->niters = 2*nnodes;
      
      // run to compute errors
      for (size_t i = 0; i <= nsteps_l - 1; ++i) {
        dt = Tend / double(nsteps[i]);
        
      // NOTE: Setting M for Legendre nodes actually only results in M-2 "real" nodes,
      // because the first and the last are used for initial and final value.      
      // Hence us nnodes+2 as argument  
        err[i] = run_scalar_sdc(nsteps[i], dt, nnodes+2, niters, lambda, nodetype);
      }
      
      for (size_t i = 0; i <= nsteps_l - 2; ++i) {
        convrate_legendre[i] = log10(err[i+1] / err[i]) / log10(double(nsteps[i]) / double(nsteps[i + 1]));
      }
            
    }

    virtual void TearDown()
    {}
};

/*
 * The test below verifies that the code approximately (up to a safety factor) reproduces
 * the theoretically expected rate of convergence for Lobatto and Legendre nodes
 */
TEST_P(ConvergenceTest, GaussNodes)
{
  for (size_t i = 0; i <= nsteps_l - 2; ++i) {

    // Lobatto nodes reproduce the convergence rate quite accurately, so use DoubleNear
    EXPECT_THAT(convrate_lobatto[i],
                testing::DoubleNear(double(2 * nnodes - 2), 0.99)) << "Convergence rate at node "
                                                                   << i
                                                                   << " not within expected range";
                                                                   
    // convergence rates for Legendre nodes should be 2*nodes but are actually better, so
    // use Ge here
    EXPECT_THAT(convrate_legendre[i],
                testing::Ge(double(2 * nnodes))) << "Convergence rate at node "
                                                 << i
                                                 << " not within expected range";                                                                

  }
}

INSTANTIATE_TEST_CASE_P(ScalarSDC, ConvergenceTest, Range<size_t>(2, 7));

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
