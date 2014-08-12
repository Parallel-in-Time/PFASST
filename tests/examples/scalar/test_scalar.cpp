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

    const complex<double> lambda = complex<double>(-1.0, 1.0);
    const double Tend = 4.0;
    const vector<size_t> nsteps = { 2, 5, 10, 15, 20 };
    size_t nsteps_l;
    vector<double> err;
    vector<double> convrate;
    double dt;
    size_t niters;
    const pfasst::QuadratureType nodetype = pfasst::QuadratureType::GaussLobatto;

  public:
    virtual void SetUp()
    {
      this->nnodes = this->GetParam();
      this->nsteps_l = this->nsteps.size();
      this->err.resize(this->nsteps.size());
      this->convrate.resize(this->nsteps.size());

      // Expect convergence rate of 2*nodes-2 from collocation formula,
      // doing an identical number of iteration should suffice to reach this as each
      // iteration should increase order by one
      this->niters = 2 * nnodes - 2;

      for (size_t i = 0; i <= nsteps_l - 1; ++i) {
        dt = Tend / double(nsteps[i]);
        err[i] = run_scalar_sdc(nsteps[i], dt, nnodes, niters, lambda, nodetype);
      }
    }

    virtual void TearDown()
    {}
};

/*
 * For Lobatto nodes, the resulting method should of order 2*M-2 with M=number of nodes
 * The test below verifies that the code approximately (up to a safety factor) reproduces
 * the theoretically expected rate of convergence 
 */
TEST_P(ConvergenceTest, GaussLobattoNodes)
{
  for (size_t i = 0; i <= nsteps_l - 2; ++i) {
    convrate[i] = log10(err[i+1] / err[i]) / log10(double(nsteps[i]) / double(nsteps[i + 1]));

    EXPECT_THAT(convrate[i],
                testing::DoubleNear(double(2 * nnodes - 2), 0.99)) << "Convergence rate at node "
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
