/*
 * Tests for the scalar example solving the test equation
 */
#include <cmath>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace ::testing;

#include <pfasst/quadrature.hpp>

#define PFASST_UNIT_TESTING
#include "../examples/scalar/scalar_sdc.cpp"
#undef PFASST_UNIT_TESTING

/*
 * parameterized test that verifies that given sufficiently many nodes and iterations,
 * SDC can reproduce the analytic solution with very high precision
 */
class HighPrecisionTest
  : public TestWithParam<pfasst::QuadratureType>
{
  protected:
    pfasst::QuadratureType nodetype; // parameter
    
    complex<double> lambda;
    double Tend = 0.5;
    size_t nsteps = 1;
    size_t niters = 30;
    size_t nnodes = 8; 
    size_t nnodes_in_call;
    
    void set_parameters()
    {
      switch (this->nodetype)
      {
        case pfasst::QuadratureType::GaussLobatto:
          break;
          
        default:
          break;
      }
      
    }
            
  public:
    virtual void SetUp()
    {
      this->nodetype = GetParam();
      this->set_parameters();
    }
    
    virtual void TearDown()
    {}
};

TEST_P(HighPrecisionTest, AllNodes)
{
}

INSTANTIATE_TEST_CASE_P(ScalarSDC, HighPrecisionTest,
                                Values(pfasst::QuadratureType::GaussLobatto,
                                       pfasst::QuadratureType::GaussLegendre,
                                       pfasst::QuadratureType::GaussRadau,
                                       pfasst::QuadratureType::ClenshawCurtis,
                                       pfasst::QuadratureType::Uniform)
);

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}