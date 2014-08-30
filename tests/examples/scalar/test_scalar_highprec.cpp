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
    
    const complex<double> lambda = complex<double>(-1.0,1.0);
    const double dt = 0.2; // = Tend for single step
    const size_t nsteps = 1;
    const size_t niters = 30;
    const size_t nnodes = 8; 
    size_t nnodes_in_call;
    double err;
    
    void set_parameters()
    {
      switch (this->nodetype)
      {
        case pfasst::QuadratureType::GaussLobatto:
          this->nnodes_in_call = this->nnodes;
          break;
        
        case pfasst::QuadratureType::GaussLegendre:
          this->nnodes_in_call = this->nnodes + 2;
          break;
          
        case pfasst::QuadratureType::GaussRadau:
          this->nnodes_in_call = this->nnodes + 1;
          break;
        
        case pfasst::QuadratureType::ClenshawCurtis:
          this->nnodes_in_call = this->nnodes + 1;
          break;
        
        case pfasst::QuadratureType::Uniform:
          this->nnodes_in_call = this->nnodes + 1;
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
      this->err = run_scalar_sdc(this->nsteps, this->dt, this->nnodes_in_call,
                                 this->niters, this->lambda, this->nodetype);
    }
    
    virtual void TearDown()
    {}
};

TEST_P(HighPrecisionTest, AllNodes)
{
  EXPECT_THAT(err, Le<double>(5e-12)) << "Failed to bring relative error below 5e-12";
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