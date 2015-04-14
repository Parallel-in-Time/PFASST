/**
 * Tests for the van der Pol oscillator for the case nu=0, where it becomes
 * the linear oscillator and an analytical solution is available
 */
#include <cmath>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include <pfasst/quadrature.hpp>

#define PFASST_UNIT_TESTING
#include "../examples/vanderpol/vdp_sdc.cpp"
#undef PFASST_UNIT_TESTING
using namespace pfasst::examples::vdp;

/*
 * parameterized test fixture with number of nodes as parameter
 */
class VdPConvergenceTest
  : public TestWithParam<tuple<size_t, pfasst::quadrature::QuadratureType>>
{
  protected:
    size_t nnodes;

    const double nu = 0.0;
    const double x0 = 1.0;
    const double y0 = 0.5;

    size_t niters;
    double end_time;
    vector<size_t> nsteps;
    vector<double> err;
    vector<double> convrate;
    pfasst::quadrature::QuadratureType nodetype;

  public:
    virtual void SetUp()
    {
      this->nnodes = get<0>(GetParam());
      this->nodetype = get<1>(GetParam());

      switch (this->nodetype)
      {
        case pfasst::quadrature::QuadratureType::GaussLegendre:
          this->niters = 2 * this->nnodes;
          this->end_time = 0.88;
          this->nsteps = { 7, 9, 11 };
          break;

        case pfasst::quadrature::QuadratureType::GaussRadau:
          this->niters = 2 * this->nnodes - 1;
          this->end_time = 0.88;
          this->nsteps = { 7, 9, 11, 13 };
          break;

        default:
          break;
      }

      // run to compute errors
      for (size_t i = 0; i < this->nsteps.size(); ++i) {
        auto dt = this->end_time / this->nsteps[i];
        this->err.push_back(run_vdp_sdc(this->nsteps[i], dt, this->nnodes,
                                        this->niters, this->nu, this->x0, this->y0, this->nodetype));
      }

      // compute convergence rates
      for (size_t i = 0; i < this->nsteps.size() - 1; ++i) {
        this->convrate.push_back(log10(this->err[i+1] / this->err[i]) /
                                 log10(double(this->nsteps[i]) / double(this->nsteps[i + 1])));
      }
    }

};

/*
 * The test below verifies that the code approximately (up to a safety factor) reproduces
 * the theoretically expected rate of convergence
 */
TEST_P(VdPConvergenceTest, AllNodes)
{
  int order;
  string quad;

  switch (this->nodetype)
  {
  case pfasst::quadrature::QuadratureType::GaussLegendre:
    order = 2 * this->nnodes;
    quad = "Gauss-Legendre";
    break;
  case pfasst::quadrature::QuadratureType::GaussRadau:
    order = 2 * this->nnodes - 1;
    quad = "Gauss-Radau";
    break;
  default:
    EXPECT_TRUE(false);
    break;
  }

  for (size_t i = 0; i < this->nsteps.size() - 1; ++i) {
    EXPECT_THAT(convrate[i], Ge<double>(0.99 * order)) << "Convergence rate for "
                                                       << this->nnodes << " " << quad << " nodes"
                                                       << " for nsteps " << this->nsteps[i]
                                                       << " not within expected range.";
  }
}

INSTANTIATE_TEST_CASE_P(VanDerPol, VdPConvergenceTest,
                        Combine(Range<size_t>(3, 4),
                                Values(pfasst::quadrature::QuadratureType::GaussLegendre,
                                       pfasst::quadrature::QuadratureType::GaussRadau)));

int main(int argc, char** argv)
{
  pfasst::init(argc, argv);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
