/*
 * Tests for the scalar example solving the test equation
 */
#include <cmath>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include <pfasst/quadrature.hpp>

#define PFASST_UNIT_TESTING
#include "../examples/scalar/scalar_sdc.cpp"
#undef PFASST_UNIT_TESTING
using namespace pfasst::examples::scalar;

/*
 * parameterized test fixture with number of nodes as parameter
 */
class ConvergenceTest
  : public TestWithParam<tuple<size_t, pfasst::quadrature::QuadratureType>>
{
  protected:
    size_t nnodes;
    size_t niters;
    double end_time;
    vector<size_t> nsteps;
    vector<double> err;
    vector<double> convrate;
    complex<double> lambda;
    pfasst::quadrature::QuadratureType nodetype;

  public:
    virtual void SetUp()
    {
      this->nnodes = get<0>(GetParam());
      this->nodetype = get<1>(GetParam());

      switch (this->nodetype)
      {
        case pfasst::quadrature::QuadratureType::GaussLobatto:
          this->niters = 2 * this->nnodes - 2;
          this->end_time = 4.0;
          this->lambda = complex<double>(-1.0, 1.0);
          this->nsteps = { 2, 5, 10, 15 };
          break;

        case pfasst::quadrature::QuadratureType::GaussLegendre:
          this->niters = 2 * this->nnodes;
          this->end_time = 6.0;
          this->lambda = complex<double>(-1.0, 2.0);
          this->nsteps = { 2, 4, 6, 8, 10 };
          break;

        case pfasst::quadrature::QuadratureType::GaussRadau:
          this->niters = 2 * this->nnodes;
          this->end_time = 5.0;
          this->lambda = complex<double>(-1.0, 2.0);
          this->nsteps = { 4, 6, 8, 10, 12 };
          break;

        case pfasst::quadrature::QuadratureType::ClenshawCurtis:
          this->niters = this->nnodes + 1;
          this->end_time = 1.0;
          this->lambda = complex<double>(-1.0, 1.0);
          this->nsteps = { 7, 9, 11, 13 };
          break;

        case pfasst::quadrature::QuadratureType::Uniform:
          this->niters = this->nnodes;
          this->end_time = 5.0;
          this->lambda = complex<double>(-1.0, 1.0);
          this->nsteps = { 9, 11, 13, 15 };
          break;

        default:
          break;
      }

      // run to compute errors
      for (size_t i = 0; i < this->nsteps.size(); ++i) {
        auto dt = this->end_time / double(this->nsteps[i]);
        this->err.push_back(run_scalar_sdc(this->nsteps[i], dt, this->nnodes,
                                           this->niters, this->lambda, this->nodetype));
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
TEST_P(ConvergenceTest, AllNodes)
{
  int order;
  string quad;
  double fudge = 0.9;

  switch (this->nodetype)
  {
  case pfasst::quadrature::QuadratureType::GaussLobatto:
    order = 2 * this->nnodes - 2;
    quad = "Gauss-Lobatto";
    break;

  case pfasst::quadrature::QuadratureType::GaussLegendre:
    order = 2 * this->nnodes;
    quad = "Gauss-Legendre";
    break;

  case pfasst::quadrature::QuadratureType::GaussRadau:
    order = 2 * this->nnodes;
    quad = "Gauss-Radau";
    break;

  case pfasst::quadrature::QuadratureType::ClenshawCurtis:
    order = this->nnodes;
    quad = "Clenshaw-Curtis";
    break;

  case pfasst::quadrature::QuadratureType::Uniform:
    order = this->nnodes;
    fudge = 0.8;
    quad = "Uniform";
    break;

  default:
    EXPECT_TRUE(false);
    break;
  }

  for (size_t i = 0; i < this->nsteps.size() - 1; ++i) {
    EXPECT_THAT(convrate[i], Ge<double>(fudge * order)) << "Convergence rate for "
                                                        << this->nnodes << " " << quad << " nodes"
                                                        << " for nsteps " << this->nsteps[i]
                                                        << " not within expected range.";
  }
}

INSTANTIATE_TEST_CASE_P(ScalarSDC, ConvergenceTest,
                        Combine(Range<size_t>(3, 7),
                                Values(pfasst::quadrature::QuadratureType::GaussLobatto,
                                       pfasst::quadrature::QuadratureType::GaussLegendre,
                                       pfasst::quadrature::QuadratureType::GaussRadau,
                                       pfasst::quadrature::QuadratureType::ClenshawCurtis,
                                       pfasst::quadrature::QuadratureType::Uniform))
);

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
