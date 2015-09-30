/*
 * Convergence tests for the advection-diffusion examples.
 */
#include <cmath>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include <pfasst/quadrature.hpp>

#define PFASST_UNIT_TESTING
#include "../examples/advection_diffusion/mpi_pfasst.cpp"
#undef PFASST_UNIT_TESTING
using namespace pfasst::examples::advection_diffusion;

/*
 * parameterized test fixture with number of nodes as parameter
 */
class ConvergenceTest
  : public TestWithParam<tuple<size_t, pfasst::quadrature::QuadratureType>>
{
  protected:
    size_t nnodes;
    size_t niters;
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
        case pfasst::quadrature::QuadratureType::GaussLobatto:
          this->niters = 2 * this->nnodes - 2;
          this->nsteps = { 4, 8, 16, 32 };
          break;

        case pfasst::quadrature::QuadratureType::GaussLegendre:
          this->niters = 2 * this->nnodes;
          this->nsteps = { 4, 8, 16, 32 };
          break;

        case pfasst::quadrature::QuadratureType::ClenshawCurtis:
          this->niters = this->nnodes;
          this->nsteps = { 4, 8, 16, 32 };
          break;

        case pfasst::quadrature::QuadratureType::Uniform:
          this->niters = this->nnodes;
          this->nsteps = { 4, 8, 16, 32 };
          break;

        default:
          break;
      }

      // run to compute errors
      for (size_t i = 0; i < this->nsteps.size(); ++i) {
        auto dt = 0.5 / double(this->nsteps[i]);

        auto errors = run_mpi_pfasst(
          0.0,                  // abs_res_tol
          0.0,                  // rel_res_tol
          this->niters,         // niters
          this->nsteps[i],      // nsteps
          dt,                   // dt
          128,                  // ndofs_f
          64,                   // ndofs_c
          this->nnodes,         // nnodes_f
          (this->nnodes+1)/2-1  // nnodes_c
                                     );

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        auto last_error = pfasst::examples::advection_diffusion::ktype(this->nsteps[i]-size+rank, this->niters-1);
        this->err.push_back(errors.at(last_error));
      }

      // compute convergence rates
      for (size_t i = 0; i < this->nsteps.size() - 1; ++i) {
        this->convrate.push_back(log10(this->err[i+1] / this->err[i]) /
                                 log10(double(this->nsteps[i]) / double(this->nsteps[i + 1])));
      }
    }
};

/*
 * The test below verifies that the code reproduces the theoretically expected rate of convergence
 * (or better) ON THE LAST PROCESSOR.
 */
TEST_P(ConvergenceTest, AllNodes)
{
  int order;
  int rank, size;
  string quad;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank != size-1) {
    EXPECT_TRUE(true);
    return;
  }

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
    order = 2 * this->nnodes - 1;
    quad = "Gauss-Radau";
    break;

  case pfasst::quadrature::QuadratureType::ClenshawCurtis:
    order = this->nnodes;
    quad = "Clenshaw-Curtis";
    break;

  case pfasst::quadrature::QuadratureType::Uniform:
    order = this->nnodes;
    quad = "Uniform";
    break;

  default:
    EXPECT_TRUE(false);
    break;
  }

  for (size_t i = 0; i < this->nsteps.size() - 1; ++i) {
    EXPECT_THAT(convrate[i], Ge<double>(order)) << "Convergence rate for "
                                                << this->nnodes << " " << quad << " nodes"
                                                << " for nsteps " << this->nsteps[i]
                                                << " not within expected range.";
  }
}

INSTANTIATE_TEST_CASE_P(AdvectionDiffusionPFASST, ConvergenceTest,
                        Combine(Range<size_t>(5, 6),
                                Values(pfasst::quadrature::QuadratureType::GaussLobatto,
                                       //                                       pfasst::quadrature::QuadratureType::GaussLegendre,
                                       //                                       pfasst::quadrature::QuadratureType::GaussRadau,
                                       pfasst::quadrature::QuadratureType::ClenshawCurtis,
                                       pfasst::quadrature::QuadratureType::Uniform))
);

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  pfasst::init(argc, argv,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_opts,
               pfasst::examples::advection_diffusion::AdvectionDiffusionSweeper<>::init_logs);
  int result = 1, max_result;  // GTest return value 1 (failure), 0 (success)
  result = RUN_ALL_TESTS();
  MPI_Allreduce(&result, &max_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Finalize();
  return max_result;
}
