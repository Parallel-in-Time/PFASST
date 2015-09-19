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

    vector<size_t> nsteps;
    size_t nsteps_l;
    vector<double> err;
    vector<double> convrate;
    double dt;
    size_t niters;
    pfasst::quadrature::QuadratureType nodetype;

    // set parameters base on node type
    void set_parameters()
    {
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
    }

  public:
    virtual void SetUp()
    {
      this->nnodes = get<0>(GetParam());
      this->nodetype = get<1>(GetParam());
      this->set_parameters();
      this->err.resize(this->nsteps.size());
      this->convrate.resize(this->nsteps.size()-1);

      // run to compute errors
      for (size_t i = 0; i < this->err.size(); ++i) {
        this->dt = 0.5 / double(this->nsteps[i]);

        auto errors = run_mpi_pfasst(
          0.0,                  // abs_res_tol
          0.0,                  // rel_res_tol
          this->niters,         // niters
          this->nsteps[i],      // nsteps
          this->dt,             // dt
          128,                  // ndofs_f
          64,                   // ndofs_c
          this->nnodes,         // nnodes_f
          (this->nnodes+1)/2-1  // nnodes_c
                                     );

        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        auto last_error = pfasst::examples::advection_diffusion::ktype(this->nsteps[i]-size+rank, this->niters-1);
        this->err[i] = errors.at(last_error);
      }

      // compute convergence rates
      for (size_t i = 0; i < this->convrate.size(); ++i) {
        this->convrate[i] = log10(this->err[i+1] / this->err[i]) /
                              log10(double(this->nsteps[i]) / double(this->nsteps[i + 1]));
      }
    }

    virtual void TearDown()
    {}
};

/*
 * The test below verifies that the code reproduces the theoretically expected rate of convergence
 * (or better) ON THE LAST PROCESSOR.
 */
TEST_P(ConvergenceTest, AllNodes)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank != size-1) {
    EXPECT_TRUE(true);
    return;
  }

  for (size_t i = 0; i < this->convrate.size(); ++i) {
    switch (this->nodetype)
    {
      case pfasst::quadrature::QuadratureType::GaussLobatto:
        // Expect convergence rate of 2*nodes-2 from collocation formula, doing an identical number
        // of iteration should suffice to reach this as each iteration should increase order by one
        EXPECT_THAT(convrate[i], Ge<double>(2 * this->nnodes - 2)) << "Convergence rate for "
                                                                   << this->nnodes
                                                                   << " Gauss-Lobatto nodes "
                                                                   << " at nsteps no. " << i
                                                                   << " not within expected range.";
        break;

      case pfasst::quadrature::QuadratureType::GaussLegendre:
        // convergence rates for Legendre nodes should be 2*nodes but are actually better, so
        // use Ge here
        EXPECT_THAT(convrate[i], Ge<double>(2 * this->nnodes)) << "Convergence rate for "
                                                               << this->nnodes
                                                               << " Gauss-Legendre nodes "
                                                               << " at nsteps no. " << i
                                                               << " not within expected range.";
        break;

      case pfasst::quadrature::QuadratureType::GaussRadau:
        EXPECT_THAT(convrate[i], Ge<double>(2 * this->nnodes - 1)) << "Convergence rate for "
                                                                         << this->nnodes
                                                                         << " Gauss-Radu nodes "
                                                                         << " at nsteps no. " << i
                                                                         << " not within expected range.";
        break;

      case pfasst::quadrature::QuadratureType::ClenshawCurtis:
        EXPECT_THAT(convrate[i], Ge<double>(this->nnodes))  << "Convergence rate for "
                                                            << this->nnodes
                                                            << " Clenshaw-Curtis nodes "
                                                            << " at nsteps no. " << i
                                                            << " not within expected range.";
        break;

      case pfasst::quadrature::QuadratureType::Uniform:
        EXPECT_THAT(convrate[i], Ge<double>(this->nnodes)) << "Convergence rate for "
                                                           << this->nnodes
                                                           << " equidistant nodes "
                                                           << " at nsteps no. " << i
                                                           << " not within expected range.";
        break;

      default:
        EXPECT_TRUE(false);
        break;
    }
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
