#include <gtest/gtest.h>

#include <pfasst/sdc.hpp>
#include "examples/advection_diffusion/fft.hpp"
#include "examples/advection_diffusion/advection_diffusion_sweeper.hpp"

using pfasst::sdc_factory;
using namespace std;


TEST(SdcTest, ConvergeAdvectionDiffusion)
{
  typedef double time_precision;
  typedef AdvectionDiffusionSweeper<time_precision> sweeperT;
  typedef pfasst::encap::VectorFactory<time_precision> dataT;

  const size_t         nsteps        = 4;
  const time_precision dt            = 0.01;
  const size_t         nnodes        = 5;
  const size_t         ndofs         = 64;
  const size_t         niters        = 4;
  const time_precision initial_value = 0.0;

  auto sdc = sdc_factory<sweeperT,
                         dataT,
                         time_precision>(niters, nsteps, dt, ndofs, nnodes, "gauss-lobatto", 
                                         initial_value);

  sdc.run();

  fftw_cleanup();

  auto level = sdc.get_finest<sweeperT>();

  EXPECT_EQ(level->get_f1evals(), 65);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
