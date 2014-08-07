#include <memory>

#include <pfasst.hpp>
#include <pfasst/sdc.hpp>

#include "particle_2d.hpp"
#include "boris_sweeper.hpp"

error_map run_boris_sdc()
{
  pfasst::SDC<> sdc;

  const size_t nsteps     = 4;
  const double dt         = 0.01;
  const size_t nnodes     = 5;
  const size_t nparticles = 1;
  const size_t niters     = 4;

  const double mass = 1.0;
  const double charge = 1.0;

  auto nodes   = pfasst::compute_nodes(nnodes, "gauss-lobatto");
  auto factory = make_shared<Particle2DFactory<double, double>>(mass, charge);
  auto sweeper = make_shared<BorisSweeper<double, double>>();

  sweeper->set_nodes(nodes);
  sweeper->set_factory(factory);

  sdc.add_level(sweeper);
  sdc.set_duration(0.0, nsteps*dt, dt, niters);
  sdc.setup();

  auto q0 = sweeper->get_state(0);
  sweeper->exact(q0, 0.0);

  sdc.run();

  return sweeper->get_errors();
}


int main(int /* argn */, char** /* argc */)
{
  run_boris_sdc();
}
