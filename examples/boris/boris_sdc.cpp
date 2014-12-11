#include <memory>

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/sdc.hpp>

#include "particle.hpp"
#include "particle_cloud.hpp"
#include "bindings/wrapper_interface.hpp"
#include "bindings/wrapper_simple_physics_solver.hpp"
#include "boris_sweeper.hpp"

namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      template<typename scalar>
      error_map<scalar> run_boris_sdc(const size_t nsteps, const scalar dt, const size_t nnodes,
                                      const size_t nparticles, const size_t niters)
      {
        SDC<> sdc;

        const double mass = 1.0;
        const double charge = 1.0;

        const double epsilon = config::get_value<scalar>("epsilon", -1.0);

        auto quad    = quadrature::quadrature_factory<double>(nnodes, quadrature::QuadratureType::GaussLobatto);
        auto factory = make_shared<ParticleCloudFactory<double>>(nparticles, 3, mass, charge);

        shared_ptr<bindings::WrapperInterface<double, double>> impl_solver = \
          make_shared<bindings::WrapperSimplePhysicsSolver<double, double>>();
        bindings::setup(impl_solver);
        auto sweeper = make_shared<BorisSweeper<double, double>>(impl_solver);

        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);

        sdc.add_level(sweeper);
        sdc.set_duration(0.0, nsteps*dt, dt, niters);
        sdc.setup();

        auto p0 = dynamic_pointer_cast<ParticleCloud<double>>(sweeper->get_state(0));
        cout << p0->size() << endl;
  //       p0[0]->pos()[0] = 10;
  //       p0[0]->vel()[0] = 100;
  //       p0[0]->vel()[2] = 100;

        sweeper->set_initial_energy();
  //       LOG(INFO) << OUT::green << "Initial Particle: " << *(dynamic_pointer_cast<ParticleCloud<double>>(sweeper->get_state(0)));
        sdc.run();

        return sweeper->get_errors();
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::examples::boris::enable_config_options<double>();
  pfasst::init(argc, argv);

  const size_t nsteps     = pfasst::config::get_value<size_t>("num_steps", 10);
  const double dt         = pfasst::config::get_value<double>("delta_step", 0.015625);
  const size_t nnodes     = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t nparticles = 1;
  const size_t niters     = pfasst::config::get_value<size_t>("num_iter", 7);

  pfasst::examples::boris::run_boris_sdc<double>(nsteps, dt, nnodes, nparticles, niters);
}
#endif
