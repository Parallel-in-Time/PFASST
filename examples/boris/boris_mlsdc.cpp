#include <memory>

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/sdc.hpp>
#include <pfasst/controller/mlsdc.hpp>

#include "pfasst/contrib/encap/particle.hpp"
#include "pfasst/contrib/encap/particle_cloud.hpp"
#include "pfasst/contrib/util/bindings/wrapper_interface.hpp"
#include "bindings/wrapper_simple_physics_solver.hpp"
#include "pfasst/contrib/sweeper/boris_sweeper.hpp"
#include "pfasst/contrib/util/injective_transfer.hpp"


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
        MLSDC<> controller;

        const double mass = 1.0;
        const double charge = 1.0;

        shared_ptr<bindings::WrapperInterface<double, double>> impl_solver = \
          make_shared<bindings::WrapperSimplePhysicsSolver<double, double>>();
        bindings::setup(dynamic_pointer_cast<bindings::WrapperSimplePhysicsSolver<double, double>>(impl_solver));

        auto quad1    = quadrature::quadrature_factory<double>(nnodes, quadrature::QuadratureType::GaussLobatto);
        auto factory1 = make_shared<ParticleCloudFactory<double>>(nparticles, 3, mass, charge);
        string data_file1 = "s" + to_string(nsteps) + "_i" + to_string(niters) + "_dt" + to_string(dt) + "_m" + to_string(nnodes) + "_p" + to_string(nparticles) + "_level1.csv";
        auto sweeper1 = make_shared<BorisSweeper<double, double>>(impl_solver, data_file1);
        auto transfer1 = make_shared<InjectiveTransfer<double, double>>();
        sweeper1->set_quadrature(quad1);
        sweeper1->set_factory(factory1);
        controller.add_level(sweeper1, transfer1);

        auto quad2    = quadrature::quadrature_factory<double>(nnodes, quadrature::QuadratureType::GaussLobatto);
        auto factory2 = make_shared<ParticleCloudFactory<double>>(nparticles, 3, mass, charge);
        string data_file2 = "s" + to_string(nsteps) + "_i" + to_string(niters) + "_dt" + to_string(dt) + "_m" + to_string(nnodes) + "_p" + to_string(nparticles) + "_level2.csv";
        auto sweeper2 = make_shared<BorisSweeper<double, double>>(impl_solver, data_file2);
        auto transfer2 = make_shared<InjectiveTransfer<double, double>>();
        sweeper2->set_quadrature(quad2);
        sweeper2->set_factory(factory2);
        controller.add_level(sweeper2, transfer2);

        controller.set_duration(0.0, nsteps*dt, dt, niters);
        controller.setup();

        shared_ptr<Particle<double>> center = make_shared<Particle<double>>();
        center->pos()[0] = 10;
        center->vel()[0] = 100;
        center->vel()[2] = 100;

        auto fine_sweeper = controller.get_finest<BorisSweeper<double, double>>();
        shared_ptr<ParticleCloud<double>> q0 = dynamic_pointer_cast<ParticleCloud<double>>(fine_sweeper->get_start_state());
        q0->distribute_around_center(center);
        CLOG(INFO, "Boris") << "Initial Particle (fine) : " << *(dynamic_pointer_cast<ParticleCloud<double>>(fine_sweeper->get_start_state()));
        fine_sweeper->set_initial_energy();

//         auto coarse_sweeper = controller.get_coarsest<BorisSweeper<double, double>>();
//         shared_ptr<ParticleCloud<double>> q0_coarse = dynamic_pointer_cast<ParticleCloud<double>>(coarse_sweeper->get_start_state());
//         q0_coarse->copy(q0);
//         LOG(INFO) << OUT::green << "Initial Particle (coarse): " << *(dynamic_pointer_cast<ParticleCloud<double>>(coarse_sweeper->get_start_state()));
//         coarse_sweeper->set_initial_energy();

        controller.run();

        return fine_sweeper->get_errors();
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::init(argc, argv, pfasst::examples::boris::init_opts<double>);
  pfasst::init(argc, argv);
  pfasst::log::add_custom_logger("Boris");
  pfasst::log::add_custom_logger("BorisTransfer");
  pfasst::log::add_custom_logger("Solver");

  const size_t nsteps     = pfasst::config::get_value<size_t>("num_steps", 1);
  const double dt         = pfasst::config::get_value<double>("delta_step", 0.015625);
  const size_t nnodes     = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t nparticles = pfasst::config::get_value<size_t>("num_particles", 1);
  const size_t niters     = pfasst::config::get_value<size_t>("num_iter", 2);

  pfasst::examples::boris::run_boris_sdc<double>(nsteps, dt, nnodes, nparticles, niters);
}
#endif
