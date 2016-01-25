#include <memory>

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/controller/sdc.hpp>

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
                                      const size_t nparticles, const size_t niters,
                                      const double abs_res_tol, const double rel_res_tol)
      {
        SDC<> sdc;

        const double mass = 1.0;
        const double charge = 1.0;

        auto quad    = quadrature::quadrature_factory<double>(nnodes, quadrature::QuadratureType::GaussLobatto);
        auto factory = make_shared<ParticleCloudFactory<double>>(nparticles, 3, mass, charge);

        shared_ptr<bindings::WrapperInterface<double, double>> impl_solver = \
          make_shared<bindings::WrapperSimplePhysicsSolver<double, double>>();
        bindings::setup(dynamic_pointer_cast<bindings::WrapperSimplePhysicsSolver<double, double>>(impl_solver));

        string data_file = "s" + to_string(nsteps) + "_i" + to_string(niters) + "_dt" + to_string(dt) + "_m" + to_string(nnodes) + "_p" + to_string(nparticles) + ".csv";
        auto sweeper = make_shared<BorisSweeper<double, double>>(impl_solver, data_file);

        sweeper->set_quadrature(quad);
        sweeper->set_factory(factory);
        sweeper->set_residual_tolerances(abs_res_tol, rel_res_tol);

        sdc.add_level(sweeper);
        sdc.set_duration(0.0, nsteps*dt, dt, niters);
        sdc.setup();

        shared_ptr<Particle<double>> center = make_shared<Particle<double>>();
        center->pos()[0] = 10;
        center->vel()[0] = 100;
        center->vel()[2] = 100;

        shared_ptr<ParticleCloud<double>> q0 = dynamic_pointer_cast<ParticleCloud<double>>(sweeper->get_start_state());
        q0->distribute_around_center(center);
        ML_CLOG(INFO, "Boris", OUT::green << "Initial Particle: "
                               << *(dynamic_pointer_cast<ParticleCloud<double>>(sweeper->get_start_state())));

        sweeper->set_initial_energy();
        sdc.run();

        return sweeper->get_errors();
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#ifndef PFASST_UNIT_TESTING
int main(int argc, char** argv)
{
  pfasst::init(argc, argv,
               pfasst::examples::boris::init_opts<>,
               pfasst::examples::boris::init_logs<>);

  const size_t nsteps      = pfasst::config::get_value<size_t>("num_steps", 1);
  const double dt          = pfasst::config::get_value<double>("delta_step", 0.015625);
  const size_t nnodes      = pfasst::config::get_value<size_t>("num_nodes", 5);
  const size_t nparticles  = pfasst::config::get_value<size_t>("num_particles", 1);
  const size_t niters      = pfasst::config::get_value<size_t>("num_iter", 2);
  const double abs_res_tol = pfasst::config::get_value<double>("abs_res_tol", 0.0);
  const double rel_res_tol = pfasst::config::get_value<double>("rel_res_tol", 0.0);

  ML_CLOG(INFO, "Boris", "nsteps=" << nsteps << ", "
                         << "dt=" << dt << ", "
                         << "nnodes=" << nnodes << ", "
                         << "nparticles=" << nparticles << ", "
                         << "niter=" << niters << ", "
                         << "abs res=" << abs_res_tol << ", "
                         << "rel res=" << rel_res_tol);

  pfasst::examples::boris::run_boris_sdc<double>(nsteps, dt, nnodes, nparticles, niters, abs_res_tol, rel_res_tol);
}
#endif
