#include "wrapper_simple_physics_solver.hpp"

#include "simple_physics_solver.hpp"
namespace solver = simple_physics_solver;

#include <pfasst/logging.hpp>

#define UNUSED(expr) (void)(expr)

namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      namespace bindings
      {
        template<typename scalar, typename time>
        size_t
        WrapperSimplePhysicsSolver<scalar, time>::vector_to_array(const vector<scalar>& vec, scalar* arr)
        {
          for (size_t i = 0; i < vec.size(); ++i) {
            arr[i] = vec[i];
          }
          return vec.size();
        }

        template<typename scalar, typename time>
        size_t
        WrapperSimplePhysicsSolver<scalar, time>::vector2d_to_array(const vector<scalar>& vec, scalar* arr)
        {
          const size_t npart = vec.size() / DIM;
          for (size_t p = 0; p < npart; ++p) {
            for (size_t d = 0; d < DIM; ++d) {
              arr[p * DIM + d] = vec[p * DIM + d];
            }
          }
          return vec.size() * DIM;
        }

        template<typename scalar, typename time>
        size_t
        WrapperSimplePhysicsSolver<scalar, time>::pack_positions(const particle_cloud_type& particles, scalar* packed)
        {
          return this->vector2d_to_array(particles->positions(), packed);
        }

        template<typename scalar, typename time>
        size_t
        WrapperSimplePhysicsSolver<scalar, time>::pack_velocities(const particle_cloud_type& particles, scalar* packed)
        {
          return this->vector2d_to_array(particles->velocities(), packed);
        }

        template<typename scalar, typename time>
        size_t
        WrapperSimplePhysicsSolver<scalar, time>::pack_charges(const particle_cloud_type& particles, scalar* packed)
        {
          return this->vector_to_array(particles->charges(), packed);
        }

        template<typename scalar, typename time>
        size_t
        WrapperSimplePhysicsSolver<scalar, time>::pack_masses(const particle_cloud_type& particles, scalar* packed)
        {
          return this->vector_to_array(particles->masses(), packed);
        }

        template<typename scalar, typename time>
        size_t
        WrapperSimplePhysicsSolver<scalar, time>::pack_all(const particle_cloud_type& particles,
                                                           scalar* packed_positions, scalar* packed_velocities,
                                                           scalar* packed_charges, scalar* packed_masses)
        {
          size_t size = this->pack_positions(particles, packed_positions);
          this->pack_velocities(particles, packed_velocities);
          this->pack_charges(particles, packed_charges);
          this->pack_masses(particles, packed_masses);
          return size;
        }

        template<typename scalar, typename time>
        vector<scalar>
        WrapperSimplePhysicsSolver<scalar, time>::unpack_1d(const scalar* packed, const size_t num_particles)
        {
          vector<scalar> out(num_particles);
          for (size_t p = 0; p < num_particles; ++p) {
            out[p] = packed[p];
          }
          return out;
        }

        template<typename scalar, typename time>
        ParticleCloudComponent<scalar>
        WrapperSimplePhysicsSolver<scalar, time>::unpack_2d(const scalar* packed,
                                                            const size_t num_particles)
        {
          return this->unpack_1d(packed, num_particles * DIM);
        }


        template<typename scalar, typename time>
        WrapperSimplePhysicsSolver<scalar, time>::WrapperSimplePhysicsSolver()
        {}

        template<typename scalar, typename time>
        WrapperSimplePhysicsSolver<scalar, time>::~WrapperSimplePhysicsSolver()
        {}

        template<typename scalar, typename time>
        ParticleCloudComponent<scalar>
        WrapperSimplePhysicsSolver<scalar, time>::external_e_field_evaluate(const particle_cloud_type& particles, const time t)
        {
          ML_CVLOG(6, "SolverBinding", "evaluating external E-Field at t=" << t);
          size_t num_particles = particles->size();
          assert(DIM == particles->dim());

          scalar* packed_positions = new scalar[num_particles * DIM];
          scalar* packed_masses = new scalar[num_particles];
          scalar* packed_charges = new scalar[num_particles];
          this->pack_positions(particles, packed_positions);
          this->pack_masses(particles, packed_masses);
          this->pack_charges(particles, packed_charges);

          scalar* packed_forces = new scalar[num_particles * DIM];
          solver::evaluate_external_e_field(packed_positions, packed_charges, packed_masses, num_particles, t,
                                            this->config.get(), packed_forces);

          delete[] packed_positions;
          delete[] packed_masses;
          delete[] packed_charges;
          auto forces = this->unpack_2d(packed_forces, num_particles);
          delete[] packed_forces;
          return forces;
        }

        template<typename scalar, typename time>
        ParticleCloudComponent<scalar>
        WrapperSimplePhysicsSolver<scalar, time>::e_field_evaluate(const particle_cloud_type& particles, const time t)
        {
          ML_CVLOG(6, "SolverBinding", "evaluating complete E-Field at t=" << t);
          size_t num_particles = particles->size();
          assert(DIM == particles->dim());

          scalar* packed_positions = new scalar[num_particles * DIM];
          scalar* packed_masses = new scalar[num_particles];
          scalar* packed_charges = new scalar[num_particles];
          this->pack_positions(particles, packed_positions);
          this->pack_masses(particles, packed_masses);
          this->pack_charges(particles, packed_charges);

          scalar* packed_forces = new scalar[num_particles * DIM];
          solver::evaluate_e_field(packed_positions, packed_charges, packed_masses, num_particles, t,
                                   this->config.get(), packed_forces);

          delete[] packed_positions;
          delete[] packed_masses;
          delete[] packed_charges;
          auto forces = this->unpack_2d(packed_forces, num_particles);
          delete[] packed_forces;
          return forces;
        }

        template<typename scalar, typename time>
        ParticleCloudComponent<scalar>
        WrapperSimplePhysicsSolver<scalar, time>::b_field_evaluate(const particle_cloud_type& particles, const time t)
        {
          ML_CVLOG(6, "SolverBinding", "evaluating B-Field at t=" << t);
          size_t num_particles = particles->size();
          assert(DIM == particles->dim());

          scalar* packed_velocities = new scalar[num_particles * DIM];
          scalar* packed_masses = new scalar[num_particles];
          scalar* packed_charges = new scalar[num_particles];
          this->pack_velocities(particles, packed_velocities);
          this->pack_masses(particles, packed_masses);
          this->pack_charges(particles, packed_charges);

          scalar* packed_forces = new scalar[num_particles * DIM];
          solver::evaluate_b_field(packed_velocities, packed_charges, packed_masses, num_particles, t,
                                   this->config.get(), packed_forces);

          delete[] packed_velocities;
          delete[] packed_masses;
          delete[] packed_charges;
          auto forces = this->unpack_2d(packed_forces, num_particles);
          delete[] packed_forces;
          return forces;
        }

        template<typename scalar, typename time>
        ParticleCloudComponent<scalar>
        WrapperSimplePhysicsSolver<scalar, time>::b_field_vecs(const particle_cloud_type& particles, const time t)
        {
          auto b_vecs = cloud_component_factory<scalar>(particles->size(), particles->dim());
          for (size_t p = 0; p < particles->size(); ++p) {
            auto bvec = this->get_b_field_vector() / particles->charges()[p] / particles->masses()[p];
            std::copy(bvec.cbegin(), bvec.cend(), b_vecs.begin() + (p * DIM));
          }
          return b_vecs;
        }

        template<typename scalar, typename time>
        ParticleCloudComponent<scalar>
        WrapperSimplePhysicsSolver<scalar, time>::force_evaluate(const particle_cloud_type& particles, const time t)
        {
          ML_CVLOG(6, "SolverBinding", "compute total force at t=" << t);
          auto e_force = this->e_field_evaluate(particles, t);
          auto b_force = this->b_field_evaluate(particles, t);
          return e_force + b_force;
        }

        template<typename scalar, typename time>
        scalar
        WrapperSimplePhysicsSolver<scalar, time>::energy(const particle_cloud_type& particles, const time t)
        {
          ML_CVLOG(6, "SolverBinding", "computing system's total energy at t=" << t);
          size_t num_particles = particles->size();
          assert(DIM == particles->dim());

          scalar* packed_velocities = new scalar[num_particles * DIM];
          scalar* packed_positions = new scalar[num_particles * DIM];
          scalar* packed_masses = new scalar[num_particles];
          scalar* packed_charges = new scalar[num_particles];
          this->pack_positions(particles, packed_positions);
          this->pack_velocities(particles, packed_velocities);
          this->pack_masses(particles, packed_masses);
          this->pack_charges(particles, packed_charges);

          auto energy = solver::compute_energy(packed_positions, packed_velocities, packed_charges, packed_masses,
                                               num_particles, t, this->config.get());
          delete[] packed_velocities;
          delete[] packed_positions;
          delete[] packed_masses;
          delete[] packed_charges;
          return energy;
        }

        template<typename scalar, typename time>
        ParticleComponent<scalar> WrapperSimplePhysicsSolver<scalar, time>::get_b_field_vector()
        {
          scalar* packed_vec = new scalar[DIM];
          solver::get_b_field_vector(this->config.get(), packed_vec);
          auto vec = this->unpack_1d(packed_vec, DIM);
          delete[] packed_vec;
          return vec;
        }

        template<typename scalar, typename time>
        void
        WrapperSimplePhysicsSolver<scalar, time>::set_config(shared_ptr<solver::SimplePhysicsSolverConfig> config)
        {
          this->config = config;
        }

        template<typename scalar, typename time>
        scalar
        WrapperSimplePhysicsSolver<scalar, time>::omega_e() const
        {
          return this->config->omega_e;
        }

        template<typename scalar, typename time>
        scalar
        WrapperSimplePhysicsSolver<scalar, time>::omega_b() const
        {
          return this->config->omega_b;
        }

        template<typename scalar, typename time>
        scalar
        WrapperSimplePhysicsSolver<scalar, time>::epsilon() const
        {
          return this->config->epsilon;
        }

        template<typename precision, typename time>
        void WrapperSimplePhysicsSolver<precision, time>::log(el::base::type::ostream_t& os) const
        {
          os << "WrapperSimplePhysicsSolver()";
        }


        template<typename scalar, typename time>
        void setup(shared_ptr<WrapperSimplePhysicsSolver<scalar, time>> wrapper)
        {
          wrapper->set_config(make_shared<solver::SimplePhysicsSolverConfig>());
        }

        template<typename scalar, typename time, typename ArgT>
        void setup(shared_ptr<WrapperSimplePhysicsSolver<scalar, time>> wrapper, ArgT arg)
        {
          UNUSED(wrapper); UNUSED(arg);
        }

        template<typename scalar, typename time, typename ArgT, typename... ArgsT>
        void setup(shared_ptr<WrapperSimplePhysicsSolver<scalar, time>> wrapper, ArgT arg, ArgsT... args)
        {
          setup(wrapper, arg);
          setup(wrapper, args...);
        }
      }  // ::pfasst::examples::boris::bindings
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst
