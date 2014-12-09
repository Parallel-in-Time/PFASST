#include "wrapper_simple_physics_solver.hpp"

#include "simple_physics_solver.hpp"
namespace solver = simple_physics_solver;

#define UNUSED(expr) (void)(expr)

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
      WrapperSimplePhysicsSolver<scalar, time>::vector2d_to_array(const vector<vector<scalar>>& vec,
                                                                  scalar* arr)
      {
        for (size_t p = 0; p < vec.size(); ++p) {
          for (size_t d = 0; d < DIM; ++d) {
            arr[p * DIM + d] = vec[p][d];
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
        ParticleCloudComponent<scalar> out(num_particles);
        for (size_t p = 0; p < num_particles; ++p) {
          out[p] = this->unpack_1d(packed + p, DIM);
        }
        return out;
      }


      template<typename scalar, typename time>
      WrapperSimplePhysicsSolver<scalar, time>::WrapperSimplePhysicsSolver()
      {}

      template<typename scalar, typename time>
      WrapperSimplePhysicsSolver<scalar, time>::~WrapperSimplePhysicsSolver()
      {}

      template<typename scalar, typename time>
      ParticleCloudComponent<scalar>
      WrapperSimplePhysicsSolver<scalar, time>::e_field_evaluate(const particle_cloud_type& particles, const time t)
      {
        size_t num_particles = particles->size();
        assert(DIM == particles->dim());

        scalar packed_positions[num_particles * DIM];
        scalar packed_masses[num_particles];
        scalar packed_charges[num_particles];
        this->pack_positions(particles, packed_positions);
        this->pack_masses(particles, packed_masses);
        this->pack_charges(particles, packed_charges);

        scalar packed_forces[num_particles * DIM];
        solver::evaluate_e_field(packed_positions, packed_charges, packed_masses, num_particles, t,
                                 this->config.get(), packed_forces);

        auto out = this->unpack_2d(packed_forces, num_particles);
//         delete[] packed_positions;
//         delete[] packed_masses;
//         delete[] packed_charges;
//         delete[] packed_forces;
        return out;
      }

      template<typename scalar, typename time>
      ParticleCloudComponent<scalar>
      WrapperSimplePhysicsSolver<scalar, time>::b_field_evaluate(const particle_cloud_type& particles, const time t)
      {
        size_t num_particles = particles->size();
        assert(DIM == particles->dim());

        scalar packed_velocities[num_particles * DIM];
        scalar packed_masses[num_particles];
        scalar packed_charges[num_particles];
        this->pack_velocities(particles, packed_velocities);
        this->pack_masses(particles, packed_masses);
        this->pack_charges(particles, packed_charges);

        scalar packed_forces[num_particles * DIM];
        solver::evaluate_b_field(packed_velocities, packed_charges, packed_masses, num_particles, t,
                                 this->config.get(), packed_forces);

        auto out = this->unpack_2d(packed_forces, num_particles);
//         delete[] packed_velocities;
//         delete[] packed_masses;
//         delete[] packed_charges;
//         delete[] packed_forces;
        return out;
      }

      template<typename scalar, typename time>
      ParticleCloudComponent<scalar>
      WrapperSimplePhysicsSolver<scalar, time>::force_evaluate(const particle_cloud_type& particles, const time t)
      {
        auto e_force = this->e_field_evaluate(particles, t);
        auto b_force = this->b_field_evaluate(particles, t);
        return e_force + b_force;
      }

      template<typename scalar, typename time>
      scalar
      WrapperSimplePhysicsSolver<scalar, time>::energy(const particle_cloud_type& particles, const time t)
      {
        size_t num_particles = particles->size();
        assert(DIM == particles->dim());

        scalar packed_velocities[num_particles * DIM];
        scalar packed_positions[num_particles * DIM];
        scalar packed_masses[num_particles];
        scalar packed_charges[num_particles];
        this->pack_positions(particles, packed_positions);
        this->pack_velocities(particles, packed_velocities);
        this->pack_masses(particles, packed_masses);
        this->pack_charges(particles, packed_charges);

        scalar energy = solver::compute_energy(packed_positions, packed_velocities, packed_charges, packed_masses,
                                               num_particles, t, this->config.get());

//         delete[] packed_positions;
//         delete[] packed_velocities;
//         delete[] packed_masses;
//         delete[] packed_charges;
        return energy;
      }

      template<typename scalar, typename time>
      ParticleComponent<scalar> WrapperSimplePhysicsSolver<scalar, time>::get_b_field_vector()
      {
        scalar packed_vec[DIM];
        solver::get_b_field_vector(this->config.get(), packed_vec);
        return this->unpack_1d(packed_vec, DIM);
      }

      template<typename scalar, typename time>
      void
      WrapperSimplePhysicsSolver<scalar, time>::set_config(shared_ptr<solver::SimplePhysicsSolverConfig>& config)
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
    }  // ::examples::boris::bindings
  }  // ::examples::boris
}  // ::examples
