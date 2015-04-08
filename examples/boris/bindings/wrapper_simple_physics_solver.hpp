/**
 * @ingroup Boris_Bindings
 */
#ifndef _EXAMPLES__BORIS__BINDINGS__WRAPPER_SIMPLE_PHYSICS_SOLVER_HPP_
#define _EXAMPLES__BORIS__BINDINGS__WRAPPER_SIMPLE_PHYSICS_SOLVER_HPP_

#include <cstdlib>
#include <memory>
#include <utility>
using namespace std;

#include "../particle_cloud.hpp"
#include "wrapper_interface.hpp"

#include "simple_physics_solver.hpp"
namespace solver = simple_physics_solver;


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      namespace bindings
      {
        /**
         * @ingroup Boris_Bindings
         */
        template<
          typename scalar,
          typename time
        >
        class WrapperSimplePhysicsSolver
          : public WrapperInterface<scalar, time>
        {
          private:
            shared_ptr<solver::SimplePhysicsSolverConfig> config;

          public:
            typedef shared_ptr<ParticleCloud<scalar>> particle_cloud_type;

          protected:
            virtual size_t vector_to_array(const vector<scalar>& vec, scalar* arr);
            virtual size_t vector2d_to_array(const vector<scalar>& vec, scalar* arr);

            virtual size_t pack_positions(const particle_cloud_type& particles, scalar* packed) override;

            virtual size_t pack_velocities(const particle_cloud_type& particles, scalar* packed) override;

            virtual size_t pack_charges(const particle_cloud_type& particles, scalar* packed) override;

            virtual size_t pack_masses(const particle_cloud_type& particles, scalar* packed) override;

            virtual size_t pack_all(const particle_cloud_type& particles,
                                    scalar* packed_positions, scalar* packed_velocities,
                                    scalar* packed_charges, scalar* packed_masses) override;

            virtual vector<scalar> unpack_1d(const scalar* packed, const size_t num_particles);
            virtual ParticleCloudComponent<scalar> unpack_2d(const scalar* packed,
                                                             const size_t num_particles);

          public:
            WrapperSimplePhysicsSolver();
            virtual ~WrapperSimplePhysicsSolver();

            virtual ParticleCloudComponent<scalar> external_e_field_evaluate(const particle_cloud_type& particles, const time t) override;
            virtual ParticleCloudComponent<scalar> e_field_evaluate(const particle_cloud_type& particles, const time t) override;
            virtual ParticleCloudComponent<scalar> b_field_evaluate(const particle_cloud_type& particles, const time t) override;
            virtual ParticleCloudComponent<scalar> force_evaluate(const particle_cloud_type& particles, const time t) override;
            virtual scalar energy(const particle_cloud_type& particles, const time t) override;

            virtual ParticleComponent<scalar> get_b_field_vector() override;
            virtual ParticleCloudComponent<scalar> b_field_vecs(const particle_cloud_type& particles, const time t) override;

            virtual void set_config(shared_ptr<solver::SimplePhysicsSolverConfig> config);
            virtual scalar omega_b() const;
            virtual scalar omega_e() const;
            virtual scalar epsilon() const;

            virtual void log(el::base::type::ostream_t& os) const;
        };


        //! @ingroup Boris_Bindings
        template<typename scalar, typename time>
        void setup(shared_ptr<WrapperSimplePhysicsSolver<scalar, time>> wrapper);

        //! @ingroup Boris_Bindings
        template<typename scalar, typename time, typename ArgT>
        void setup(shared_ptr<WrapperSimplePhysicsSolver<scalar, time>> wrapper, ArgT arg);

        //! @ingroup Boris_Bindings
        template<typename scalar, typename time, typename ArgT, typename... ArgsT>
        void setup(shared_ptr<WrapperSimplePhysicsSolver<scalar, time>> wrapper, ArgT arg, ArgsT... args);
      }  // ::pfasst::examples::boris::bindings
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#include "wrapper_simple_physics_solver_impl.hpp"

#endif  // _EXAMPLES__BORIS__BINDINGS__WRAPPER_SIMPLE_PHYSICS_SOLVER_HPP_
