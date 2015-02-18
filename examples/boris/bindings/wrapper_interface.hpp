/**
 * @defgroup BorisBindings Bindings
 * @ingroup Boris
 *
 * @file examples/boris/particle_util.hpp
 */
#ifndef _EXAMPLES__BORIS__BINDINGS__WRAPPER_INTERFACE_HPP_
#define _EXAMPLES__BORIS__BINDINGS__WRAPPER_INTERFACE_HPP_

#include <cstdlib>
#include <memory>
#include <utility>
using namespace std;

#include <pfasst/logging.hpp>

#include "../particle_cloud.hpp"

#define UNUSED(expr) (void)(expr)


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      /**
       * @ingroup BorisBindings
       */
      namespace bindings
      {
        /**
         * @ingroup BorisBindings
         */
        template<
          typename scalar,
          typename time
        >
        class WrapperInterface
          : public el::Loggable
        {
          public:
            typedef shared_ptr<ParticleCloud<scalar>> particle_cloud_type;

          protected:
            virtual size_t pack_positions(const particle_cloud_type& particles, scalar* packed) = 0;
            virtual size_t pack_velocities(const particle_cloud_type& particles, scalar* packed) = 0;
            virtual size_t pack_charges(const particle_cloud_type& particles, scalar* packed) = 0;
            virtual size_t pack_masses(const particle_cloud_type& particles, scalar* packed) = 0;
            virtual size_t pack_all(const particle_cloud_type& particles,
                                    scalar* packed_positions, scalar* packed_velocities,
                                    scalar* packed_charges, scalar* packed_masses) = 0;

          public:
            virtual ~WrapperInterface() {}

            virtual ParticleCloudComponent<scalar>
            external_e_field_evaluate(const particle_cloud_type& particles, const time t) = 0;

            virtual ParticleCloudComponent<scalar>
            e_field_evaluate(const particle_cloud_type& particles, const time t) = 0;

            virtual ParticleCloudComponent<scalar>
            b_field_evaluate(const particle_cloud_type& particles, const time t) = 0;

            virtual ParticleCloudComponent<scalar>
            force_evaluate(const particle_cloud_type& particles, const time t) = 0;

            virtual ParticleComponent<scalar> get_b_field_vector() = 0;
            virtual ParticleCloudComponent<scalar> b_field_vecs(const particle_cloud_type& particles, const time t) = 0;

            virtual scalar energy(const particle_cloud_type& particles, const time t) = 0;

            virtual scalar omega_b() const = 0;
            virtual scalar omega_e() const = 0;
            virtual scalar epsilon() const = 0;

            virtual void log(el::base::type::ostream_t& os) const = 0;
        };

        //! @ingroup BorisBindings
        template<typename scalar, typename time>
        void setup(shared_ptr<WrapperInterface<scalar, time>> wrapper)
        {
          UNUSED(wrapper);
        }

        //! @ingroup BorisBindings
        template<typename scalar, typename time, typename ArgT>
        void setup(shared_ptr<WrapperInterface<scalar, time>> wrapper, ArgT arg)
        {
          UNUSED(wrapper); UNUSED(arg);
        }

        //! @ingroup BorisBindings
        template<typename scalar, typename time, typename ArgT, typename... ArgsT>
        void setup(shared_ptr<WrapperInterface<scalar, time>> wrapper, ArgT arg, ArgsT... args)
        {
          setup(wrapper, arg);
          setup(wrapper, args...);
        }
      }  // ::pfasst::examples::boris::bindings
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__BORIS__BINDINGS__WRAPPER_INTERFACE_HPP_
