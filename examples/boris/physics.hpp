#ifndef _EXAMPLES__BORIS__PHYSICS__HPP_
#define _EXAMPLES__BORIS__PHYSICS__HPP_

#include <memory>
#include <functional>
using namespace std;

#include <pfasst/globals.hpp>
#include <pfasst/interfaces.hpp>

#include "particle.hpp"

namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      template<
        typename scalar,
        typename time,
        template <typename, typename> class ParticleT
      >
      class PhysicalField
      {
        public:
          //! @{
          typedef ParticleT<scalar, time> particle_type;
          //! @}

        public:
          //! @{
          virtual ~PhysicalField()
          {}
          //! @}
      };


      template<
        typename scalar,
        typename time,
        template <typename, typename> class ParticleT
      >
      class ElectricField
        : public PhysicalField<scalar, time, ParticleT>
      {
        protected:
          typedef PhysicalField<scalar, time, ParticleT> parent_type;

        public:
          //! @{
          scalar omega_e;
          //! @}

        public:
          //! @{
          ElectricField(scalar omega_e)
            : omega_e(omega_e)
          {}

          ElectricField()
            : ElectricField(scalar(1.0))
          {}

          virtual ~ElectricField()
          {}
          //! @}

          //! @{
          /**
           * Calculates the potential energy \\( E(\\vec{x}_m,t) \\) of a particle \\( m \\) at a given 
           * time \\( t \\) with respect to all particles.
           */
          virtual typename parent_type::particle_type::acceleration_type
          evaluate(vector<shared_ptr<typename parent_type::particle_type>> particles, size_t m,
                   time t) const
          {
            UNUSED(particles); UNUSED(m); UNUSED(t);
            throw NotImplementedYet("evaluate of ElectricField");
          }
          //! @}
      };


      template<
        typename scalar,
        typename time,
        template <typename, typename> class ParticleT
      >
      class MagneticField
        : public PhysicalField<scalar, time, ParticleT>
      {
        protected:
          typedef PhysicalField<scalar, time, ParticleT> parent_type;

        public:
          //! @{
          scalar omega_b;
          //! @}

        public:
          //! @{
          MagneticField(scalar omega_b)
            : omega_b(omega_b)
          {}

          MagneticField()
            : MagneticField(scalar(1.0))
          {}

          virtual ~MagneticField()
          {}
          //! @}

          //! @{
          /**
           * Calculates the kinetic energy \\( \\vec{v}_m \\times B(t) \\) of a particle \\( m \\) at a 
           * given time \\( t \\) with respect to all particles.
           */
          virtual typename parent_type::particle_type::acceleration_type
          evaluate(vector<shared_ptr<typename parent_type::particle_type>> particles, size_t m,
                   time t) const
          {
            UNUSED(particles); UNUSED(m); UNUSED(t);
            throw NotImplementedYet("evaluate of MagneticField");
          }
          //! @}
      };


      template<
        typename scalar,
        typename time,
        template <typename, typename> class ParticleT,
        template <typename, typename, template <typename, typename> class> class EFieldT,
        template <typename, typename, template <typename, typename> class> class BFieldT
      >
      class EnergyOperator
      {
        private:
          typedef EnergyOperator<scalar, time, ParticleT, EFieldT, BFieldT> this_type;

        public:
          //! @{
          typedef ParticleT<scalar, time> particle_type;
          typedef EFieldT<scalar, time, ParticleT> e_field_type;
          typedef BFieldT<scalar, time, ParticleT> b_field_type;
          //! @}

        protected:
          e_field_type e_field;
          b_field_type b_field;

        public:
          //! @{
          EnergyOperator()
          {}
          EnergyOperator(const e_field_type& e_field, const b_field_type& b_field)
            :   e_field(e_field)
              , b_field(b_field)
          {}
          //! @}

          //! @{
          virtual const e_field_type& get_e_field() const
          {
            return this->e_field;
          }
          virtual const b_field_type& get_b_field() const
          {
            return this->b_field;
          }
          //! @}

          //! @{
          /**
           * Computes the energy of a particle.
           *
           * The energy of the particles \\( P \\) in a given electric field \\( E(P,t) \\) and 
           * magnetic field \\( B(t) \\) at a given time \\( t \\) is usually given by the sum of the 
           * potential and kinetic energy.
           */
          virtual scalar evaluate(const vector<shared_ptr<particle_type>>& particles, const time t) const
          {
            UNUSED(particles); UNUSED(t);
            throw NotImplementedYet("evaluate of Energy Operator");
            return scalar(0.0);
          }
          //! @}
      };
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__BORIS__PHYSICS__HPP_
