#ifndef _EXAMPLES__BORIS__PHYSICS__HPP_
#define _EXAMPLES__BORIS__PHYSICS__HPP_

#include <memory>
#include <functional>

#include <pfasst/globals.hpp>
#include <pfasst/interfaces.hpp>

#include "particle.hpp"

using namespace std;
using namespace pfasst;


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
    scalar omega_z;
    //! @}

  public:
    //! @{
    ElectricField(scalar omega_z)
      : omega_z(omega_z)
    {}

    ElectricField()
      : ElectricField(scalar(1.0))
    {}

    virtual ~ElectricField()
    {}
    //! @}

    //! @{
    /**
     * Calculates the potential energy \\( E(\\vec{x},t) \\) of a particle at a given time \\( t \\).
     */
    virtual typename parent_type::particle_type::acceleration_type
    evaluate(shared_ptr<const typename parent_type::particle_type> particle, time t)
    {
      UNUSED(particle); UNUSED(t);
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
    scalar omega_c;
    //! @}

  public:
    //! @{
    MagneticField(scalar omega_c)
      : omega_c(omega_c)
    {}

    MagneticField()
      : MagneticField(scalar(1.0))
    {}

    virtual ~MagneticField()
    {}
    //! @}

    //! @{
    /**
     * Calculates the kinetic energy \\( \\vec{v} \\times B(t) \\) of a particle at a given time 
     * \\( t \\).
     */
    virtual typename parent_type::particle_type::acceleration_type
    evaluate(shared_ptr<const typename parent_type::particle_type> particle, time t)
    {
      UNUSED(particle); UNUSED(t);
      throw NotImplementedYet("evaluate of MagneticField");
    }
    //! @}
};


template<
  typename scalar,
  typename time,
  typename ParticleT,
  typename EFieldT,
  typename BFieldT
>
class EnergyOperator
{
  public:
    //! @{
    typedef ParticleT particle_type;
    typedef EFieldT e_field_type;
    typedef BFieldT b_field_type;
    //! @}

  public:
    //! @{
    EnergyOperator()
    {}

    virtual ~EnergyOperator()
    {}
    //! @}

    //! @{
    /**
     * Computes the energy of a particle.
     *
     * The energy of the particle \\( p \\) in a given electric field \\( E(p,t) \\) and 
     * magnetic field \\( B(t) \\) at a given time \\( t \\) is usually given by the sum of the 
     * potential and kinetic energy.
     */
    virtual scalar evaluate(shared_ptr<const particle_type> particle,
                            time t,
                            shared_ptr<const e_field_type> E,
                            shared_ptr<const b_field_type> B)
    {
      UNUSED(particle); UNUSED(t); UNUSED(E); UNUSED(B);
      throw NotImplementedYet("evaluate of Energy Operator");
      return scalar(0.0);
    }
    //! @}
};


#endif  // _EXAMPLES__BORIS__PHYSICS__HPP_
