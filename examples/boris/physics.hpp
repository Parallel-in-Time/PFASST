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

    //! @{
    virtual scalar evaluate(shared_ptr<const particle_type> particle, time t) = 0;
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
  private:
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
    virtual scalar evaluate(shared_ptr<const typename parent_type::particle_type> particle, time t)
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
  private:
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
    virtual scalar evaluate(shared_ptr<const typename parent_type::particle_type> particle, time t)
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
    virtual scalar evaluate(shared_ptr<const particle_type> p,
                            time t,
                            shared_ptr<const e_field_type> E,
                            shared_ptr<const b_field_type> B)
    {
      UNUSED(p); UNUSED(t); UNUSED(E); UNUSED(B);
      throw NotImplementedYet("evaluate of Energy Operator");
      return scalar(0.0);
    }
    //! @}
};


#endif  // _EXAMPLES__BORIS__PHYSICS__HPP_
