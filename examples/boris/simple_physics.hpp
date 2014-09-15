#ifndef _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_
#define _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_

#include "physics.hpp"
#include "particle_3d.hpp"

template<
  typename scalar,
  typename time,
  template <typename, typename> class ParticleT
>
class IdealQuadrupolePotential
  : public ElectricField<scalar, time, ParticleT>
{
  private:
    typedef ElectricField<scalar, time, ParticleT> parent_type;

  public:
    virtual ~IdealQuadrupolePotential()
    {}

    /**
     * Calculates the potential energy of a particle in a ideal quadrupole potential.
     *
     * @f[
     *   E(\vec{x}_m) = - \frac{\varepsilon \omega_z^2}{\alpha} \begin{pmatrix}1 & 0 \\ 0 & 1\end{pmatrix} \vec{x}
     * @f]
     */
    virtual scalar evaluate(vector<shared_ptr<typename parent_type::particle_type>> particles,
                            size_t m,
                            time t) override
    {}
};


template<
  typename scalar,
  typename time,
  template <typename, typename> class ParticleT
>
class ConstantMagneticField
  : public MagneticField<scalar, time, ParticleT>
{
  private:
    typedef MagneticField<scalar, time, ParticleT> parent_type;

  public:
    virtual ~ConstantMagneticField()
    {}

    /**
     * Calculates the kinetic energy of a particle in a constant magnetic field.
     *
     * As the magnetic field is assumed constant in space and time, evaluation can be simplified to
     * @f[
     *   \vec{v}_m \times B = \frac{\omega_c}{\alpha} \begin{pmatrix}0 & 1 \\ -1 & 0\end{pmatrix} \vec{v}
     * @f]
     */
    virtual scalar evaluate(vector<shared_ptr<typename parent_type::particle_type>> particles,
                            size_t m,
                            time t) override
    {}
};

#endif  // _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_
