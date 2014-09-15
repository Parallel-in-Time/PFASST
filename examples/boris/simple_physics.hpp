#ifndef _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_
#define _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_

#include <Eigen/Dense>

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

  protected:
    Eigen::Matrix<scalar, 3, 3> matrix;

  public:
    /**
     * Using \\( \omega_E = 4.9 \\) as in WSR-Boris-SDC-Paper
     */
    IdealQuadrupolePotential()
      : ElectricField<scalar, time, ParticleT>(scalar(4.9))
    {
      this->matrix << scalar(1.0), scalar(0.0), scalar(0.0),
                      scalar(0.0), scalar(1.0), scalar(0.0),
                      scalar(0.0), scalar(0.0), scalar(-2.0);
    }
    virtual ~IdealQuadrupolePotential()
    {}

    /**
     * Calculates the potential energy of a particle in a ideal quadrupole potential.
     *
     * @f[
     *   E(\vec{x}_m) = - \frac{\varepsilon \omega_z^2}{\alpha} \begin{pmatrix}1 & 0 \\ 0 & 1\end{pmatrix} \vec{x}
     * @f]
     */
    virtual typename parent_type::particle_type::acceleration_type
    evaluate(vector<shared_ptr<typename parent_type::particle_type>> particles,
             size_t m,
             time t) override
    {
      UNUSED(t);
      typename parent_type::particle_type::position_type pos_m = particles[m]->pos();
      scalar a = this->matrix(0,0) * pos_m.x + this->matrix(0,1) * pos_m.y + this->matrix(0,2) * pos_m.z;
      scalar b = this->matrix(1,0) * pos_m.x + this->matrix(1,1) * pos_m.y + this->matrix(1,2) * pos_m.z;
      scalar c = this->matrix(2,0) * pos_m.x + this->matrix(2,1) * pos_m.y + this->matrix(2,2) * pos_m.z;
      scalar factor = (this->omega_e * this->omega_e) / particles[m]->alpha();
      return typename parent_type::particle_type::acceleration_type(factor * a, factor * b, factor * c);
    }
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

  protected:
    Eigen::Matrix<scalar, 3, 3> matrix;

  public:
    /**
     * Using \\( \omega_B = 25.0 \\) as in WSR-Boris-SDC-Paper
     */
    ConstantMagneticField()
      : MagneticField<scalar, time, ParticleT>(scalar(25.0))
    {
      this->matrix << scalar(0.0),  scalar(1.0), scalar(0.0),
                      scalar(-1.0), scalar(0.0), scalar(0.0),
                      scalar(0.0),  scalar(0.0), scalar(0.0);
    }
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
    virtual typename parent_type::particle_type::acceleration_type
    evaluate(vector<shared_ptr<typename parent_type::particle_type>> particles,
             size_t m,
             time t) override
    {
      UNUSED(t);
      typename parent_type::particle_type::velocity_type vel_m = particles[m]->vel();
      scalar a = this->matrix(0,0) * vel_m.u + this->matrix(0,1) * vel_m.v + this->matrix(0,2) * vel_m.w;
      scalar b = this->matrix(1,0) * vel_m.u + this->matrix(1,1) * vel_m.v + this->matrix(1,2) * vel_m.w;
      scalar c = this->matrix(2,0) * vel_m.u + this->matrix(2,1) * vel_m.v + this->matrix(2,2) * vel_m.w;
      scalar factor = this->omega_b / particles[m]->alpha();
      return typename parent_type::particle_type::acceleration_type(factor * a, factor * b, factor * c);
    }
};

#endif  // _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_
