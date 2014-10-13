#ifndef _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_
#define _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_

#include <Eigen/Dense>

#include <pfasst.hpp>

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
    scalar epsilon;

  public:
    /**
     * Using \\( \omega_E = 4.9 \\) as in WSR-Boris-SDC-Paper
     */
    IdealQuadrupolePotential()
      : ElectricField<scalar, time, ParticleT>(scalar(-4.9))
    {
      this->matrix << scalar(1.0), scalar(0.0), scalar(0.0),
                      scalar(0.0), scalar(1.0), scalar(0.0),
                      scalar(0.0), scalar(0.0), scalar(-2.0);
      this->epsilon = scalar(-1.0);
    }
    IdealQuadrupolePotential(const scalar epsilon)
      :   IdealQuadrupolePotential()
    {
      this->epsilon = epsilon;
    }
    virtual ~IdealQuadrupolePotential()
    {}

    /**
     * Calculates the potential energy of a particle in a ideal quadrupole potential.
     *
     * @f[
     *   E(\vec{x}_m) = - \frac{\varepsilon \omega_z^2}{\alpha} \begin{pmatrix}1 & 0 \\ 0 & 1\end{pmatrix} \vec{x}
     * @f]
     * 
     * @note The inter-particle Coulomb interactions (i.e. the inner electric field) is not
     *     implemented yet.
     */
    virtual typename parent_type::particle_type::acceleration_type
    evaluate(vector<shared_ptr<typename parent_type::particle_type>> particles,
             size_t m,
             time t) const override
    {
      UNUSED(t);
      Matrix<scalar> pos_m(particles[m]->pos().as_matrix().transpose());
      scalar factor = (- this->epsilon) * (this->omega_e * this->omega_e) / particles[m]->alpha();
      Matrix<scalar> accel(factor * this->matrix * pos_m);
      return typename parent_type::particle_type::acceleration_type(accel);
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

    virtual Matrix<time> get_field_vector() const
    {
      Matrix<time> vec(1,3);
      vec.fill(time(0.0));
      vec(0, 2) = time(1.0);
      vec *= this->omega_b;
      return vec;
    }

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
             time t) const override
    {
      UNUSED(t);
      Matrix<scalar> vel_m(particles[m]->vel().as_matrix().transpose());
      scalar factor = this->omega_b / particles[m]->alpha();
      Matrix<scalar> accel(factor * this->matrix * vel_m);
      return typename parent_type::particle_type::acceleration_type(accel);
    }
};


template<
  typename scalar,
  typename time,
  template <typename, typename> class ParticleT,
  template <typename, typename, template <typename, typename> class> class EFieldT,
  template <typename, typename, template <typename, typename> class> class BFieldT
>
class SimplePhysicsEnergyOperator
  : public EnergyOperator<scalar, time, ParticleT, EFieldT, BFieldT>
{
  private:
    typedef SimplePhysicsEnergyOperator<scalar, time, ParticleT, EFieldT, BFieldT> this_type;
    typedef EnergyOperator<scalar, time, ParticleT, EFieldT, BFieldT> parent_type;

  public:
    //! @{
    typedef ParticleT<scalar, time> particle_type;
    typedef EFieldT<scalar, time, ParticleT> e_field_type;
    typedef BFieldT<scalar, time, ParticleT> b_field_type;
    //! @}

  protected:
    Matrix<time> op;
    scalar epsilon;

  public:
    //! @{
    SimplePhysicsEnergyOperator()
    {}
    SimplePhysicsEnergyOperator(const e_field_type& e_field,
                                const b_field_type& b_field,
                                scalar episolon = -1.0)
      :   parent_type(e_field, b_field)
        , epsilon(episolon)
    {
      this->op = Matrix<time>(6, 6);
      this->op.fill(time(0.0));
      this->op(0, 0) = this->epsilon * pow(this->e_field.omega_e, 2);
      this->op(1, 1) = this->epsilon * pow(this->e_field.omega_e, 2);
      this->op(2, 2) = - 2 * this->epsilon * pow(this->e_field.omega_e, 2);
      this->op(3, 3) = 1.0;
      this->op(4, 4) = 1.0;
      this->op(5, 5) = 1.0;
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
      scalar energy(0.0);
      Matrix<time> u_vec = Matrix<time>(6, 1);
      for (size_t i = 0; i < particles.size(); ++i) {
        auto p = particles[i];
        u_vec.block(0, 0, 3, 1) = p->pos().as_matrix().transpose();
        u_vec.block(3, 0, 3, 1) = p->vel().as_matrix().transpose();
        auto e = Matrix<time>(u_vec.transpose() * (p->mass() / time(2.0) * this->op) * u_vec);
        energy += e(0);
      }
      return energy;
    }
    //! @}
};

static_assert(std::is_move_constructible<SimplePhysicsEnergyOperator<double, double, Particle3DEncapsulation, IdealQuadrupolePotential, ConstantMagneticField>>::value, "");
static_assert(std::is_move_assignable<SimplePhysicsEnergyOperator<double, double, Particle3DEncapsulation, IdealQuadrupolePotential, ConstantMagneticField>>::value, "");


#endif  // _EXAMPLES__BORIS__SIMPLE_PHYSICS_HPP_
