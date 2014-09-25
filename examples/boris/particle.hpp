#ifndef _EXAMPLES__BORIS__PARTICLE__HPP_
#define _EXAMPLES__BORIS__PARTICLE__HPP_

#include <cstdlib>
#include <stdexcept>
#include <memory>
#include <cassert>

#include <Eigen/Dense>

#include <pfasst/globals.hpp>
#include <pfasst/interfaces.hpp>
#include <pfasst/encap/encapsulation.hpp>

using namespace std;
using namespace pfasst;
using namespace pfasst::encap;

template<typename coeff>
using Matrix = Eigen::Matrix<coeff, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

typedef false_type deactivated_function;


/**
 * Type to denote a scalar with physical unit \\( t \\).
 *
 * This is used where we use a plain scalar to be multiplied with a unit of \\( \\frac{\cdot}{t^p} \\)
 * to be *converted* into a value with unit \\( \\frac{\cdot}{t^{p-1}} \\).
 *
 * @tparam scalar underlying fundamental type
 */
template<typename scalar>
struct dt
{
  scalar v;

  dt() : dt(scalar(0.0)) {}
  dt(const scalar value) : v(value) {}
  virtual ~dt() {}
};


/**
 * Type to denote a scalar with physical unit \\( t^2 \\).
 *
 * This is used where we use a plain scalar to be multiplied with a unit of \\( \\frac{\cdot}{t^p} \\)
 * to be *converted* into a value with unit \\( \\frac{\cdot}{t^{p-2}} \\).
 *
 * @tparam scalar underlying fundamental type
 */
template<typename scalar>
struct dtdt
{
  scalar v;

  dtdt() : dtdt(scalar(0.0)) {}
  dtdt(const scalar value) : v(value) {}
  virtual ~dtdt() {}
};


template<
  typename scalar,
  typename time = time_precision
>
class ParticleComponentEncapsulation
  : public Encapsulation<time>
{
  public:
    //! @{
    const size_t DIM = 0;
    //! @}

  public:
    //! @{
    virtual ~ParticleComponentEncapsulation()
    {}
    //! @}

    //! @{
    virtual void copy(shared_ptr<const Encapsulation<time>> other) override
    {
      UNUSED(other);
      assert(deactivated_function::value);
    }

    virtual Matrix<scalar> as_matrix() const = 0;
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      UNUSED(a); UNUSED(x);
      throw NotImplementedYet("ax+y for Encapsulation onto ParticleComponent");
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, Matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      throw NotImplementedYet("aA*x for Encapsulation onto ParticleComponent");
    }
    //! @}

    //! @{
    virtual scalar operator[](const size_t index) const
    {
      UNUSED(index);
      throw NotImplementedYet("operator[] for ParticleComponent");
    }
    //! @}
};


template<
  typename scalar,
  typename time = time_precision
>
class PositionEncapsulation
  : public ParticleComponentEncapsulation<scalar, time>
{
  protected:
    typedef ParticleComponentEncapsulation<scalar, time> parent_type;

  public:
    //! @{
    virtual ~PositionEncapsulation()
    {}
    //! @}
};


template<
  typename scalar,
  typename time = time_precision
>
class VelocityEncapsulation
  : public ParticleComponentEncapsulation<scalar, time>
{
  protected:
    typedef ParticleComponentEncapsulation<scalar, time> parent_type;

  public:
    //! @{
    virtual ~VelocityEncapsulation()
    {}
    //! @}
};


template<
  typename scalar,
  typename time = time_precision
>
class AccelerationEncapsulation
  : public ParticleComponentEncapsulation<scalar, time>
{
  protected:
    typedef ParticleComponentEncapsulation<scalar, time> parent_type;

  public:
    //! @{
    virtual ~AccelerationEncapsulation()
    {}
    //! @}
};


template<
  typename scalar,
  typename time = time_precision,
  template <typename, typename> class PositionT = PositionEncapsulation,
  template <typename, typename> class VelocityT = VelocityEncapsulation,
  template <typename, typename> class AccelerationT = AccelerationEncapsulation
>
class ParticleEncapsulation
  : public Encapsulation<time>
{
  public:
    //! @{
    const size_t DIM = 0;
    typedef PositionT<scalar, time> position_type;
    typedef VelocityT<scalar, time> velocity_type;
    typedef AccelerationT<scalar, time> acceleration_type;
    typedef ParticleEncapsulation<scalar, time, PositionT, VelocityT, AccelerationT> this_type;
    //! @}

  protected:
    //! @{
    scalar m_mass;
    scalar m_charge;
    scalar m_alpha;
    position_type m_pos;
    velocity_type m_vel;
    acceleration_type m_accel;
    //! @}

  public:
    //! @{
    ParticleEncapsulation()
      : ParticleEncapsulation(scalar(1.0), scalar(1.0))
    {}

    ParticleEncapsulation(const scalar mass, const scalar charge)
      :   m_mass(mass)
        , m_charge(charge)
        , m_alpha(mass / charge)
    {}

    ParticleEncapsulation(const scalar mass, const scalar charge,
                          const position_type& pos, const velocity_type& vel,
                          const acceleration_type& accel)
      :   m_mass(mass)
        , m_charge(charge)
        , m_alpha(mass / charge)
        , m_pos(pos)
        , m_vel(vel)
        , m_accel(accel)
    {}

    ParticleEncapsulation(const this_type& other)
      : ParticleEncapsulation(other.mass(), other.charge(),
                              other.pos(), other.vel(), other.accel())
    {}

    ParticleEncapsulation(this_type&& other)
      : ParticleEncapsulation(std::move(other.m_mass), std::move(other.m_charge),
                              std::move(other.m_pos), std::move(other.m_vel),
                              std::move(other.m_accel))
    {}

    virtual ~ParticleEncapsulation()
    {}
    //! @}

    //! @{
    this_type& operator=(const this_type& other)
    {
      assert(this != &other);
      this_type tmp(other);
      std::swap(this->m_mass, tmp.m_mass);
      std::swap(this->m_charge, tmp.m_charge);
      std::swap(this->m_pos, tmp.m_pos);
      std::swap(this->m_vel, tmp.m_vel);
      std::swap(this->m_accel, tmp.m_accel);
      return *this;
    }

    this_type& operator=(this_type&& other)
    {
      std::swap(this->m_mass, other.m_mass);
      std::swap(this->m_charge, other.m_charge);
      std::swap(this->m_pos, other.m_pos);
      std::swap(this->m_vel, other.m_vel);
      std::swap(this->m_accel, other.m_accel);
    }
    //! @}

    //! @{
    virtual void copy(shared_ptr<const Encapsulation<time>> other) override
    {
      UNUSED(other);
      assert(deactivated_function::value);
    }

    virtual Matrix<time> as_matrix() const = 0;
    virtual Matrix<time> as_vector() const = 0;
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override final
    {
      UNUSED(a); UNUSED(x);
      throw NotImplementedYet("ax+y not viable for a Particle.");
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, Matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override final
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      throw NotImplementedYet("aA*x not viable for a Particle.");
    }
    //! @}

    //! @{
          scalar&            mass()         { return this->m_mass; }
    const scalar&            mass() const   { return this->m_mass; }
          scalar&            charge()       { return this->m_charge; }
    const scalar&            charge() const { return this->m_charge; }
          scalar&            alpha()        { return this->m_alpha; }
    const scalar&            alpha() const  { return this->m_alpha; }
          position_type&     pos()          { return this->m_pos; }
    const position_type&     pos() const    { return this->m_pos; }
          velocity_type&     vel()          { return this->m_vel; }
    const velocity_type&     vel() const    { return this->m_vel; }
          acceleration_type& accel()        { return this->m_accel; }
    const acceleration_type& accel() const  { return this->m_accel; }
    //! @}
};


#endif  // _EXAMPLES__BORIS__PARTICLE__HPP_
