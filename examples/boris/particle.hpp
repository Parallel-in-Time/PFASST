#ifndef _EXAMPLES__BORIS__PARTICLE__HPP_
#define _EXAMPLES__BORIS__PARTICLE__HPP_

#include <cstdlib>
#include <stdexcept>
#include <memory>
#include <cassert>

#include <boost/numeric/ublas/matrix.hpp>

#include <pfasst/globals.hpp>
#include <pfasst/interfaces.hpp>
#include <pfasst/encap/encapsulation.hpp>

using namespace std;
using namespace pfasst;
using namespace pfasst::encap;


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
      throw NotImplementedYet("copy of Encapsulation into ParticleComponent");
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      UNUSED(a); UNUSED(x);
      throw NotImplementedYet("ax+y for Encapsulation onto ParticleComponent");
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      throw NotImplementedYet("aA*x for Encapsulation onto ParticleComponent");
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
    AccelerationEncapsulation()
    {}

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
    //! @}

    //! @{
    scalar mass;
    scalar charge;
    shared_ptr<position_type> pos;
    shared_ptr<velocity_type> vel;
    shared_ptr<acceleration_type> accel;
    //! @}

  public:
    //! @{
    ParticleEncapsulation()
      : ParticleEncapsulation(scalar(1.0), scalar(1.0))
    {}

    ParticleEncapsulation(scalar mass, scalar charge)
      :   mass(mass)
        , charge(charge)
        , pos(make_shared<position_type>())
        , vel(make_shared<velocity_type>())
        , accel(make_shared<acceleration_type>())
    {}

    virtual ~ParticleEncapsulation()
    {}
    //! @}

    //! @{
    virtual void copy(shared_ptr<const Encapsulation<time>> other) override
    {
      shared_ptr<const ParticleEncapsulation<scalar, time>> other_cast = \
        dynamic_pointer_cast<const ParticleEncapsulation<scalar, time>>(other);
      assert(other_cast);
      this->copy(other_cast);
    }

    virtual void copy(shared_ptr<const ParticleEncapsulation<scalar, time>> other)
    {
      assert(other->DIM == this->DIM);
      this->mass = other->mass;
      this->charge = other->charge;
      this->pos->copy(other->const_pos());
      this->vel->copy(other->const_vel());
      this->accel->copy(other->const_accel());
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override final
    {
      UNUSED(a); UNUSED(x);
      throw NotImplementedYet("ax+y not viable for a Particle.");
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override final
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      throw NotImplementedYet("aA*x not viable for a Particle.");
    }
    //! @}

    //! @{
    virtual shared_ptr<const position_type> const_pos() const
    {
      return const_pointer_cast<const position_type>(this->pos);
    }
    virtual shared_ptr<const velocity_type> const_vel() const
    {
      return const_pointer_cast<const velocity_type>(this->vel);
    }
    virtual shared_ptr<const acceleration_type> const_accel() const
    {
      return const_pointer_cast<const acceleration_type>(this->accel);
    }
    //! @}
};


#endif  // _EXAMPLES__BORIS__PARTICLE__HPP_
