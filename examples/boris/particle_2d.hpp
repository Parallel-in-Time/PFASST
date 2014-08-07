#ifndef _EXAMPLES__BORIS__PARTICLE_2D__HPP_
#define _EXAMPLES__BORIS__PARTICLE_2D__HPP_

#include <memory>
#include <stdexcept>
#include <cassert>
#include <type_traits>

#include "particle.hpp"

using namespace pfasst;
using namespace pfasst::encap;


template<
  typename scalar,
  typename time = time_precision
>
class Position2DEncapsulation
  : public PositionEncapsulation<scalar, time>
{
  public:
    const size_t DIM = 2;
    scalar x;
    scalar y;

    //! @{
    Position2DEncapsulation()
      : Position2DEncapsulation(scalar(0.0), scalar(0.0))
    {}

    Position2DEncapsulation(const scalar x, const scalar y)
      :   x(x)
        , y(y)
    {}

    virtual ~Position2DEncapsulation()
    {}
    //! @}

    //! @{
    virtual void zero() override
    {
      this->x = scalar(0.0);
      this->y = scalar(0.0);
    }

    virtual void copy(shared_ptr<const Encapsulation<time>> other) override
    {
      shared_ptr<const Position2DEncapsulation<scalar, time>> other_cast = \
        dynamic_pointer_cast<const Position2DEncapsulation<scalar, time>>(other);
      assert(other_cast);
      this->copy(other_cast);
    }

    virtual void copy(shared_ptr<const Position2DEncapsulation<scalar, time>> other)
    {
      assert(this->DIM == other->DIM);
      this->x = other->x;
      this->y = other->y;
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      shared_ptr<const Position2DEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const Position2DEncapsulation<scalar, time>>(x);
      assert(x_cast);
      this->saxpy(a, x_cast);
    }

    virtual void saxpy(time a, shared_ptr<const Position2DEncapsulation<scalar, time>> x)
    {
      this->x += a * x->x;
      this->y += a * x->y;
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<Position2DEncapsulation<scalar, time>>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<Position2DEncapsulation<scalar, time>>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<Position2DEncapsulation<scalar, time>>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    virtual void mat_apply(vector<shared_ptr<Position2DEncapsulation<scalar, time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Position2DEncapsulation<scalar, time>>> src, 
                           bool zero = true)
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      // TODO: implement aA*x for 2D-Position
    }
    //! @}
};


template<
  typename scalar,
  typename time = time_precision
>
class Velocity2DEncapsulation
  : public VelocityEncapsulation<scalar, time>
{
  public:
    const size_t DIM = 2;
    scalar u;
    scalar v;

    //! @{
    Velocity2DEncapsulation()
      : Velocity2DEncapsulation(scalar(0.0), scalar(0.0))
    {}

    Velocity2DEncapsulation(const scalar u, const scalar v)
      :   u(u)
        , v(v)
    {}

    virtual ~Velocity2DEncapsulation()
    {}
    //! @}

    //! @{
    virtual void zero() override
    {
      this->u = scalar(0.0);
      this->v = scalar(0.0);
    }

    virtual void copy(shared_ptr<const Encapsulation<time>> other) override
    {
      shared_ptr<const Velocity2DEncapsulation<scalar, time>> other_cast = \
        dynamic_pointer_cast<const Velocity2DEncapsulation<scalar, time>>(other);
      assert(other_cast);
      this->copy(other_cast);
    }

    virtual void copy(shared_ptr<const Velocity2DEncapsulation<scalar, time>> other) 
    {
      assert(this->DIM == other->DIM);
      this->u = other->u;
      this->v = other->v;
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      shared_ptr<const Velocity2DEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const Velocity2DEncapsulation<scalar, time>>(x);
      assert(x_cast);
      this->saxpy(a, x_cast);
    }

    virtual void saxpy(time a, shared_ptr<const Velocity2DEncapsulation<scalar, time>> x)
    {
      this->u += a * x->u;
      this->v += a * x->v;
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<Velocity2DEncapsulation<scalar, time>>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<Velocity2DEncapsulation<scalar, time>>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<Velocity2DEncapsulation<scalar, time>>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    virtual void mat_apply(vector<shared_ptr<Velocity2DEncapsulation<scalar, time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Velocity2DEncapsulation<scalar, time>>> src, 
                           bool zero = true)
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      // TODO: implement aA*x for 2D-Velocity
    }
    //! @}
};


template<
  typename scalar,
  typename time = time_precision
>
class Acceleration2DEncapsulation
  : public AccelerationEncapsulation<scalar, time>
{
  public:
    const size_t DIM = 2;
    scalar a;
    scalar b;

    //! @{
    Acceleration2DEncapsulation()
      : Acceleration2DEncapsulation(scalar(0.0), scalar(0.0))
    {}

    Acceleration2DEncapsulation(const scalar a, const scalar b)
      :   a(a)
        , b(b)
    {}

    virtual ~Acceleration2DEncapsulation()
    {}
    //! @}

    //! @{
    virtual void zero() override
    {
      this->a = scalar(0.0);
      this->b = scalar(0.0);
    }

    virtual void copy(shared_ptr<const Encapsulation<time>> other) override
    {
      shared_ptr<const Acceleration2DEncapsulation<scalar, time>> other_cast = \
        dynamic_pointer_cast<const Acceleration2DEncapsulation<scalar, time>>(other);
      assert(other_cast);
      this->copy(other_cast);
    }

    virtual void copy(shared_ptr<const Acceleration2DEncapsulation<scalar, time>> other) 
    {
      assert(this->DIM == other->DIM);
      this->a = other->a;
      this->b = other->b;
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      shared_ptr<const Acceleration2DEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const Acceleration2DEncapsulation<scalar, time>>(x);
      assert(x_cast);
      this->saxpy(a, x_cast);
    }

    virtual void saxpy(time a, shared_ptr<const Acceleration2DEncapsulation<scalar, time>> x)
    {
      this->a += a * x->a;
      this->b += a * x->b;
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<Acceleration2DEncapsulation<scalar, time>>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<Acceleration2DEncapsulation<scalar, time>>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<Acceleration2DEncapsulation<scalar, time>>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    virtual void mat_apply(vector<shared_ptr<Acceleration2DEncapsulation<scalar, time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Acceleration2DEncapsulation<scalar, time>>> src, 
                           bool zero = true)
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      // TODO: implement aA*x for 2D-Acceleration
    }
    //! @}
};


template<
  typename scalar,
  typename time = time_precision
>
class Particle2DEncapsulation
  : public ParticleEncapsulation<scalar,
                                 time,
                                 Position2DEncapsulation,
                                 Velocity2DEncapsulation,
                                 Acceleration2DEncapsulation>
{
  private:
    typedef ParticleEncapsulation<scalar,
                                  time,
                                  Position2DEncapsulation,
                                  Velocity2DEncapsulation,
                                  Acceleration2DEncapsulation> parent_type;

  public:
    //! @{
    const size_t DIM = 2;
    //! @}

  public:
    //! @{
    Particle2DEncapsulation()
      : Particle2DEncapsulation(scalar(1.0), scalar(1.0))
    {}

    Particle2DEncapsulation(scalar mass, scalar charge)
      : parent_type(mass, charge)
    {}

    virtual ~Particle2DEncapsulation()
    {}
    //! @}

    //! @{
    virtual void copy(shared_ptr<const Encapsulation<time>> other) override
    {
      shared_ptr<const Particle2DEncapsulation<scalar, time>> other_cast = \
        dynamic_pointer_cast<const Particle2DEncapsulation<scalar, time>>(other);
      assert(other_cast);
      this->copy(other_cast);
    }

    virtual void copy(shared_ptr<const ParticleEncapsulation<scalar, time>> other) override
    {
      shared_ptr<const Particle2DEncapsulation<scalar, time>> other_cast = \
        dynamic_pointer_cast<const Particle2DEncapsulation<scalar, time>>(other);
      assert(other_cast);
      this->copy(other_cast);
    }

    virtual void copy(shared_ptr<const Particle2DEncapsulation<scalar, time>> other)
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
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x)
    {
      shared_ptr<const Particle2DEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const Particle2DEncapsulation<scalar, time>>(x);
      assert(x_cast);
      this->saxpy(a, x_cast);
    }
    virtual void saxpy(time a, shared_ptr<const Particle2DEncapsulation<scalar, time>> x)
    {
      UNUSED(a); UNUSED(x);
      // TODO: implement ax+y for 2D-Particle
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true)
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<Particle2DEncapsulation<scalar, time>>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<Particle2DEncapsulation<scalar, time>>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<Particle2DEncapsulation<scalar, time>>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    virtual void mat_apply(vector<shared_ptr<Particle2DEncapsulation<scalar, time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Particle2DEncapsulation<scalar, time>>> src, 
                           bool zero = true)
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      // TODO: implement aA*x for 2D-Particle
    }

    /**
     * gives the particle's energy as its maximum norm
     */
    scalar norm0() const
    {
      // TODO: implement maximum norm for 2D-Particle
    }
    //! @}
};

// XXX: not sure whether this factory is actually useful
template<
  typename scalar,
  typename time = time_precision
>
class Particle2DFactory
  : public EncapFactory<time>
{
  private:
    scalar mass;
    scalar charge;
  public:
    int dofs() { return 1; }
    Particle2DFactory(scalar mass, scalar charge)
      : mass(mass), charge(charge)
    {}
    virtual shared_ptr<Encapsulation<time>> create(const EncapType)
    {
      return make_shared<Particle2DEncapsulation<scalar, time>>(this->mass, this->charge);
    }
};

#endif  // _EXAMPLES__BORIS__PARTICLE_2D__HPP_
