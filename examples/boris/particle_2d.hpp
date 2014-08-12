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
  private:
    typedef Position2DEncapsulation<scalar, time> this_type;

  public:
    const size_t DIM = 2;
    scalar x;
    scalar y;

    //! @{
    //! Default CTor
    Position2DEncapsulation()
      : Position2DEncapsulation(scalar(0.0), scalar(0.0))
    {}

    Position2DEncapsulation(const scalar x, const scalar y)
      :   x(x)
        , y(y)
    {}

    //! Copy CTor
    Position2DEncapsulation(const this_type& other)
      :   x(other.x)
        , y(other.y)
    {}

    //! Move CTor
    Position2DEncapsulation(this_type&& other)
      :   x(std::move(other.x))
        , y(std::move(other.y))
    {}

    //! Copy CTor for shared_ptr
    Position2DEncapsulation(const shared_ptr<const this_type> other)
      :   x(other->x)
        , y(other->y)
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
      *this += a * x;
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

    //! @{
    this_type& operator=(const this_type& other)
    {
      if (this != &other) {
        assert(this->DIM == other.DIM);
        this->x = other.x;
        this->y = other.y;
      }
      return *this;
    }

    this_type& operator=(this_type&& other)
    {
      assert(this != &other);
      assert(this->DIM == other.DIM);
      this->x = other.x;
      this->y = other.y;
      return *this;
    }

    this_type& operator=(shared_ptr<const this_type> other)
    {
      if (this != other.get()) {
        assert(this->DIM == other->DIM);
        this->x = other->x;
        this->y = other->y;
      }
      return *this;
    }

    template<typename value_type>
    this_type& operator*=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      this->x *= value;
      this->y *= value;
      return *this;
    }

    template<typename value_type>
    friend this_type operator*(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      return this_type(first.x * value, first.y * value);
    }

    template<typename value_type>
    friend this_type operator*(const this_type& first, const value_type value)
    {
      return value * first;
    }

    template<typename value_type>
    friend this_type operator*(const value_type value, const shared_ptr<const this_type> first)
    {
      return value * *(first.get());
    }

    template<typename value_type>
    friend this_type operator*(const shared_ptr<const this_type> first, const value_type value)
    {
      return value * *(first.get());
    }

    this_type& operator+=(const shared_ptr<const this_type> second)
    {
      return this += second.get();
    }

    this_type& operator+=(const this_type& second)
    {
      assert(this->DIM == second.DIM);
      this->x += second.x;
      this->y += second.y;
      return *this;
    }

    template<typename value_type>
    this_type& operator+=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Position.");
      this->x += value;
      this->y += value;
      return *this;
    }

    friend this_type operator+(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.x + second.x, first.y + second.y);
    }

    friend this_type operator+(const shared_ptr<const this_type> first,
                               const shared_ptr<const this_type> second)
    {
      return first.get() + second.get();
    }

    friend this_type operator+(const this_type& first, const shared_ptr<const this_type> second)
    {
      return first + second.get();
    }

    friend this_type operator+(const shared_ptr<const this_type> first, const this_type& second)
    {
      return second + first.get();
    }

    template<typename value_type>
    friend this_type operator+(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Position.");
      return this_type(first.x + value, first.y + value);
    }

    template<typename value_type>
    friend this_type operator+(const this_type& first, const value_type value)
    {
      return value + first;
    }

    template<typename value_type>
    friend this_type operator+(const value_type value, const shared_ptr<const this_type> first)
    {
      return value + first.get();
    }

    template<typename value_type>
    friend this_type operator+(const shared_ptr<const this_type> first, const value_type value)
    {
      return value + first.get();
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
  private:
    typedef Velocity2DEncapsulation<scalar, time> this_type;

  public:
    const size_t DIM = 2;
    scalar u;
    scalar v;

    //! @{
    //! Default CTor
    Velocity2DEncapsulation()
      : Velocity2DEncapsulation(scalar(0.0), scalar(0.0))
    {}

    Velocity2DEncapsulation(const scalar u, const scalar v)
      :   u(u)
        , v(v)
    {}

    //! Copy CTor
    Velocity2DEncapsulation(const this_type& other)
      :   u(other.u)
        , v(other.v)
    {}

    //! Move CTor
    Velocity2DEncapsulation(this_type&& other)
      :   u(std::move(other.u))
        , v(std::move(other.v))
    {}

    //! Copy CTor for shared_ptr
    Velocity2DEncapsulation(const shared_ptr<const this_type> other)
      :   u(other->v)
        , v(other->u)
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

    //! @{
    template<typename dt_precision>
    Position2DEncapsulation<scalar, time> convert(const dt<dt_precision> dt)
    {
      return Position2DEncapsulation<scalar, time>(this->u * dt.v, this->v * dt.v);
    }
    //! @}

    //! @{
    this_type& operator=(const this_type& other)
    {
      if (this != &other) {
        assert(this->DIM == other.DIM);
        this->u = other.u;
        this->v = other.v;
      }
      return *this;
    }

    this_type& operator=(this_type&& other)
    {
      assert(this != &other);
      assert(this->DIM == other.DIM);
      this->u = other.u;
      this->v = other.v;
      return *this;
    }

    this_type& operator=(shared_ptr<const this_type> other)
    {
      if (this != other.get()) {
        assert(this->DIM == other->DIM);
        this->u = other->u;
        this->v = other->v;
      }
      return *this;
    }

    template<typename value_type>
    this_type& operator*=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Velocity can only be multiplied by arithmetic types.");
      this->u *= value;
      this->v *= value;
      return *this;
    }

    template<typename value_type>
    friend this_type operator*(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Velocity can only be multiplied by arithmetic types.");
      return this_type(first.u * value, first.v * value);
    }

    template<typename value_type>
    friend this_type operator*(const this_type& first, const value_type value)
    {
      return value * first;
    }

    template<typename value_type>
    friend this_type operator*(const value_type value, const shared_ptr<const this_type> first)
    {
      return value * *(first.get());
    }

    template<typename value_type>
    friend this_type operator*(const shared_ptr<const this_type> first, const value_type value)
    {
      return value * *(first.get());
    }

    this_type& operator+=(const shared_ptr<const this_type> second)
    {
      return this += second.get();
    }

    this_type& operator+=(const this_type& second)
    {
      assert(this->DIM == second.DIM);
      this->u += second.u;
      this->v += second.v;
      return *this;
    }

    template<typename value_type>
    this_type& operator+=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Velocity.");
      this->u += value;
      this->v += value;
      return *this;
    }

    friend this_type operator+(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.u + second.u, first.v + second.v);
    }

    friend this_type operator+(const shared_ptr<const this_type> first,
                               const shared_ptr<const this_type> second)
    {
      return first.get() + second.get();
    }

    friend this_type operator+(const this_type& first, const shared_ptr<const this_type> second)
    {
      return first + second.get();
    }

    friend this_type operator+(const shared_ptr<const this_type> first, const this_type& second)
    {
      return second + first.get();
    }

    template<typename value_type>
    friend this_type operator+(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Velocity.");
      return this_type(first.u + value, first.v + value);
    }

    template<typename value_type>
    friend this_type operator+(const this_type& first, const value_type value)
    {
      return value + first;
    }

    template<typename value_type>
    friend this_type operator+(const value_type value, const shared_ptr<const this_type> first)
    {
      return value + first.get();
    }

    template<typename value_type>
    friend this_type operator+(const shared_ptr<const this_type> first, const value_type value)
    {
      return value + first.get();
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
  private:
    typedef Acceleration2DEncapsulation<scalar, time> this_type;

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

    //! Copy CTor
    Acceleration2DEncapsulation(const this_type& other)
      :   a(other.a)
        , b(other.b)
    {}

    //! Move CTor
    Acceleration2DEncapsulation(this_type&& other)
      :   a(std::move(other.a))
        , b(std::move(other.b))
    {}

    //! Copy CTor for shared_ptr
    Acceleration2DEncapsulation(const shared_ptr<const this_type> other)
      :   a(other->a)
        , b(other->b)
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
      *this += a * x;
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

    //! @{
    template<typename dtdt_precision>
    Position2DEncapsulation<scalar, time> convert(const dtdt<dtdt_precision> dtdt)
    {
      return Position2DEncapsulation<scalar, time>(this->a * dtdt.v, this->b * dtdt.v);
    }

    template<typename dt_precision>
    Velocity2DEncapsulation<scalar, time> convert(const dt<dt_precision> dt)
    {
      return Velocity2DEncapsulation<scalar, time>(this->a * dt.v, this->b * dt.v);
    }
    //! @}

    //! @{
    this_type& operator=(const this_type& other)
    {
      if (this != &other) {
        assert(this->DIM == other.DIM);
        this->a = other.a;
        this->b = other.b;
      }
      return *this;
    }

    this_type& operator=(this_type&& other)
    {
      assert(this != &other);
      assert(this->DIM == other.DIM);
      this->a = other.a;
      this->b = other.b;
      return *this;
    }

    this_type& operator=(shared_ptr<const this_type> other)
    {
      if (this != other.get()) {
        assert(this->DIM == other->DIM);
        this->a = other->a;
        this->b = other->b;
      }
      return *this;
    }

    template<typename value_type>
    this_type& operator*=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Acceleration can only be multiplied by arithmetic types.");
      this->a *= value;
      this->b *= value;
      return *this;
    }

    template<typename value_type>
    friend this_type operator*(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Acceleration can only be multiplied by arithmetic types.");
      return this_type(first.a * value, first.b * value);
    }

    template<typename value_type>
    friend this_type operator*(const this_type& first, const value_type value)
    {
      return value * first;
    }

    template<typename value_type>
    friend this_type operator*(const value_type value, const shared_ptr<const this_type> first)
    {
      return value * *(first.get());
    }

    template<typename value_type>
    friend this_type operator*(const shared_ptr<const this_type> first, const value_type value)
    {
      return value * *(first.get());
    }

    this_type& operator+=(const shared_ptr<const this_type> second)
    {
      return this += second.get();
    }

    this_type& operator+=(const this_type& second)
    {
      assert(this->DIM == second.DIM);
      this->a += second.a;
      this->b += second.b;
      return *this;
    }

    template<typename value_type>
    this_type& operator+=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Acceleration.");
      this->a += value;
      this->b += value;
      return *this;
    }

    friend this_type operator+(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.a + second.a, first.b + second.b);
    }

    friend this_type operator+(const shared_ptr<const this_type> first,
                               const shared_ptr<const this_type> second)
    {
      return first.get() + second.get();
    }

    friend this_type operator+(const this_type& first, const shared_ptr<const this_type> second)
    {
      return first + second.get();
    }

    friend this_type operator+(const shared_ptr<const this_type> first, const this_type& second)
    {
      return second + first.get();
    }

    template<typename value_type>
    friend this_type operator+(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Acceleration.");
      return this_type(first.a + value, first.b + value);
    }

    template<typename value_type>
    friend this_type operator+(const this_type& first, const value_type value)
    {
      return value + first;
    }

    template<typename value_type>
    friend this_type operator+(const value_type value, const shared_ptr<const this_type> first)
    {
      return value + first.get();
    }

    template<typename value_type>
    friend this_type operator+(const shared_ptr<const this_type> first, const value_type value)
    {
      return value + first.get();
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
    /**
     * gives the particle's energy as its maximum norm
     */
    virtual scalar norm0() const
    {
      // TODO: implement maximum norm for 2D-Particle
      throw NotImplementedYet("how to compute a norm of a Particle?");
      return scalar(-1.0);
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
