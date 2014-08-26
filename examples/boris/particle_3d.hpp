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
class Position3DEncapsulation
  : public PositionEncapsulation<scalar, time>
{
  private:
    typedef Position3DEncapsulation<scalar, time> this_type;

  public:
    const size_t DIM = 3;
    scalar x;
    scalar y;
    scalar z;

    //! @{
    //! Default CTor
    Position3DEncapsulation()
      : Position3DEncapsulation(scalar(0.0), scalar(0.0), scalar(0.0))
    {}

    Position3DEncapsulation(const scalar x, const scalar y, const scalar z)
      :   x(x)
        , y(y)
        , z(z)
    {}

    //! Copy CTor
    Position3DEncapsulation(const this_type& other)
      :   x(other.x)
        , y(other.y)
        , z(other.z)
    {}

    //! Move CTor
    Position3DEncapsulation(this_type&& other)
      :   x(std::move(other.x))
        , y(std::move(other.y))
        , z(std::move(other.z))
    {}

    virtual ~Position3DEncapsulation()
    {}
    //! @}

    //! @{
    virtual void zero() override
    {
      this->x = scalar(0.0);
      this->y = scalar(0.0);
      this->z = scalar(0.0);
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      shared_ptr<const this_type> x_cast = dynamic_pointer_cast<const this_type>(x);
      assert(x_cast);
      this->saxpy(a, x_cast);
    }

    virtual void saxpy(time a, shared_ptr<const this_type> x)
    {
      assert(this->DIM == x->DIM);
      this->saxpy(a, *(x.get()));
    }

    virtual void saxpy(time a, const this_type& x)
    {
      assert(this->DIM == x.DIM);
      this->operator+=(a * x);
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<this_type>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<this_type>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<this_type>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    virtual void mat_apply(vector<shared_ptr<this_type>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<this_type>> src, 
                           bool zero = true)
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      // TODO: implement aA*x for 3D-Position
      throw NotImplementedYet("aA*x for 3D-Position");
    }
    //! @}

    //! @{
    this_type& operator=(const this_type& other)
    {
      assert(this != &other);
      assert(this->DIM == other.DIM);
      this_type tmp(other);
      std::swap(this->x, tmp.x);
      std::swap(this->y, tmp.y);
      std::swap(this->z, tmp.z);
      return *this;
    }

    this_type& operator=(this_type&& other)
    {
      assert(this != &other);
      assert(this->DIM == other.DIM);
      this->x = std::move(other.x);
      this->y = std::move(other.y);
      this->z = std::move(other.z);
      return *this;
    }

    //! this *= value
    template<typename value_type>
    this_type& operator*(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      this->x *= value;
      this->y *= value;
      this->z *= value;
      return *this;
    }

    //! newpos = value * pos
    template<typename value_type>
    friend this_type operator*(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      return this_type(value * first.x, value * first.y, value * first.z);
    }

    //! newpos = pos * value
    template<typename value_type>
    friend this_type operator*(const this_type& first, const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      return this_type(first.x * value, first.y * value, first.z * value);
    }

    //! this += value
    template<typename value_type>
    this_type& operator+(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      this->x += value;
      this->y += value;
      this->z += value;
      return *this;
    }

    //! newpos = value + pos
    template<typename value_type>
    friend this_type operator+(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      return this_type(first.x + value, first.y + value, first.z + value);
    }

    //! newpos = pos + value
    template<typename value_type>
    friend this_type operator+(const this_type& first, const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Position can only be multiplied by arithmetic types.");
      return this_type(first.x + value, first.y + value, first.z + value);
    }

    //! newpos = pos + other
    friend this_type operator+(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.x + second.x, first.y + second.y, first.z + second.z);
    }

    //! this += other
    this_type& operator+=(const this_type& other)
    {
      assert(this->DIM == other.DIM);
      this->x += other.x;
      this->y += other.y;
      this->z += other.z;
      return *this;
    }

    //! newpos = pos - other
    friend this_type operator-=(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.x - second.x, first.y - second.y, first.z - second.y);
    }

    //! this -= other
    this_type& operator-(const this_type& other)
    {
      assert(this->DIM == other.DIM);
      this->x -= other.x;
      this->y -= other.y;
      this->z -= other.z;
      return *this;
    }
    //! @}
};


template<
  typename scalar,
  typename time = time_precision
>
class Velocity3DEncapsulation
  : public VelocityEncapsulation<scalar, time>
{
  private:
    typedef Velocity3DEncapsulation<scalar, time> this_type;

  public:
    const size_t DIM = 3;
    scalar u;
    scalar v;
    scalar w;

    //! @{
    //! Default CTor
    Velocity3DEncapsulation()
      : Velocity3DEncapsulation(scalar(0.0), scalar(0.0), scalar(0.0))
    {}

    Velocity3DEncapsulation(const scalar u, const scalar v, const scalar w)
      :   u(u)
        , v(v)
        , w(w)
    {}

    //! Copy CTor
    Velocity3DEncapsulation(const this_type& other)
      :   u(other.u)
        , v(other.v)
        , w(other.w)
    {}

    //! Move CTor
    Velocity3DEncapsulation(this_type&& other)
      :   u(std::move(other.u))
        , v(std::move(other.v))
        , w(std::move(other.w))
    {}

    //! Copy CTor for shared_ptr
    Velocity3DEncapsulation(const shared_ptr<const this_type> other)
      :   u(other->v)
        , v(other->u)
        , w(other->w)
    {}

    virtual ~Velocity3DEncapsulation()
    {}
    //! @}

    //! @{
    virtual void zero() override
    {
      this->u = scalar(0.0);
      this->v = scalar(0.0);
      this->w = scalar(0.0);
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      shared_ptr<const this_type> x_cast = dynamic_pointer_cast<const this_type>(x);
      assert(x_cast);
      this->saxpy(a, x_cast);
    }

    virtual void saxpy(time a, shared_ptr<const this_type> x)
    {
      assert(this->DIM == x->DIM);
      this->saxpy(a, *(x.get()));
    }

    virtual void saxpy(time a, const this_type& x)
    {
      assert(this->DIM == x.DIM);
      this->operator+=(a * x);
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<this_type>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<this_type>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<this_type>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    virtual void mat_apply(vector<shared_ptr<this_type>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<this_type>> src, 
                           bool zero = true)
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      // TODO: implement aA*x for 3D-Velocity
      throw NotImplementedYet("aA*x for 3D-Velocity");
    }
    //! @}

    //! @{
    template<typename dt_precision>
    Position3DEncapsulation<scalar, time> convert(const dt<dt_precision> dt) const
    {
      return Position3DEncapsulation<scalar, time>(this->u * dt.v, this->v * dt.v, this->w * dt.v);
    }
    //! @}

    //! @{
    this_type& operator=(const this_type& other)
    {
      if (this != &other) {
        assert(this->DIM == other.DIM);
        this->u = other.u;
        this->v = other.v;
        this->w = other.w;
      }
      return *this;
    }

    this_type& operator=(this_type&& other)
    {
      assert(this != &other);
      assert(this->DIM == other.DIM);
      this->u = other.u;
      this->v = other.v;
      this->w = other.w;
      return *this;
    }

    this_type& operator=(const shared_ptr<const this_type> other)
    {
      if (this != other.get()) {
        assert(this->DIM == other->DIM);
        this->u = other->u;
        this->v = other->v;
        this->w = other->w;
      }
      return *this;
    }

    //! vel += value
    template<typename value_type>
    this_type& operator+=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Velocity.");
      this->u += value;
      this->v += value;
      this->w += value;
      return *this;
    }

    //! vel += other
    this_type& operator+=(const this_type& second)
    {
      assert(this->DIM == second.DIM);
      this->u += second.u;
      this->v += second.v;
      this->w += second.w;
      return *this;
    }

    //! newvel = vel + value
    template<typename value_type>
    friend this_type operator+(const this_type& first, const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Velocity.");
      return this_type(first.u + value, first.v + value, first.w + value);
    }

    //! newvel = value + vel
    template<typename value_type>
    friend this_type operator+(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Velocity.");
      return this_type(first.u + value, first.v + value, first.w + value);
    }

    //! newvel = vel + other
    friend this_type operator+(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.u + second.u, first.v + second.v, first.w + second.w);
    }

    //! vel *= value
    template<typename value_type>
    this_type& operator*=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Velocity can only be multiplied by arithmetic types.");
      this->u *= value;
      this->v *= value;
      this->w *= value;
      return *this;
    }

    //! newvel = value * vel
    template<typename value_type>
    friend this_type operator*(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Velocity can only be multiplied by arithmetic types.");
      return this_type(first.u * value, first.v * value, first.w * value);
    }

    //! newvel = vel * value
    template<typename value_type>
    friend this_type operator*(const this_type& first, const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Velocity can only be multiplied by arithmetic types.");
      return this_type(first.u * value, first.v * value, first.w * value);
    }
};


template<typename scalar, typename time>
//! newvel = vel x vel
Velocity3DEncapsulation<scalar, time>
dot(const Velocity3DEncapsulation<scalar, time>& first,
    const Velocity3DEncapsulation<scalar, time>& second)
{
  assert(first.DIM == second.DIM);
  return Velocity3DEncapsulation<scalar, time>(first.v * second.w - first.w * second.v,
                                               first.w * second.u - first.u * second.w,
                                               first.u * second.v - first.v * second.u);
}


template<
  typename scalar,
  typename time = time_precision
>
class Acceleration3DEncapsulation
  : public AccelerationEncapsulation<scalar, time>
{
  private:
    typedef Acceleration3DEncapsulation<scalar, time> this_type;

  public:
    const size_t DIM = 3;
    scalar a;
    scalar b;
    scalar c;

    //! @{
    Acceleration3DEncapsulation()
      : Acceleration3DEncapsulation(scalar(0.0), scalar(0.0), scalar(0.0))
    {}

    Acceleration3DEncapsulation(const scalar a, const scalar b, const scalar c)
      :   a(a)
        , b(b)
        , c(c)
    {}

    //! Copy CTor
    Acceleration3DEncapsulation(const this_type& other)
      :   a(other.a)
        , b(other.b)
        , c(other.c)
    {}

    //! Move CTor
    Acceleration3DEncapsulation(this_type&& other)
      :   a(std::move(other.a))
        , b(std::move(other.b))
        , c(std::move(other.c))
    {}

    //! Copy CTor for shared_ptr
    Acceleration3DEncapsulation(const shared_ptr<const this_type> other)
      :   a(other->a)
        , b(other->b)
        , c(other->c)
    {}

    virtual ~Acceleration3DEncapsulation()
    {}
    //! @}

    //! @{
    virtual void zero() override
    {
      this->a = scalar(0.0);
      this->b = scalar(0.0);
      this->c = scalar(0.0);
    }
    //! @}

    //! @{
    virtual void saxpy(time a, shared_ptr<const Encapsulation<time>> x) override
    {
      shared_ptr<const this_type> x_cast = dynamic_pointer_cast<const this_type>(x);
      assert(x_cast);
      this->saxpy(a, x_cast);
    }

    virtual void saxpy(time a, shared_ptr<const this_type> x)
    {
      assert(this->DIM == x->DIM);
      this->saxpy(a, *(x.get()));
    }

    virtual void saxpy(time a, const this_type& x)
    {
      assert(this->DIM == x.DIM);
      this->operator+=(a * x);
    }

    virtual void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<Encapsulation<time>>> src, 
                           bool zero = true) override
    {
      size_t ndst = dst.size();
      size_t nsrc = src.size();

      vector<shared_ptr<this_type>> dst_cast(ndst), src_cast(nsrc);
      for (size_t n = 0; n < ndst; n++) {
        dst_cast[n] = dynamic_pointer_cast<this_type>(dst[n]);
        assert(dst_cast[n]);
      }
      for (size_t m = 0; m < nsrc; m++) {
        src_cast[m] = dynamic_pointer_cast<this_type>(src[m]);
        assert(src_cast[m]);
      }

      dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
    }

    virtual void mat_apply(vector<shared_ptr<this_type>> dst, 
                           time a, matrix<time> mat,
                           vector<shared_ptr<this_type>> src, 
                           bool zero = true)
    {
      UNUSED(dst); UNUSED(a); UNUSED(mat); UNUSED(src); UNUSED(zero);
      // TODO: implement aA*x for 3D-Acceleration
      throw NotImplementedYet("aA*x for 3D-Acceleration");
    }
    //! @}

    //! @{
    template<typename dtdt_precision>
    Position3DEncapsulation<scalar, time> convert(const dtdt<dtdt_precision> dtdt) const
    {
      return Position3DEncapsulation<scalar, time>(this->a * dtdt.v, this->b * dtdt.v, this->c * dtdt.v);
    }

    template<typename dt_precision>
    Velocity3DEncapsulation<scalar, time> convert(const dt<dt_precision> dt) const
    {
      return Velocity3DEncapsulation<scalar, time>(this->a * dt.v, this->b * dt.v, this->c * dt.v);
    }
    //! @}

    //! @{
    this_type& operator=(const this_type& other)
    {
      if (this != &other) {
        assert(this->DIM == other.DIM);
        this->a = other.a;
        this->b = other.b;
        this->c = other.c;
      }
      return *this;
    }

    this_type& operator=(this_type&& other)
    {
      assert(this != &other);
      assert(this->DIM == other.DIM);
      this->a = other.a;
      this->b = other.b;
      this->c = other.c;
      return *this;
    }

    this_type& operator=(const shared_ptr<const this_type> other)
    {
      if (this != other.get()) {
        assert(this->DIM == other->DIM);
        this->a = other->a;
        this->b = other->b;
        this->c = other->c;
      }
      return *this;
    }

    //! newaccel = value * accel
    template<typename value_type>
    friend this_type operator*(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Acceleration can only be multiplied by arithmetic types.");
      return this_type(first.a * value, first.b * value, first.c * value);
    }

    //! newaccel = accel * value
    template<typename value_type>
    friend this_type operator*(const this_type& first, const value_type value)
    {
      return this_type(value * first.a, value * first.b, value * first.c);
    }

    //! newaccel = accel + value
    template<typename value_type>
    friend this_type operator+(const this_type& first, const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Acceleration can only be multiplied by arithmetic types.");
      return this_type(first.a + value, first.b + value, first.c + value);
    }

    //! newaccel = value + accel
    template<typename value_type>
    friend this_type operator+(const value_type value, const this_type& first)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Acceleration can only be multiplied by arithmetic types.");
      return this_type(value + first.a, value + first.b, value + first.c);
    }

    //! newaccel = accel + other
    friend this_type operator+(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.a + second.a, first.b + second.b, first.c + second.c);
    }

    //! newaccel = accel - other
    friend this_type operator-(const this_type& first, const this_type& second)
    {
      assert(first.DIM == second.DIM);
      return this_type(first.a - second.a, first.b - second.b, first.c - second.c);
    }

    //! accel *= value
    template<typename value_type>
    this_type& operator*=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Acceleration can only be multiplied by arithmetic types.");
      this->a *= value;
      this->b *= value;
      this->c *= value;
      return *this;
    }

    //! accel += value
    template<typename value_type>
    this_type& operator+=(const value_type value)
    {
      static_assert(is_arithmetic<value_type>::value,
                    "Only arithmetic types can be added to Acceleration.");
      this->a += value;
      this->b += value;
      this->c += value;
      return *this;
    }

    //! accel -= value
    template<typename value_type>
    this_type& operator-=(const value_type value)
    {
      this->a -= value;
      this->b -= value;
      this->c -= value;
      return *this;
    }

    //! accel += other
    this_type& operator+=(const this_type& second)
    {
      assert(this->DIM == second.DIM);
      this->a += second.a;
      this->b += second.b;
      this->c += second.c;
      return *this;
    }

    //! accel -= other
    this_type& operator-=(const this_type& second)
    {
      assert(this->DIM == second.DIM);
      this->a -= second.a;
      this->b -= second.b;
      this->c -= second.c;
      return *this;
    }
};



template<
  typename scalar,
  typename time = time_precision
>
class Particle3DEncapsulation
  : public ParticleEncapsulation<scalar,
                                 time,
                                 Position3DEncapsulation,
                                 Velocity3DEncapsulation,
                                 Acceleration3DEncapsulation>
{
  private:
    typedef ParticleEncapsulation<scalar,
                                  time,
                                  Position3DEncapsulation,
                                  Velocity3DEncapsulation,
                                  Acceleration3DEncapsulation> parent_type;
    typedef Particle3DEncapsulation<scalar, time> this_type;

  public:
    //! @{
    const size_t DIM = 3;
    //! @}

  public:
    //! @{
    Particle3DEncapsulation()
      : Particle3DEncapsulation(scalar(1.0), scalar(1.0))
    {}

    Particle3DEncapsulation(const scalar mass, const scalar charge)
      : parent_type(mass, charge)
    {}

    Particle3DEncapsulation(const scalar mass, const scalar charge,
                            const typename parent_type::position_type& pos,
                            const typename parent_type::velocity_type& vel,
                            const typename parent_type::acceleration_type& accel)
      : parent_type(mass, charge, pos, vel, accel)
    {}

    Particle3DEncapsulation(const this_type& other)
      : Particle3DEncapsulation(other.mass(), other.charge(),
                                other.pos(), other.vel(), other.accel())
    {}

    Particle3DEncapsulation(this_type&& other)
      : Particle3DEncapsulation(std::move(other.m_mass), std::move(other.m_charge),
                                std::move(other.m_pos), std::move(other.m_vel),
                                std::move(other.m_accel))
    {}

    Particle3DEncapsulation(shared_ptr<const this_type> other)
      : Particle3DEncapsulation(*(other.get()))
    {}

    virtual ~Particle3DEncapsulation()
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
    /**
     * gives the particle's energy as its maximum norm
     */
    virtual scalar norm0() const
    {
      // TODO: implement maximum norm for 3D-Particle
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
class Particle3DFactory
  : public EncapFactory<time>
{
  private:
    scalar mass;
    scalar charge;
  public:
    int dofs() { return 1; }
    Particle3DFactory(scalar mass, scalar charge)
      : mass(mass), charge(charge)
    {}
    virtual shared_ptr<Encapsulation<time>> create(const EncapType)
    {
      return make_shared<Particle3DEncapsulation<scalar, time>>(this->mass, this->charge);
    }
};

#endif  // _EXAMPLES__BORIS__PARTICLE_2D__HPP_
