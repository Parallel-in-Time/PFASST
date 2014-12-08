#include "particle.hpp"
#include "particle_util.hpp"


#include <cassert>
#include <iomanip>
using namespace std;

namespace examples
{
  namespace boris
  {
    template<typename precision>
    Particle<precision>::Particle(const size_t dim)
      : Particle(dim, precision(1.0), precision(1.0))
    {}

    template<typename precision>
    Particle<precision>::Particle(const size_t dim, const precision charge, const precision mass)
      :   _dim(dim)
        , _charge(charge)
        , _mass(mass)
        , _pos(dim)
        , _vel(dim)
    {
      assert(this->_pos.size() == dim);
      assert(this->_vel.size() == dim);
    }

    template<typename precision>
    Particle<precision>::~Particle()
    {}

    template<typename precision>
    size_t Particle<precision>::DIM() const
    {
      return this->_dim;
    }

    template<typename precision>
    const ParticleComponent<precision>& Particle<precision>::pos() const
    {
      return this->_pos;
    }
    template<typename precision>
    ParticleComponent<precision>& Particle<precision>::pos()
    {
      return this->_pos;
    }

    template<typename precision>
    const ParticleComponent<precision>& Particle<precision>::vel() const
    {
      return this->_vel;
    }
    template<typename precision>
    ParticleComponent<precision>& Particle<precision>::vel()
    {
      return this->_vel;
    }

    template<typename precision>
    const precision Particle<precision>::charge() const
    {
      return this->_charge;
    }

    template<typename precision>
    const precision Particle<precision>::mass() const
    {
      return this->_mass;
    }

    template<typename precision>
    const precision Particle<precision>::alpha() const
    {
      return this->_charge / this->_mass;
    }

    template<typename precision>
    void Particle<precision>::set_charge(const precision& charge)
    {
      this->_charge = charge;
    }

    template<typename precision>
    void Particle<precision>::set_mass(const precision& mass)
    {
      this->_mass = mass;
    }

    template<typename precision>
    void Particle<precision>::log(el::base::type::ostream_t& os) const
    {
      os << fixed << setprecision(LOG_PRECISION);
//       os << "Particle(q=" << this->_charge << ", m=" << this->_mass << ", pos=" << this->_pos << ", vel=" << this->_vel << ")";
      os.unsetf(ios_base::floatfield);
    }
  }  // ::examples::boris
}  // ::examples
