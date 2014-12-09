#ifndef _EXAMPLES__BORIS__PARTICLE_HPP_
#define _EXAMPLES__BORIS__PARTICLE_HPP_

#ifndef LOG_PRECISION
  #define LOG_PRECISION 5
#endif


#include <memory>
#include <vector>
using namespace std;

// #include <pfasst/easylogging++.h>


namespace examples
{
  namespace boris
  {
    template<
      typename precision
    >
    using ParticleComponent = vector<precision>;


    template<
      typename precision
    >
    class Particle
//       : public el::Loggable
    {
      protected:
        size_t _dim;
        precision _charge;
        precision _mass;
        ParticleComponent<precision> _pos;
        ParticleComponent<precision> _vel;

      public:
        explicit Particle(const size_t dim = 3);
        Particle(const size_t dim, const precision charge, const precision mass);
        virtual ~Particle();

              size_t DIM() const;
              ParticleComponent<precision>& pos();
        const ParticleComponent<precision>& pos() const;
              ParticleComponent<precision>& vel();
        const ParticleComponent<precision>& vel() const;
        const precision charge() const;
        const precision mass() const;
        const precision alpha() const;

//         void set_pos(const ParticleComponent<precision>& pos);
//         void set_vel(const ParticleComponent<precision>& vel);
        void set_charge(const precision& charge);
        void set_mass(const precision& mass);

//         virtual void log(el::base::type::ostream_t& os) const;
    };
  }  // ::examples::boris
}  // ::examples


#include "particle_impl.hpp"

#endif  // _EXAMPLES__BORIS__PARTICLE_HPP_
