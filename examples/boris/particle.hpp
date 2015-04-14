/**
 * @file examples/boris/particle.hpp
 * @ingroup BorisFiles
 */
#ifndef _EXAMPLES__BORIS__PARTICLE_HPP_
#define _EXAMPLES__BORIS__PARTICLE_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/logging.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      /**
       * @ingroup Boris
       */
      template<
        typename precision
      >
      using ParticleComponent = vector<precision>;

      /**
       * @ingroup Boris
       */
      template<typename T>
      inline el::base::type::ostream_t& operator<<(el::base::type::ostream_t& os, const vector<T>& vec);


      /**
       * @ingroup Boris
       */
      template<
        typename precision
      >
      class Particle
        : public el::Loggable
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

          inline size_t dim() const;
                ParticleComponent<precision>& pos();
          const ParticleComponent<precision>& pos() const;
                ParticleComponent<precision>& vel();
          const ParticleComponent<precision>& vel() const;
          const precision charge() const;
          const precision mass() const;
          const precision alpha() const;

          void set_charge(const precision& charge);
          void set_mass(const precision& mass);

          virtual void log(el::base::type::ostream_t& os) const;
      };


      /**
       * @ingroup Boris
       */
      template<typename precision>
      inline el::base::type::ostream_t& operator<<(el::base::type::ostream_t& os,
                                                   const shared_ptr<Particle<precision>>& sp_particle);
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

// #include "particle_util.hpp"
#include "particle_impl.hpp"

#endif  // _EXAMPLES__BORIS__PARTICLE_HPP_
