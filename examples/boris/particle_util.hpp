#ifndef _EXAMPLES__BORIS__PARTICLE_UTIL_HPP_
#define _EXAMPLES__BORIS__PARTICLE_UTIL_HPP_

#include <memory>
#include <vector>
using namespace std;

#include "particle.hpp"
#include "particle_cloud.hpp"


namespace examples
{
  namespace boris
  {
    template<typename precision>
    static ParticleComponent<precision>      cross_prod(const ParticleComponent<precision>&      first,
                                                        const ParticleComponent<precision>&      second);
    template<typename precision>
    static ParticleCloudComponent<precision> cross_prod(const ParticleCloudComponent<precision>& first,
                                                        const ParticleCloudComponent<precision>& second);
    template<typename precision>
    static ParticleCloudComponent<precision> cross_prod(const ParticleCloudComponent<precision>& first,
                                                        const ParticleComponent<precision>&      second);


    template<typename precision>
    static precision distance(const Particle<precision>& first,
                              const Particle<precision>& second);
    template<typename precision>
    static precision distance(const shared_ptr<Particle<precision>>& first,
                              const shared_ptr<Particle<precision>>& second);

    template<typename precision>
    static vector<precision> distance_to_reference(const ParticleCloud<precision>& cloud,
                                                   const Particle<precision>&      reference);
    template<typename precision>
    static vector<precision> distance_to_reference(const shared_ptr<ParticleCloud<precision>>& cloud,
                                                   const shared_ptr<Particle<precision>>&      reference);


    template<typename precision>
    ParticleComponent<precision>       operator+ (const ParticleComponent<precision>&      first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision>
    ParticleCloudComponent<precision>  operator+ (const ParticleCloudComponent<precision>& first,
                                                  const ParticleCloudComponent<precision>& second);
    template<typename precision>
    ParticleCloudComponent<precision>  operator+ (const ParticleCloudComponent<precision>& first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision, typename ValueT>
    ParticleComponent<precision>       operator+ (const ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>  operator+ (const ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleComponent<precision>       operator+ (const ValueT&                            value,
                                                  const ParticleComponent<precision>&      vec   );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>  operator+ (const ValueT&                            value,
                                                  const ParticleCloudComponent<precision>& vec   );

    template<typename precision>
    ParticleComponent<precision>&      operator+=(      ParticleComponent<precision>&      first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision>
    ParticleCloudComponent<precision>& operator+=(      ParticleCloudComponent<precision>& first,
                                                  const ParticleCloudComponent<precision>& second);
    template<typename precision>
    ParticleCloudComponent<precision>& operator+=(      ParticleCloudComponent<precision>& first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision, typename ValueT>
    ParticleComponent<precision>&      operator+=(      ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>& operator+=(      ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );

    template<typename precision>
    ParticleComponent<precision>       operator- (const ParticleComponent<precision>&      first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision>
    ParticleCloudComponent<precision>  operator- (const ParticleCloudComponent<precision>& first,
                                                  const ParticleCloudComponent<precision>& second);
    template<typename precision>
    ParticleCloudComponent<precision>  operator- (const ParticleCloudComponent<precision>& first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision, typename ValueT>
    ParticleComponent<precision>       operator- (const ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>  operator- (const ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleComponent<precision>       operator- (const ValueT&                            value,
                                                  const ParticleComponent<precision>&      vec   );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>  operator- (const ValueT&                            value,
                                                  const ParticleCloudComponent<precision>& vec   );

    template<typename precision>
    ParticleComponent<precision>&      operator-=(      ParticleComponent<precision>&      first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision>
    ParticleCloudComponent<precision>& operator-=(      ParticleCloudComponent<precision>& first,
                                                  const ParticleCloudComponent<precision>& second);
    template<typename precision>
    ParticleCloudComponent<precision>& operator-=(      ParticleCloudComponent<precision>& first,
                                                  const ParticleComponent<precision>&      second);
    template<typename precision, typename ValueT>
    ParticleComponent<precision>&      operator-=(      ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>& operator-=(      ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );

    template<typename precision, typename ValueT>
    ParticleComponent<precision>       operator* (const ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>  operator* (const ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleComponent<precision>       operator* (const ValueT&                            value,
                                                  const ParticleComponent<precision>&      vec   );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>  operator* (const ValueT&                            value,
                                                  const ParticleCloudComponent<precision>& vec   );

    template<typename precision, typename ValueT>
    ParticleComponent<precision>&      operator*=(      ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>& operator*=(      ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );

    template<typename precision, typename ValueT>
    ParticleComponent<precision>       operator/ (const ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>  operator/ (const ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );

    template<typename precision, typename ValueT>
    ParticleComponent<precision>&      operator/=(      ParticleComponent<precision>&      vec,
                                                  const ValueT&                            value );
    template<typename precision, typename ValueT>
    ParticleCloudComponent<precision>& operator/=(      ParticleCloudComponent<precision>& vec,
                                                  const ValueT&                            value );
  }  // ::examples::boris
}  // ::examples


#include "particle_util_impl.hpp"

#endif  // _EXAMPLES__BORIS__PARTICLE_UTIL_HPP_
