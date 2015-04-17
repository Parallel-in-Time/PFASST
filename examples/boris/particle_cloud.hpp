/**
 * @file examples/boris/particle_cloud.hpp
 * @ingroup BorisFiles
 */
#ifndef _EXAMPLES__BORIS__PARTICLE_CLOUD_HPP_
#define _EXAMPLES__BORIS__PARTICLE_CLOUD_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/logging.hpp>
#include <pfasst/encap/encapsulation.hpp>
#ifdef WITH_MPI
  #include <pfasst/mpi_communicator.hpp>
#endif

#include "particle.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      /**
       * @ingroup Boris
       */
      template<typename precision>
      using ParticleCloudComponent = vector<precision>;


      /**
       * @ingroup Boris
       */
      template<typename precision>
      class ParticleCloud
        :   public encap::Encapsulation<precision>
          , public el::Loggable
      {
        private:
          size_t _dim;
          size_t _num_particles;
          ParticleCloudComponent<precision> _positions;
          ParticleCloudComponent<precision> _velocities;
          vector<precision> _charges;
          vector<precision> _masses;

          precision _default_charge;
          precision _default_mass;

#ifdef WITH_MPI
          //! @{
          vector<MPI_Request> recv_request;
          vector<MPI_Request> send_request;
          //! @}

          //! @{
          inline mpi::MPICommunicator& as_mpi(ICommunicator* comm);
          //! @}
#endif

        public:
          explicit ParticleCloud(const size_t num_particles = 0,
                                 const size_t dim = 3,
                                 const precision default_charge = precision(1.0),
                                 const precision default_mass = precision(1.0));
          virtual ~ParticleCloud();

          virtual void zero() override;
          virtual void copy(shared_ptr<const encap::Encapsulation<precision>> other);

          inline size_t size() const;
          inline size_t dim() const;
                ParticleCloudComponent<precision>& positions();
          const ParticleCloudComponent<precision>& positions() const;
                ParticleCloudComponent<precision>& velocities();
          const ParticleCloudComponent<precision>& velocities() const;
                vector<precision>& charges();
          const vector<precision>& charges() const;
                vector<precision>& masses();
          const vector<precision>& masses() const;

          ParticleComponent<precision> center_of_mass() const;
          // !! EXPENSIVE !!
          shared_ptr<Particle<precision>> operator[](const size_t index) const;
          // !! EXPENSIVE !!
          shared_ptr<Particle<precision>> at(const size_t index) const;
          void set_at(const size_t index, const shared_ptr<Particle<precision>>& particle);

          // !! VERY !! EXPENSIVE !! (i.e. never use for production code)
          vector<shared_ptr<Particle<precision>>> particles() const;

          void distribute_around_center(const shared_ptr<Particle<precision>>& center);

          // TODO: unify behaviour with particle_util::norm0 (e.g. norm_max vs. norm0 (==sqrt(^2))
          virtual precision norm0() const;

#ifdef WITH_MPI
          //! @{
          virtual void post(ICommunicator* comm, int tag) override;
          virtual void recv(ICommunicator* comm, int tag, bool blocking) override;
          virtual void send(ICommunicator* comm, int tag, bool blocking) override;
          virtual void broadcast(ICommunicator* comm) override;
          //! @}
#endif

          virtual void log(el::base::type::ostream_t& os) const;
      };


      /**
       * @ingroup BorisUtilities
       */
      template<typename precision>
      static precision distance(const Particle<precision>& first,
                                const Particle<precision>& second);
      /**
       * @ingroup BorisUtilities
       */
      template<typename precision>
      static precision distance(const shared_ptr<Particle<precision>>& first,
                                const shared_ptr<Particle<precision>>& second);

      /**
       * @ingroup BorisUtilities
       */
      template<typename precision>
      static vector<precision> distance_to_reference(const ParticleCloud<precision>& cloud,
                                                     const Particle<precision>&      reference);
      /**
       * @ingroup BorisUtilities
       */
      template<typename precision>
      static vector<precision> distance_to_reference(const shared_ptr<ParticleCloud<precision>>& cloud,
                                                     const shared_ptr<Particle<precision>>&      reference);


      /**
       * @ingroup BorisUtilities
       */
      template<typename precision>
      inline MAKE_LOGGABLE(shared_ptr<ParticleCloud<precision>>, sp_cloud, os);
      /**
       * @ingroup BorisUtilities
       */
      template<typename precision>
      inline MAKE_LOGGABLE(shared_ptr<const ParticleCloud<precision>>, sp_cloud, os);


      /**
       * @ingroup Boris
       */
      template<typename precision>
      class ParticleCloudFactory
        : public encap::EncapFactory<precision>
      {
        private:
          size_t _num_particles;
          size_t _dim;
          precision _default_charge;
          precision _default_mass;

        public:
          ParticleCloudFactory(const size_t num_particles, const size_t dim, const precision default_charge,
                               const precision default_mass);
          inline size_t num_particles() const;
          inline size_t dim() const;
          virtual shared_ptr<encap::Encapsulation<precision>> create(const encap::EncapType);
      };
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst


// #include "particle_util.hpp"
#include "particle_cloud_impl.hpp"

#endif  // _EXAMPLES__BORIS__PARTICLE_CLOUD_HPP_
