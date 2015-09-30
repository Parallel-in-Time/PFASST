#include "particle_cloud.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <random>
#include <string>
using namespace std;

#include <boost/format.hpp>

#include <pfasst/globals.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#ifdef WITH_MPI
  #include <mpi.h>
#endif

#include "particle_util.hpp"

#define PFASST_RANDOM_SEED 42


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
#ifdef WITH_MPI
      template<typename precision>
      inline mpi::MPICommunicator& ParticleCloud<precision>::as_mpi(ICommunicator* comm)
      {
        auto mpi = dynamic_cast<mpi::MPICommunicator*>(comm);
        assert(mpi);
        return *mpi;
      }
#endif


      template<typename precision>
      ParticleCloud<precision>::ParticleCloud(const size_t num_particles, const size_t dim,
                                              const precision default_charge, const precision default_mass)
        :   _dim(dim)
          , _num_particles(num_particles)
          , _positions(num_particles * dim)
          , _velocities(num_particles * dim)
          , _charges(num_particles, default_charge)
          , _masses(num_particles, default_mass)
          , _default_charge(default_charge)
          , _default_mass(default_mass)
#ifdef WITH_MPI
          , recv_request(2)
          , send_request(2)
#endif
      {
        this->zero();
      }

      template<typename precision>
      ParticleCloud<precision>::~ParticleCloud()
      {}

      template<typename precision>
      void ParticleCloud<precision>::zero()
      {
        fill(this->_positions.begin(), this->_positions.end(), precision(0.0));
        fill(this->_velocities.begin(), this->_velocities.end(), precision(0.0));
        fill(this->_charges.begin(), this->_charges.end(), this->_default_charge);
        fill(this->_masses.begin(), this->_masses.end(), this->_default_mass);
#ifdef WITH_MPI
        fill(this->recv_request.begin(), this->recv_request.end(), MPI_REQUEST_NULL);
        fill(this->send_request.begin(), this->send_request.end(), MPI_REQUEST_NULL);
#endif
      }

      template<typename precision>
      void ParticleCloud<precision>::copy(shared_ptr<const encap::Encapsulation<precision>> other)
      {
        shared_ptr<const ParticleCloud<precision>> other_c = dynamic_pointer_cast<const ParticleCloud<precision>>(other);
        assert(other_c);
//         this->copy(other_c);
//       }

//       template<typename precision>
//       void ParticleCloud<precision>::copy(shared_ptr<const ParticleCloud<precision>> other)
//       {
        this->_dim = other_c->dim();
        this->_num_particles = other_c->size();
        this->_positions = other_c->positions();
        this->_velocities = other_c->velocities();
        this->_charges = other_c->charges();
        this->_masses = other_c->masses();
        this->_default_charge = other_c->_default_charge;
        this->_default_mass = other_c->_default_mass;
      }

      /*
      template<typename precision>
      void ParticleCloud<precision>::extend(const size_t new_size)
      {
        if (new_size <= this->_num_particles) {
          throw invalid_argument("new size (" + to_string(new_size) + ") is not larger current size (" 
                                 + to_string(this->_num_particles) + ")");
        }

        this->_positions.resize(new_size * this->dim());
        fill(this->_positions.begin() + this->_num_particles, this->_positions.end(), precision(0.0));
        this->_velocities.resize(new_size * this->dim());
        fill(this->_velocities.begin() + this->_num_particles, this->_velocities.end(), precision(0.0));
        this->_masses.resize(new_size);
        fill(this->_masses.begin() + this->_num_particles, this->_masses.end(), this->_default_mass);
        this->_charges.resize(new_size);
        fill(this->_charges.begin() + this->_num_particles, this->_charges.end(), this->_default_charge);

        this->_num_particles = new_size;
      }

      template<typename precision>
      void ParticleCloud<precision>::erase(const size_t pos)
      {
        if (pos >= this->_num_particles) {
          throw out_of_range("index to erase (" + to_string(pos) + ") is out of bounds ("
                             + to_string(this->_num_particles) + ")");
        }

        this->_positions.erase(this->_positions.cbegin() + pos, this->_positions.cbegin() + pos + 1);
        this->_velocities.erase(this->_velocities.cbegin() + pos, this->_velocities.cbegin() + pos + 1);
        this->_charges.erase(this->_charges.cbegin() + pos, this->_charges.cbegin() + pos + 1);
        this->_masses.erase(this->_masses.cbegin() + pos, this->_masses.cbegin() + pos + 1);
        this->_num_particles--;
      }

      template<typename precision>
      void ParticleCloud<precision>::push_back(const shared_ptr<Particle<precision>>& particle)
      {
        if (particle->DIM() != this->_dim) {
          throw invalid_argument("particles dimension (" + to_string(particle->DIM())
                                 + ") does not match cloud dimension ("+to_string(this->_dim)+")");
        }

        this->_positions.push_back(particle->pos());
        this->_velocities.push_back(particle->vel());
        this->_charges.push_back(particle->charge());
        this->_masses.push_back(particle->mass());
        this->_num_particles++;
      }

      template<typename precision>
      void ParticleCloud<precision>::insert(const size_t pos, const shared_ptr<Particle<precision>>& particle)
      {
        if (particle->DIM() != this->_dim) {
          throw invalid_argument("particles dimension (" + to_string(particle->DIM())
                                 + ") does not match cloud dimension ("+to_string(this->_dim)+")");
        }

        this->_positions.insert(pos, particle->pos());
        this->_velocities.insert(pos, particle->vel());
        this->_charges.insert(pos, particle->charge());
        this->_masses.insert(pos, particle->mass());
        this->_num_particles++;
      }
      */

      template<typename precision>
      inline size_t ParticleCloud<precision>::size() const
      {
        return this->_num_particles;
      }

      template<typename precision>
      inline size_t ParticleCloud<precision>::dim() const
      {
        return this->_dim;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& ParticleCloud<precision>::positions()
      {
        return this->_positions;
      }

      template<typename precision>
      const ParticleCloudComponent<precision>& ParticleCloud<precision>::positions() const
      {
        return this->_positions;
      }

      template<typename precision>
      ParticleCloudComponent<precision>& ParticleCloud<precision>::velocities()
      {
        return this->_velocities;
      }

      template<typename precision>
      const ParticleCloudComponent<precision>& ParticleCloud<precision>::velocities() const
      {
        return this->_velocities;
      }

      template<typename precision>
      vector<precision>& ParticleCloud<precision>::charges()
      {
        return this->_charges;
      }

      template<typename precision>
      const vector<precision>& ParticleCloud<precision>::charges() const 
      {
        return this->_charges;
      }

      template<typename precision>
      vector<precision>& ParticleCloud<precision>::masses()
      {
        return this->_masses;
      }

      template<typename precision>
      const vector<precision>& ParticleCloud<precision>::masses() const
      {
        return this->_masses;
      }

      template<typename precision>
      vector<shared_ptr<Particle<precision>>> ParticleCloud<precision>::particles() const
      {
        vector<shared_ptr<Particle<precision>>> particles(this->size());
        for (size_t index = 0; index < this->size(); ++index) {
          particles[index] = this->operator[](index);
        }
        return particles;
      }

      template<typename precision>
      ParticleComponent<precision> ParticleCloud<precision>::center_of_mass() const
      {
        ParticleComponent<precision> center(this->dim());
        for (size_t p = 0; p < this->size(); ++p) {
          // conter += position[p]
          std::transform(center.cbegin(), center.cend(), this->positions().cbegin() + (p * this->dim()),
                         center.begin(), std::plus<precision>());
        }
        return center / this->size();
      }

      template<typename precision>
      shared_ptr<Particle<precision>> ParticleCloud<precision>::operator[](const size_t index) const
      {
        shared_ptr<Particle<precision>> particle = make_shared<Particle<precision>>(this->dim());
        copy_n(this->positions().cbegin() + (index * this->dim()), this->dim(), particle->pos().begin());
        copy_n(this->velocities().cbegin() + (index * this->dim()), this->dim(), particle->vel().begin());
        particle->set_charge(this->_charges[index]);
        particle->set_mass(this->_masses[index]);
        return particle;
      }

      template<typename precision>
      shared_ptr<Particle<precision>> ParticleCloud<precision>::at(const size_t index) const
      {
        assert(this->size() > index);
        return this->operator[](index);
      }

      template<typename precision>
      void ParticleCloud<precision>::set_at(const size_t index, const shared_ptr<Particle<precision>>& particle)
      {
        assert(this->size() > index);
        assert(particle->dim() == this->dim());
        copy_n(particle->pos().cbegin(), this->dim(), this->positions().begin() + (index * this->dim()));
        copy_n(particle->vel().cbegin(), this->dim(), this->velocities().begin() + (index * this->dim()));
        this->_masses[index] = particle->mass();
        this->_charges[index] = particle->charge();
      }

      template<typename precision>
      void ParticleCloud<precision>::distribute_around_center(const shared_ptr<Particle<precision>>& center)
      {
        VLOG_FUNC_START("ParticleCloud") << " center:" << center;
        VLOG(3) << LOG_INDENT << "distributing " << this->size() << " particles around center " << center;
        assert(this->size() > 0);

        precision scale = 1000.0;

        default_random_engine rd_gen(PFASST_RANDOM_SEED);
        precision max_pos = max(center->pos());
        precision max_vel = max(center->vel());
        uniform_real_distribution<precision> dist_pos(- max_pos / scale, max_pos / scale);
        uniform_real_distribution<precision> dist_vel(- max_vel / scale, max_vel / scale);
        VLOG(4) << LOG_INDENT << "random displacement range for";
        VLOG(4) << LOG_INDENT << " ... position: " << boost::format("[%.4f, %.4f]") % dist_pos.min() % dist_pos.max();
        VLOG(4) << LOG_INDENT << " ... velocity: " << boost::format("[%.4f, %.4f]") % dist_vel.min() % dist_vel.max();

        size_t p = 0;

        if (this->size() % 2 == 1) {
          this->set_at(p, center);
          VLOG(5) << LOG_INDENT << "setting p=1 to center";
          p++;
        }
        for (;p < this->size(); ++p) {
          ParticleComponent<precision> pos_rand(this->dim());
          ParticleComponent<precision> vel_rand(this->dim());
          std::transform(center->pos().cbegin(), center->pos().cend(), pos_rand.begin(),
                         [&](const precision& c) { return c + dist_pos(rd_gen); });
          std::transform(center->pos().cbegin(), center->pos().cend(), vel_rand.begin(),
                         [&](const precision& c) { return c + dist_vel(rd_gen); });
          std::copy(pos_rand.cbegin(), pos_rand.cend(), this->positions().begin() + (p * this->dim()));
          std::copy(vel_rand.cbegin(), vel_rand.cend(), this->velocities().begin() + (p * this->dim()));
          VLOG(5) << LOG_INDENT << "p=" << (p+1) << ": " << this->at(p);
        }
        assert(p == this->size());
        VLOG(3) << LOG_INDENT << "center after distribute: " << this->center_of_mass();
      }

      template<typename precision>
      precision ParticleCloud<precision>::norm0() const
      {
        return std::max(max_abs(this->positions()), max_abs(this->velocities()));
      }

#ifdef WITH_MPI
      template<typename precision>
      void ParticleCloud<precision>::post(ICommunicator* comm, int tag)
      {
        auto& mpi = as_mpi(comm);
        if (mpi.size() == 1) { return; }
        if (mpi.rank() == 0) { return; }

        int err = MPI_Irecv(this->positions().data(),
                            sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                            (mpi.rank() - 1) % mpi.size(), tag, mpi.comm, &(this->recv_request[0]));
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        err = MPI_Irecv(this->velocities().data(),
                        sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                        (mpi.rank() - 1) % mpi.size(), tag + 1, mpi.comm, &(this->recv_request[1]));
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
      }

      template<typename precision>
      void ParticleCloud<precision>::recv(ICommunicator* comm, int tag, bool blocking)
      {
        auto& mpi = as_mpi(comm);
        if (mpi.size() == 1) { return; }
        if (mpi.rank() == 0) { return; }

        int err;
        MPI_Status stat;
        if (blocking) {
          err = MPI_Recv(this->positions().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         (mpi.rank() - 1) % mpi.size(), tag, mpi.comm, &stat);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }

          err = MPI_Recv(this->velocities().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         (mpi.rank() - 1) % mpi.size(), tag + 1, mpi.comm, &stat);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        } else {
          for (auto req : this->recv_request) {
            err = MPI_Wait(&req, &stat);
            if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          }
        }
      }

      template<typename precision>
      void ParticleCloud<precision>::send(ICommunicator* comm, int tag, bool blocking)
      {
        auto& mpi = as_mpi(comm);
        if (mpi.size() == 1) { return; }
        if (mpi.rank() == mpi.size() - 1) { return; }

        int err = MPI_SUCCESS;
        if (blocking) {
          err = MPI_Send(this->positions().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         (mpi.rank() + 1) % mpi.size(), tag, mpi.comm);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          err = MPI_Send(this->velocities().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                         (mpi.rank() + 1) % mpi.size(), tag + 1, mpi.comm);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        } else {
          MPI_Status stat;
          for (auto req : this->send_request) {
            err = MPI_Wait(&req, &stat);
            if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          }

          err = MPI_Isend(this->positions().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                          (mpi.rank() + 1) % mpi.size(), tag, mpi.comm, &send_request[0]);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
          err = MPI_Isend(this->velocities().data(),
                         sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                          (mpi.rank() + 1) % mpi.size(), tag + 1, mpi.comm, &send_request[1]);
          if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        }
      }

      template<typename precision>
      void ParticleCloud<precision>::broadcast(ICommunicator* comm)
      {
        auto& mpi = as_mpi(comm);
        int err = MPI_Bcast(this->positions().data(),
                            sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                            comm->size() - 1, mpi.comm);
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
        err = MPI_Bcast(this->velocities().data(),
                        sizeof(precision) * this->size() * this->dim(), MPI_CHAR,
                        comm->size() - 1, mpi.comm);
        if (err != MPI_SUCCESS) { throw mpi::MPIError(); }
      }
#endif

      template<typename precision>
      void ParticleCloud<precision>::log(el::base::type::ostream_t& os) const
      {
        os << fixed << setprecision(LOG_PRECISION);
        os << "ParticleCloud(n=" << this->size() << ", d=" << this->dim() << ", particles=" << this->particles() << ")";
        os.unsetf(ios_base::floatfield);
      }


      template<typename precision>
      static precision distance(const Particle<precision>& first, const Particle<precision>& second)
      {
        assert(first.dim() == second.dim());
        vector<precision> diff = first.pos() - second.pos();
        return norm0(diff);
      }

      template<typename precision>
      static precision distance(const shared_ptr<Particle<precision>> first, const shared_ptr<Particle<precision>> second)
      {
        return distance(*(first.get()), *(second.get()));
      }

      template<typename precision>
      static vector<precision> distance_to_reference(const ParticleCloud<precision>& cloud,
                                                     const Particle<precision>&      reference)
      {
        vector<precision> distances(cloud.size());
        for (size_t i = 0; i < distances.size(); ++i) {
          distances[i] = distance(*cloud[i], reference);
        }
        return distances;
      }

      template<typename precision>
      static vector<precision> distance_to_reference(const shared_ptr<ParticleCloud<precision>>& cloud,
                                                     const shared_ptr<Particle<precision>>&      reference)
      {
        return distance_to_reference(*cloud, *reference);
      }


      template<typename precision>
      inline MAKE_LOGGABLE(shared_ptr<ParticleCloud<precision>>, sp_cloud, os)
      {
        os << "<" << addressof(sp_cloud) << ">";
        sp_cloud->log(os);
        return os;
      }

      template<typename precision>
      inline MAKE_LOGGABLE(shared_ptr<const ParticleCloud<precision>>, sp_cloud, os)
      {
        os << "<" << addressof(sp_cloud) << ">";
        sp_cloud->log(os);
        return os;
      }


      template<typename precision>
      ParticleCloudFactory<precision>::ParticleCloudFactory(const size_t num_particles, const size_t dim,
                                                            const precision default_charge,
                                                            const precision default_mass)
        :   _num_particles(num_particles)
          , _dim(dim)
          , _default_charge(default_charge)
          , _default_mass(default_mass)
      {}

      template<typename precision>
      inline size_t ParticleCloudFactory<precision>::num_particles() const
      {
        return this->_num_particles;
      }

      template<typename precision>
      inline size_t ParticleCloudFactory<precision>::dim() const
      {
        return this->_dim;
      }

      template<typename precision>
      shared_ptr<encap::Encapsulation<precision>>
      ParticleCloudFactory<precision>::create(const encap::EncapType)
      {
        return make_shared<ParticleCloud<precision>>(this->num_particles(), this->dim(),
                                                     this->_default_charge, this->_default_mass);
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst
