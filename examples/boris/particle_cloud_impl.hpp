#include "particle_cloud.hpp"

#include <cassert>
#include <string>
using namespace std;

namespace examples
{
  namespace boris
  {
    template<typename precision>
    ParticleCloud<precision>::ParticleCloud(const size_t num_particles, const size_t dim,
                                            const precision default_charge, const precision default_mass)
      :   _dim(dim)
        , _num_particles(num_particles)
        , _positions(num_particles)
        , _velocities(num_particles)
        , _charges(num_particles, default_charge)
        , _masses(num_particles, default_mass)
        , _default_charge(default_charge)
        , _default_mass(default_mass)
    {}

    template<typename precision>
    ParticleCloud<precision>::~ParticleCloud()
    {}

    template<typename precision>
    void ParticleCloud<precision>::extend(const size_t new_size)
    {
      if (new_size <= this->_num_particles) {
        throw invalid_argument("new size (" + to_string(new_size) + ") is not larger current size (" 
                               + to_string(this->_num_particles) + ")");
      }

      for (size_t new_elem = 0; new_elem < (new_size - this->_num_particles); ++new_elem) {
        this->_positions.push_back(vector<precision>(this->_dim));
        this->_velocities.push_back(vector<precision>(this->_dim));
        this->_charges.push_back(this->_default_charge);
        this->_masses.push_back(this->_default_mass);
      }
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

    template<typename precision>
    size_t ParticleCloud<precision>::size() const
    {
      return this->_num_particles;
    }

    template<typename precision>
    ParticleCloudComponent<precision>& ParticleCloud<precision>::positions()
    {
      return this->_positions;
    }

    template<typename precision>
    ParticleCloudComponent<precision>& ParticleCloud<precision>::velocities()
    {
      return this->_velocities;
    }

    template<typename precision>
    vector<precision>& ParticleCloud<precision>::charges()
    {
      return this->_charges;
    }

    template<typename precision>
    vector<precision>& ParticleCloud<precision>::masses()
    {
      return this->_masses;
    }

    template<typename precision>
    vector<shared_ptr<Particle<precision>>> ParticleCloud<precision>::particles() const
    {
      vector<shared_ptr<Particle<precision>>> particles(this->_num_particles);
      for (size_t index = 0; index < this->_num_particles; ++index) {
        particles[index] = this->operator[](index);
      }
      return particles;
    }

    template<typename precision>
    ParticleComponent<precision> ParticleCloud<precision>::center_of_mass() const
    {
      ParticleComponent<precision> center(this->_dim);
      for (auto elem : this->_positions) {
        center += elem;
      }
    }

    template<typename precision>
    shared_ptr<Particle<precision>> ParticleCloud<precision>::operator[](const size_t index) const
    {
      shared_ptr<Particle<precision>> particle = make_shared<Particle<precision>>(this->_dim);
      particle->pos() = this->_positions[index];
      particle->vel() = this->_velocities[index];
      particle->set_charge(this->_charges[index]);
      particle->set_mass(this->_masses[index]);
      return particle;
    }
  }  // ::examples::boris
}  // ::examples
