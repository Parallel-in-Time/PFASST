#include "simple_physics_solver.hpp"

#include <cassert>
#include <cmath>
using namespace std;

#define UNUSED(expr) (void)(expr)

namespace simple_physics_solver
{
  SimplePhysicsSolverConfig::SimplePhysicsSolverConfig(const double omega_e, const double omega_b,
                                                       const double epsilon, const double sigma)
    :   omega_e(omega_e)
      , omega_b(omega_b)
      , epsilon(epsilon)
      , sigma(sigma)
      , sigma2(sigma * sigma)
  {
    this->external_e_field_matrix[0][0] = double(1.0);
    this->external_e_field_matrix[0][1] = double(0.0);
    this->external_e_field_matrix[0][2] = double(0.0);
    this->external_e_field_matrix[1][0] = double(0.0);
    this->external_e_field_matrix[1][1] = double(1.0);
    this->external_e_field_matrix[1][2] = double(0.0);
    this->external_e_field_matrix[2][0] = double(0.0);
    this->external_e_field_matrix[2][1] = double(0.0);
    this->external_e_field_matrix[2][2] = double(-2.0);

    this->b_field_matrix[0][0] = double(0.0);
    this->b_field_matrix[0][1] = double(1.0);
    this->b_field_matrix[0][2] = double(0.0);
    this->b_field_matrix[1][0] = double(-1.0);
    this->b_field_matrix[1][1] = double(0.0);
    this->b_field_matrix[1][2] = double(0.0);
    this->b_field_matrix[2][0] = double(0.0);
    this->b_field_matrix[2][1] = double(0.0);
    this->b_field_matrix[2][2] = double(0.0);
  }

  SimplePhysicsSolverConfig::~SimplePhysicsSolverConfig()
  {}


  void evaluate_external_e_field(const double* positions, const double* charges, const double* masses,
                                 const size_t num_particles, const double t,
                                 const SimplePhysicsSolverConfig* config,
                                 double* forces)
  {
    UNUSED(t);
    double pre_factor = (- config->epsilon) * (config->omega_e * config->omega_e);
    double factor = double(0.0);
    for (size_t i = 0; i < num_particles; ++i) {
      factor = pre_factor / (charges[i] / masses[i]);
      internal::scale_mat_mul_vec(config->external_e_field_matrix, positions + i, factor, forces + i);
    }
  }

  void evaluate_internal_e_field(const double* positions, const double* charges, const double* masses,
                                 const size_t num_particles, const double t,
                                 const SimplePhysicsSolverConfig* config,
                                 double* exyz, double* phis)
  {
    UNUSED(masses); UNUSED(t);
    double r = double(0.0),
           r3 = double(0.0),
           dist = double(0.0);

    // computing forces on particle i
    for (size_t i = 0; i < num_particles; ++i) {

      // null result values
      phis[i] = double(0.0);
      for (size_t d = 0; d < DIM; ++d) { exyz[i * DIM + d] = double(0.0); }

      // effects of particle j on particle i
      for (size_t j = 0; j < num_particles; ++j) {
        dist = double(0.0);
        for (size_t d = 0; d < DIM; ++d) {
          dist += (positions[i * DIM + d] - positions[j * DIM + d]) * (positions[i * DIM + d] - positions[j * DIM + d]);
        }
        r = sqrt(dist * dist + config->sigma2);
        phis[i] += charges[j] / r;

        r3 = r * r * r;
        for (size_t d = 0; d < DIM; ++d) {
          exyz[i * DIM + d] += positions[j * DIM + d] / r3 * charges[j];
        }
      }
    }
  }


  void evaluate_e_field(const double* positions, const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config,
                        double* forces)
  {
    double external[num_particles * DIM];
    double internal[num_particles * DIM];
    double phis[num_particles];
    evaluate_external_e_field(positions, charges, masses, num_particles, t, config, external);
    evaluate_internal_e_field(positions, charges, masses, num_particles, t, config, internal, phis);
    for (size_t i = 0; i < num_particles; ++i) {
      for (size_t d = 0; d < DIM; ++d) {
        forces[i * DIM + d] = external[i * DIM + d] + internal[i * DIM + d];
      }
    }
  }


  void get_b_field_vector(const SimplePhysicsSolverConfig* config, double* b_field_vector)
  {
    for (size_t i = 0; i < DIM; ++i) {
      b_field_vector[i] = double(0.0);
    }
    b_field_vector[2] = double(1.0);
    for (size_t i = 0; i < DIM; ++i) {
      b_field_vector[i] *= config->omega_b;
    }
  }

  void evaluate_b_field(const double* velocities, const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config, double* forces)
  {
    UNUSED(t);
    for (size_t i = 0; i < num_particles; ++i) {
      internal::scale_mat_mul_vec(config->b_field_matrix, velocities + (i * DIM),
                                  config->omega_b / (charges[i] / masses[i]),
                                  forces + (i * DIM));
    }
  }


  double compute_energy(const double* positions, const double* velocities,
                        const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config)
  {
    double e_kin = double(0.0);
    double e_pot = double(0.0);

    double cross_prod[DIM];

    double exyz[num_particles * DIM];
    double phis[num_particles];

    evaluate_internal_e_field(positions, charges, masses, num_particles, t, config, exyz, phis);

    for (size_t i = 0; i < num_particles; ++i) {
      e_pot += charges[i] * phis[i];
      internal::cross_prod(velocities + (i * DIM), velocities + (i * DIM), cross_prod);
      for (size_t m = 0; m < DIM; ++m) {
        e_kin += masses[i] / double(2.0) * cross_prod[m];
      }
    }

    return e_kin + e_pot;
  }

  namespace internal
  {
    static void cross_prod(const double* first, const double* second, double* cross_prod)
    {
      cross_prod[0] = first[1] * second[2] - first[2] * second[1];
      cross_prod[1] = first[2] * second[0] - first[0] * second[2];
      cross_prod[2] = first[0] * second[1] - first[1] * second[0];
    }

    static void scale_mat_mul_vec(const double mat[DIM][DIM], const double vec[DIM],
                                  const double factor, double* prod)
    {
      for (size_t i = 0; i < DIM; ++i) {
        prod[i] = double(0.0);
        for (size_t j = 0; j < DIM; ++j) {
          prod[i] += factor * mat[i][j] * vec[j];
        }
      }
    }
  }  // ::simple_physics_solver::internal
}  // ::simple_physics_solver
