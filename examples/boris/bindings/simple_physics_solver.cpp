#include "simple_physics_solver.hpp"

#include <cassert>
#include <cmath>
using namespace std;

#define UNUSED(expr) (void)(expr)

namespace simple_physics_solver
{
  template<typename scalar, typename time>
  SimplePhysicsSolverConfig<scalar, time>::SimplePhysicsSolverConfig(const scalar omega_e, const scalar omega_b,
                                                                     const scalar epsilon, const scalar sigma)
    :   omega_e(omega_e)
      , omega_b(omega_b)
      , epsilon(epsilon)
      , sigma(sigma)
      , sigma2(sigma * sigma)
  {
    this->external_e_field_matrix[0][0] = scalar(1.0);
    this->external_e_field_matrix[0][1] = scalar(0.0);
    this->external_e_field_matrix[0][2] = scalar(0.0);
    this->external_e_field_matrix[1][0] = scalar(0.0);
    this->external_e_field_matrix[1][1] = scalar(1.0);
    this->external_e_field_matrix[1][2] = scalar(0.0);
    this->external_e_field_matrix[2][0] = scalar(0.0);
    this->external_e_field_matrix[2][1] = scalar(0.0);
    this->external_e_field_matrix[2][2] = scalar(-2.0);

    this->b_field_matrix[0][0] = scalar(0.0);
    this->b_field_matrix[0][1] = scalar(1.0);
    this->b_field_matrix[0][2] = scalar(0.0);
    this->b_field_matrix[1][0] = scalar(-1.0);
    this->b_field_matrix[1][1] = scalar(0.0);
    this->b_field_matrix[1][2] = scalar(0.0);
    this->b_field_matrix[2][0] = scalar(0.0);
    this->b_field_matrix[2][1] = scalar(0.0);
    this->b_field_matrix[2][2] = scalar(0.0);
  }

  template<typename scalar, typename time>
  SimplePhysicsSolverConfig<scalar, time>::~SimplePhysicsSolverConfig()
  {}


  template<typename scalar, typename time>
  static void evaluate_external_e_field(const scalar** positions, const scalar* charges, const scalar* masses,
                                        const size_t num_particles, const size_t dim, const time t,
                                        const SimplePhysicsSolverConfig<scalar, time>& config,
                                        scalar** forces)
  {
    UNUSED(t);
    scalar pre_factor = (- config.epsilon) * (config.omega_e * config.omega_e);
    scalar factor;
    for (size_t i = 0; i < num_particles; ++i) {
      factor = pre_factor / (charges[i] / masses[i]);
      internal::scale_mat_mul_vec(config.external_e_field_matrix, positions[i], dim, dim, factor, forces[i]);
    }
  }

  template<typename scalar, typename time>
  static void evaluate_internal_e_field(const scalar** positions, const scalar* charges, const scalar* masses,
                                        const size_t num_particles, const size_t dim, const time t,
                                        const SimplePhysicsSolverConfig<scalar, time>& config,
                                        scalar** exyz, scalar* phis)
  {
    UNUSED(masses); UNUSED(t);
    scalar r = scalar(0.0),
           r3 = scalar(0.0),
           dist = scalar(0.0);

    // computing forces on particle i
    for (size_t i = 0; i < num_particles; ++i) {

      // null result values
      phis[i] = scalar(0.0);
      for (size_t d = 0; d < dim; ++d) { exyz[i][d] = scalar(0.0); }

      // effects of particle j on particle i
      for (size_t j = 0; j < num_particles; ++j) {
        dist = scalar(0.0);
        for (size_t d = 0; d < dim; ++d) {
          dist += (positions[i * dim + d] - positions[j * dim + d]) * (positions[i * dim + d] - positions[j * dim + d]);
        }
        r = sqrt(dist * dist + config.sigma2);
        phis[i] += charges[j] / r;

        r3 = r * r * r;
        for (size_t d = 0; d < dim; ++d) {
          exyz[i][d] += positions[j][d] / r3 * charges[j];
        }
      }
    }
  }


  template<typename scalar, typename time>
  static void get_b_field_vector(const SimplePhysicsSolverConfig<scalar, time>& config, const size_t dim,
                                 scalar* b_field_vector)
  {
    for (size_t i = 0; i < dim; ++i) {
      b_field_vector[i] = time(0.0);
    }
    b_field_vector[2] = time(1.0);
    for (size_t i = 0; i < dim; ++i) {
      b_field_vector[i] *= config.omega_b;
    }
  }

  template<typename scalar, typename time>
  static void evaluate_b_field(const scalar** velocities, const scalar* masses, const scalar* charges,
                               const size_t num_particles, const size_t dim, const time t,
                               const SimplePhysicsSolverConfig<scalar, time>& config, scalar** forces)
  {
    UNUSED(t);
    for (size_t i = 0; i < num_particles; ++i) {
      internal::scale_mat_mul_vec(config.b_field_matrix, velocities[i], dim, dim,
                                  config.omega_b / (charges[i] / masses[i]),
                                  forces[i]);
    }
  }


  template<typename scalar, typename time>
  static scalar compute_energy(const scalar** positions, const scalar** velocities, const scalar* masses,
                               const scalar* charges,
                               const size_t num_particles, const size_t dim, const time t,
                               const SimplePhysicsSolverConfig<scalar, time>& config)
  {
    assert(dim == 3);
    scalar e_kin = scalar(0.0);
    scalar e_pot = scalar(0.0);

    scalar phi_i = scalar(0.0),
           r = scalar(0.0),
           dist = scalar(0.0);
    scalar cross_prod[dim];

    scalar exyz[num_particles][dim];
    scalar phis[num_particles];

    evaluate_internal_e_field(positions, charges, masses, num_particles, dim, t, config, exyz, phis);

    for (size_t i = 0; i < num_particles; ++i) {
      e_pot += charges[i] * phis[i];
      internal::cross_prod(velocities[i], velocities[i], cross_prod);
      for (size_t m = 0; m < dim; ++m) {
        e_kin += masses[i] / scalar(2.0) * cross_prod[m];
      }
    }

    return e_kin + e_pot;
  }

  namespace internal
  {
    template<typename scalar>
    static void cross_prod(const scalar* first, const scalar* second, scalar* cross_prod)
    {
      cross_prod[0] = first[1] * second[2] - first[2] * second[1];
      cross_prod[1] = first[2] * second[0] - first[0] * second[2];
      cross_prod[2] = first[0] * second[1] - first[1] * second[0];
    }

    template<typename scalar>
    static void scale_mat_mul_vec(const scalar** mat, const scalar* vec, const size_t rows, const size_t cols,
                                  const scalar factor, scalar* prod)
    {
      for (size_t i = 0; i < rows; ++i) {
        prod[i] = scalar(0.0);
        for (size_t j = 0; j < cols; ++j) {
          prod[i] += factor * mat[i][j] * vec[j];
        }
      }
    }
  }  // ::simple_physics_solver::internal
}  // ::simple_physics_solver
