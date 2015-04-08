#include "simple_physics_solver.hpp"

#include <cassert>
#include <iostream>
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
//     cout << "SimplePhysicsSolver::evaluate_external_e_field(t=" << t << ")" << endl;
    double pre_factor = (- config->epsilon) * (config->omega_e * config->omega_e);
    double factor = double(0.0);
    for (size_t i = 0; i < num_particles; ++i) {
//       cout << "  particle " << i << endl;
      factor = pre_factor / (charges[i] / masses[i]);
//       cout << "    position: ";
//       internal::print_vec(positions + (i * DIM));
//       cout << endl << "    factor: " << factor << endl;
      internal::scale_mat_mul_vec(config->external_e_field_matrix, positions + (i * DIM), factor, forces + (i * DIM));
//       cout << "    force: ";
//       internal::print_vec(forces + (i * DIM));
//       cout << endl;
    }
//     cout << "DONE SimplePhysicsSolver::evaluate_external_e_field(t=" << t << ")" << endl;
  }

  void evaluate_internal_e_field(const double* positions, const double* charges, const double* masses,
                                 const size_t num_particles, const double t,
                                 const SimplePhysicsSolverConfig* config,
                                 double* exyz, double* phis)
  {
    UNUSED(masses); UNUSED(t);
//     cout << "SimplePhysicsSolver::evaluate_internal_e_field(t=" << t << ")" << endl;
    double r = double(0.0),
           r3 = double(0.0),
           dist2 = double(0.0);
    double* dist = new double[DIM];

    // computing forces on particle i
    for (size_t i = 0; i < num_particles; ++i) {
//       cout << "  particle " << i << endl;

      // null result values
      phis[i] = double(0.0);
      for (size_t d = 0; d < DIM; ++d) { exyz[i * DIM + d] = double(0.0); }

      // effects of particle j on particle i
      for (size_t j = 0; j < num_particles; ++j) {
//         cout << "    particle " << j << endl;
        if (j == i) {
//           cout << "      skipping" << endl;
          continue;
        }
        dist2 = double(0.0);
        for (size_t d = 0; d < DIM; ++d) {
          dist[d] = positions[i * DIM + d] - positions[j * DIM + d];
          dist2 += dist[d] * dist[d];
        }
//         cout << "      dist = " << dist << endl;
        r = sqrt(dist2 + config->sigma2);
//         cout << "      r = " << r << " (= sqrt(dist^2+" << config->sigma2 << ")" << endl;
        phis[i] += charges[j] / r;

        r3 = r * r * r;
        for (size_t d = 0; d < DIM; ++d) {
          exyz[i * DIM + d] += dist[d] / r3 * charges[j];
//           cout << "      exyz[" << i * DIM + d << "] += " << positions[j * DIM + d] << " / " << r3 << " * " << charges[j] << "(q) => " << exyz[i * DIM + d] << endl;
        }
      }

//       cout << "    exyz = ";
//       internal::print_vec(exyz+(i*DIM));
//       cout << endl
//            << "    phi_i = " << phis[i] << endl;
    }
    delete[] dist;
//     cout << "DONE SimplePhysicsSolver::evaluate_internal_e_field(t=" << t << ")" << endl;
  }


  void evaluate_e_field(const double* positions, const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config,
                        double* forces)
  {
//     cout << "SimplePhysicsSolver::evaluate_e_field(t=" << t << ")" << endl;
    double* external = new double[num_particles * DIM];
    double* internal = new double[num_particles * DIM];
    double* phis = new double[num_particles];
    evaluate_external_e_field(positions, charges, masses, num_particles, t, config, external);
    evaluate_internal_e_field(positions, charges, masses, num_particles, t, config, internal, phis);
//     cout << "  -> e_forces = [ " << endl;
    for (size_t i = 0; i < num_particles; ++i) {
//       cout << "                  ";
      for (size_t d = 0; d < DIM; ++d) {
        forces[i * DIM + d] = external[i * DIM + d] + internal[i * DIM + d];
      }
//       internal::print_vec(forces + (i * DIM));
//       cout << endl;
    }
//     cout << "                ]" << endl;
    delete[] external;
    delete[] internal;
    delete[] phis;
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
//     cout << "SimplePhysicsSolver::evaluate_b_field(t=" << t << ")" << endl;
    for (size_t i = 0; i < num_particles; ++i) {
//       cout << "  particle " << i << endl;
      internal::scale_mat_mul_vec(config->b_field_matrix, velocities + (i * DIM),
                                  config->omega_b / (charges[i] / masses[i]),
                                  forces + (i * DIM));
//       cout << "    force = ";
//       internal::print_vec(forces + (i * DIM));
//       cout << endl;
    }
  }


  double compute_energy(const double* positions, const double* velocities,
                        const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config)
  {
//     cout << "SimplePhysicsSolver::compute_energy(t=" << t << ")" << endl;
    double e_kin = double(0.0);
    double e_pot = double(0.0);

    double* exyz = new double[num_particles * DIM];
    double* phis = new double[num_particles];
    double* temp = new double[DIM];
    double v2 = double(0.0);

    evaluate_internal_e_field(positions, charges, masses, num_particles, t, config, exyz, phis);

    for (size_t i = 0; i < num_particles; ++i) {
//       cout << "  particle " << i << endl;
      // potential energy (induced by external electric field, position and internal electric field (i.e. phis))
      internal::scale_mat_mul_vec(config->external_e_field_matrix, positions + (i * DIM),
                                  (- config->epsilon * config->omega_e * config->omega_e / double(2.0) * (charges[i] / masses[i])),
                                  temp);
      e_pot += (charges[i] * phis[i]) - internal::scalar_prod(temp, positions + (i * DIM));
//       cout << "    e_pot += " << charges[i] << " * " << phis[i] << " - <";
//       internal::print_vec(temp);
//       cout << ", ";
//       internal::print_vec(positions + (i * DIM));
//       cout << "> (=" << internal::scalar_prod(temp, positions + (i * DIM)) << ")" << endl;

      // kinetic energy (induced by velocity)
      v2 = internal::scalar_prod(velocities + (i * DIM), velocities + (i * DIM));
      e_kin += masses[i] / double(2.0) * v2;
//       cout << "    e_kin += " << masses[i] << "(m) / 2.0" << " * <";
//       internal::print_vec(velocities + (i * DIM));
//       cout << " , ";
//       internal::print_vec(velocities + (i * DIM));
//       cout << "> (=" << v2 << ")" << endl;
    }

    delete[] exyz;
    delete[] phis;
    delete[] temp;
//     cout << "  -> energy = " << e_kin << " + " << e_pot << endl;
//     cout << "DONE SimplePhysicsSolver::compute_energy(t=" << t << ")" << endl;
    return e_kin + e_pot;
  }

  namespace internal
  {
    inline void cross_prod(const double first[DIM], const double second[DIM], double cross_prod[DIM])
    {
      cross_prod[0] = first[1] * second[2] - first[2] * second[1];
      cross_prod[1] = first[2] * second[0] - first[0] * second[2];
      cross_prod[2] = first[0] * second[1] - first[1] * second[0];
    }

    inline double scalar_prod(const double first[DIM], const double second[DIM])
    {
      double prod = 0.0;
      for (size_t m = 0; m < DIM; ++m) {
        prod += first[m] * second[m];
      }
      return prod;
    }

    inline void scale_mat_mul_vec(const double mat[DIM][DIM], const double vec[DIM],
                                  const double factor, double prod[DIM])
    {
      for (size_t i = 0; i < DIM; ++i) {
        prod[i] = double(0.0);
        for (size_t j = 0; j < DIM; ++j) {
          prod[i] += factor * mat[i][j] * vec[j];
        }
      }
    }

    inline void print_vec(const double vec[DIM])
    {
      cout << "[";
      for (size_t i = 0; i < DIM; ++i) {
        cout << vec[i];
        if (i < DIM - 1) {
          cout << " , ";
        }
      }
      cout << "]";
    }
  }  // ::simple_physics_solver::internal
}  // ::simple_physics_solver
