/**
 * @defgroup ExternalSolver External Solver
 * @ingroup BorisBindings
 */
#ifndef _SIMPLE_PHYSICS_SOLVER_HPP_
#define _SIMPLE_PHYSICS_SOLVER_HPP_

#include <cstdlib>
using namespace std;

#ifndef DIM
  #define DIM 3
#endif


//! @ingroup ExternalSolver
namespace simple_physics_solver
{
  //! @ingroup ExternalSolver
  class SimplePhysicsSolverConfig
  {
    public:
      double omega_e;
      double omega_b;
      double epsilon;
      double sigma;
      double sigma2;
      double external_e_field_matrix[3][3];
      double b_field_matrix[3][3];

      SimplePhysicsSolverConfig(const double omega_e = double(-4.9), const double omega_b = double(25.0),
                                const double epsilon = double(-1.0), const double sigma = double(0.1));
      virtual ~SimplePhysicsSolverConfig();
  };

  //! @ingroup ExternalSolver
  void evaluate_external_e_field(const double* positions, const double* charges, const double* masses,
                                 const size_t num_particles, const double t,
                                 const SimplePhysicsSolverConfig* config,
                                 double* forces);

  //! @ingroup ExternalSolver
  void evaluate_internal_e_field(const double* positions, const double* charges, const double* masses,
                                 const size_t num_particles, const double t,
                                 const SimplePhysicsSolverConfig* config,
                                 double* exyz, double* phis);

  //! @ingroup ExternalSolver
  void evaluate_e_field(const double* positions, const double* charges, const double* masses,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config,
                        double* forces);

  //! @ingroup ExternalSolver
  void get_b_field_vector(const SimplePhysicsSolverConfig* config,
                                 double* b_field_vector);

  //! @ingroup ExternalSolver
  void evaluate_b_field(const double* velocities, const double* masses, const double* charges,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config,
                        double* forces);

  //! @ingroup ExternalSolver
  double compute_energy(const double* positions, const double* velocities, const double* masses,
                        const double* charges,
                        const size_t num_particles, const double t,
                        const SimplePhysicsSolverConfig* config);

  //! @ingroup ExternalSolver
  namespace internal
  {
    inline void cross_prod(const double first[DIM], const double second[DIM], double cross_prod[DIM]);

    inline double scalar_prod(const double first[DIM], const double second[DIM]);

    inline void scale_mat_mul_vec(const double mat[DIM][DIM], const double vec[DIM], const double factor, double prod[DIM]);

    inline void print_vec(const double vec[DIM]);
  }  // ::simple_physics_solver::internal
}  // ::simple_physics_solver

#endif
