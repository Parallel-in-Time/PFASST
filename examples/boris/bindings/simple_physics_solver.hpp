#ifndef _SIMPLE_PHYSICS_SOLVER_HPP_
#define _SIMPLE_PHYSICS_SOLVER_HPP_

#include <cstdlib>
using namespace std;

namespace simple_physics_solver
{
  template<
    typename scalar,
    typename time
  >
  class SimplePhysicsSolverConfig
  {
    public:
      scalar omega_e;
      scalar omega_b;
      scalar epsilon;
      scalar sigma;
      scalar sigma2;
      scalar external_e_field_matrix[3][3];
      scalar b_field_matrix[3][3];

      SimplePhysicsSolverConfig(const scalar omega_e = scalar(4.9), const scalar omega_b = scalar(25.0),
                                const scalar epsilon = scalar(-1.0), const scalar sigma = scalar(0.1));
      virtual ~SimplePhysicsSolverConfig();
  };

  template<typename scalar, typename time>
  static void evaluate_external_e_field(const scalar** positions, const scalar* charges, const scalar* masses,
                                        const size_t num_particles, const size_t dim, const time t,
                                        const SimplePhysicsSolverConfig<scalar, time>& config,
                                        scalar** forces);

  template<typename scalar, typename time>
  static void evaluate_internal_e_field(const scalar** positions, const scalar* charges, const scalar* masses,
                                        const size_t num_particles, const size_t dim, const time t,
                                        const SimplePhysicsSolverConfig<scalar, time>& config,
                                        scalar** exyz, scalar* phis);

  template<typename scalar, typename time>
  static void get_b_field_vector(const SimplePhysicsSolverConfig<scalar, time>& config, const size_t dim,
                                 scalar* b_field_vector);

  template<typename scalar, typename time>
  static void evaluate_b_field(const scalar** velocities, const scalar* masses, const scalar* charges,
                               const size_t num_particles, const size_t dim, const time t,
                               const SimplePhysicsSolverConfig<scalar, time>& config,
                               scalar** forces);

  template<typename scalar, typename time>
  static scalar compute_energy(const scalar** positions, const scalar** velocities, const scalar* masses,
                               const scalar* charges,
                               const size_t num_particles, const size_t dim, const time t,
                               const SimplePhysicsSolverConfig<scalar, time>& config);

  namespace internal
  {
    template<typename scalar>
    static void cross_prod(const scalar* first, const scalar* second, scalar* cross_prod);

    template<typename scalar>
    static void scale_mat_mul_vec(const scalar** mat, const scalar* vec, const size_t rows, const size_t cols,
                                  const scalar factor, scalar* prod);
  }  // ::simple_physics_solver::internal
}  // ::simple_physics_solver

#endif
