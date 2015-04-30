/**
 * @defgroup BorisFiles Files
 * @ingroup Boris
 * @file examples/boris/boris_sweeper.hpp
 */
#ifndef _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
#define _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <complex>
#include <cstdlib>
#include <ios>
#include <iostream>
#include <map>
#include <vector>
using namespace std;

#include <Eigen/Core>

#include <boost/format.hpp>

#include <pfasst.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/encap/encap_sweeper.hpp>

#include "particle.hpp"
#include "particle_cloud.hpp"
#include "particle_util.hpp"
#include "bindings/wrapper_interface.hpp"


namespace pfasst
{
  namespace examples
  {
    /**
     * @defgroup Boris Boris
     * @ingroup Examples
     *
     * This directory contains an implementations of a modification to the Boris method to solve 
     * second order ODEs with the Velocity-Verlet scheme using the PFASST framework.
     *
     * The sweeper with it's _Boris magic_ is implemented in `boris_sweeper.hpp`.
     * The physical properties of a testbed example with a panning trap are defined in `physics.hpp`
     * and `simple_physics.hpp`, while the data structure for the particles are defined in
     * `particle.hpp` and `particle_3d.hpp`.
     */
    namespace boris
    {
      using namespace pfasst::encap;

      template<typename coeff>
      using Vector3d = Eigen::Array<coeff, 3, 1>;

      template<typename coeff>
      using Matrix3d = Eigen::Matrix<coeff, 3, 3, Eigen::RowMajor>;

      template<typename coeff>
      using Matrix = Eigen::Matrix<coeff, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

      typedef pair<size_t, size_t> error_index;

      template<typename scalar>
      class ParticleError :
        public el::Loggable
      {
        public:
          scalar x;
          scalar y;
          scalar z;
          scalar u;
          scalar v;
          scalar w;

          virtual void log(el::base::type::ostream_t& os) const override;
      };

      template<typename scalar>
      struct ErrorTuple
      {
        ParticleError<scalar> p_err;
        scalar e_drift;
        scalar res;
      };

      template<typename scalar>
      using error_map = map<error_index, ErrorTuple<scalar>>;


      template<typename precision = pfasst::time_precision>
      static void init_opts();

      template<typename precision = pfasst::time_precision>
      static void init_logs();


      class LogIndent
      {
        private:
          array<size_t, 9> vlog_levels;

        public:
          LogIndent();
          virtual ~LogIndent();

          void increment(const size_t vlevel);
          void decrement(const size_t vlevel);
          const string indent(const size_t vlevel) const;
      };


      /**
       * @ingroup Boris
       */
      template<
        typename scalar,
        typename time
      >
      class BorisSweeper
        : public EncapSweeper<time>
      {
        public:
          typedef ParticleCloud<scalar> encap_type;
          typedef ParticleCloudComponent<scalar> position_type;
          typedef ParticleCloudComponent<scalar> velocity_type;
          typedef ParticleCloudComponent<scalar> acceleration_type;

        private:
          shared_ptr<bindings::WrapperInterface<scalar, time>> impl_solver;
          error_map<scalar> errors;
          bool exact_updated;
          shared_ptr<encap_type> exact_cache;
          shared_ptr<LogIndent> log_indent;

        protected:
          vector<shared_ptr<encap_type>> particles;
          vector<shared_ptr<encap_type>> saved_particles;
          shared_ptr<encap_type> start_particles;
          shared_ptr<encap_type> end_particles;

          vector<shared_ptr<acceleration_type>> tau_q_corrections;
          vector<shared_ptr<acceleration_type>> tau_qq_corrections;
          vector<acceleration_type> forces;
          vector<acceleration_type> saved_forces;
          vector<acceleration_type> b_vecs;
          vector<acceleration_type> saved_b_vecs;

          scalar initial_energy;
          vector<scalar> energy_evals;
          size_t f_evals;

          bool coarse;

          vector<velocity_type> s_integrals;
          vector<position_type> ss_integrals;

          //! delta_nodes[m] = nodes[m] - nodes[m-1]
          vector<time> delta_nodes;

          Matrix<time> s_mat;
          Matrix<time> ss_mat;
          Matrix<time> sx_mat;
          Matrix<time> st_mat;
          Matrix<time> q_mat;
          Matrix<time> qq_mat;
          Matrix<time> qx_mat;
          Matrix<time> qt_mat;

          ofstream data_stream;

          boost::format log_fmt;
          string DATA_STREAM_FORMAT_STR;
          boost::format data_stream_fmt;

          acceleration_type build_rhs(const size_t m, const bool previous = false) const;
          scalar compute_residual_max();
          void write_center_to_file(const size_t iter, const size_t sweep, const ParticleComponent<scalar>& center,
                                    const scalar energy, const scalar drift, const scalar residual);
          void write_particle_cloud_to_file(const size_t iter, const size_t sweep, const shared_ptr<encap_type>& cloud,
                                            const scalar energy, const scalar drift, const scalar residual,
                                            const bool with_center = true);
          void update_position(const size_t m, const time dt, const time ds);
          void update_velocity(const size_t m, const time ds, const vector<time>& nodes);

        public:
          //! @{
          BorisSweeper(shared_ptr<bindings::WrapperInterface<scalar, time>>& impl_solver,
                       const string& data_file);
          BorisSweeper(const BorisSweeper<scalar, time>& other) = delete;
          BorisSweeper(BorisSweeper<scalar, time>&& other) = delete;

          virtual ~BorisSweeper();
          //! @}

          //! @{
          virtual void set_state(shared_ptr<const encap_type> u0, size_t m);
          virtual void set_start_state(shared_ptr<const encap_type> u0);
          virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const override;
          virtual shared_ptr<Encapsulation<time>> get_start_state() const;
          virtual shared_ptr<acceleration_type> get_tau_q_as_force(size_t m) const;
          virtual shared_ptr<acceleration_type> get_tau_qq_as_force(size_t m) const;
          virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const override;
          virtual void set_initial_energy();
          //! @}

          //! @{
          virtual void exact(shared_ptr<Encapsulation<time>> q, time t);
          virtual void exact(shared_ptr<encap_type> q, const time t);
          virtual void exact(encap_type& q, const time t);
          virtual void echo_error(const time t, bool predict = false);
          virtual error_map<scalar> get_errors() const;
          //! @}

          //! @{
          virtual void setup(bool coarse = false) override;
          virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override;
          virtual void integrate(time dt, vector<shared_ptr<acceleration_type>> dst_q,
                                 vector<shared_ptr<acceleration_type>> dst_qq) const;
          virtual void residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const override;
          virtual void advance() override;
          virtual void evaluate(size_t m);
          virtual void predict(bool initial) override;
          virtual void sweep() override;
          virtual void save(bool initial_only=false) override;
          virtual void spread() override;
          //! @}

          //! @{
          virtual void post_sweep() override;
          virtual void post_predict() override;
          virtual void post_step() override;
          //! @}

          //! @{
          virtual void post(ICommunicator* comm, int tag) override;
          virtual void send(ICommunicator* comm, int tag, bool blocking) override;
          virtual void recv(ICommunicator* comm, int tag, bool blocking) override;
          virtual void broadcast(ICommunicator* comm) override;
          //! @}

          //! @{
          virtual void boris_solve(const time tm, const time t_next, const time ds, const size_t m,
                                   const velocity_type& c_k_term);
          //! @}
      };
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst


#include "boris_sweeper_impl.hpp"

#endif  // _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
