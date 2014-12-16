#ifndef _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
#define _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>
#include <map>
using namespace std;

#include <Eigen/Core>

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

          virtual void log(el::base::type::ostream_t& os) const override {
            os << "pos: [" << x << " " << y << " " << z << "]\tvel: [" << y << " " << v << " " << w << "]";
          }
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


      template<typename scalar>
      static void init_config_options(po::options_description& opts)
      {
        opts.add_options()
          ("num_particles", po::value<size_t>(), "number of particles in the cloud")
          ("epsilon", po::value<scalar>(), "Boris' epsilon")
          ("omega_e", po::value<scalar>(), "E-field constant")
          ("omega_b", po::value<scalar>(), "B-field constant")
          ;
      }

      template<typename scalar>
      static void enable_config_options(size_t index = -1)
      {
        config::Options::get_instance()
          .register_init_function("Boris-SDC",
                                  std::function<void(po::options_description&)>(pfasst::examples::boris::init_config_options<scalar>),
                                  index);
      }


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

        protected:
          vector<shared_ptr<encap_type>> particles;
          vector<shared_ptr<encap_type>> previous_particles;
          vector<shared_ptr<encap_type>> tau_corrections;
          vector<acceleration_type> forces;
          vector<acceleration_type> previous_forces;
          vector<acceleration_type> b_vecs;
          vector<acceleration_type> previous_b_vecs;

          scalar initial_energy;
          vector<scalar> energy_evals;
          size_t f_evals;

          vector<velocity_type> s_integrals;
          vector<position_type> ss_integrals;

          //! delta_nodes[m] = nodes[m] - nodes[m-1]
          vector<time> delta_nodes;

          Matrix<time> s_mat;
          Matrix<time> ss_mat;
          Matrix<time> sx_mat;
          Matrix<time> st_mat;
          Matrix<time> q_mat;
          Matrix<time> qx_mat;
          Matrix<time> qt_mat;

          acceleration_type build_rhs(const size_t m, const bool previous = false) const
          {
//             VLOG_FUNC_START("BorisSweeper") << "m=" << m << ", previous=" << boolalpha << previous;
            acceleration_type rhs = (previous) ? this->previous_forces[m] : this->forces[m];
//             VLOG_INDENT(3) << "rhs:" << rhs;
            if (previous) {
              rhs += cross_prod(this->previous_particles[m]->velocities(), this->previous_b_vecs[m]);
//               VLOG_INDENT(3) << "  +=" << this->previous_particles[m]->velocities() << " x " << this->previous_b_vecs[m];
            } else {
              rhs += cross_prod(this->particles[m]->velocities(), this->b_vecs[m]);
//               VLOG_INDENT(3) << "  +=" << this->particles[m]->velocities() << " x " << this->b_vecs[m];
            }
//             VLOG_INDENT(3) << "  =>" << rhs;
//             VLOG_FUNC_END("BorisSweeper");
            return rhs;
          }

          scalar compute_residual()
          {
            const auto   nodes  = this->get_nodes();
            const size_t nnodes = nodes.size();
            time dt = this->get_controller()->get_time_step();

            scalar max_residual = 0.0;

            for (size_t m = 1; m < nnodes; ++m) {
              position_type pos;
              velocity_type vel;
              for (size_t j = 0; j < nnodes; ++j) {
                pos += this->q_mat(m, j) * dt * this->particles[j]->velocities();
                vel += this->q_mat(m, j) * dt * this->forces[j];
              }
              pos += this->particles[0]->positions() - this->particles[m]->positions();
              vel += this->particles[0]->velocities() - this->particles[m]->velocities();

              for (size_t j = 0; j < pos.size(); ++j) {
                for (size_t d = 0; d < pos[j].size(); ++d) {
                  max_residual = max(max_residual, abs(pos[j][d]));
                  max_residual = max(max_residual, abs(vel[j][d]));
                }
              }
            }

            return max_residual;
          }


        public:
          //! @{
          BorisSweeper(shared_ptr<bindings::WrapperInterface<scalar, time>>& impl_solver)
            :   impl_solver(impl_solver)
              , errors()
              , exact_updated(false)
              , f_evals(0)
          {
            VLOG_FUNC_START("BorisSweeper");
          }

          BorisSweeper(const BorisSweeper<scalar, time>& other) = delete;
          BorisSweeper(BorisSweeper<scalar, time>&& other) = delete;

          virtual ~BorisSweeper()
          {
            LOG(INFO) << "number force computations:" << this->f_evals;
            VLOG_FUNC_END("BorisSweeper");
          }
          //! @}

          //! @{
          virtual void set_state(shared_ptr<const Encapsulation<time>> u0, size_t m)
          {
            shared_ptr<const encap_type> u0_cast = dynamic_pointer_cast<const encap_type>(u0);
            assert(u0_cast);
            this->set_state(u0_cast, m);
          }

          virtual void set_state(shared_ptr<const encap_type> u0, size_t m)
          {
            this->particles[m]->operator=(*(u0.get()));
          }

          virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const override
          {
            return this->particles[m];
          }

          virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const override
          {
            return this->tau_corrections[m];
          }

          virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const override
          {
            return this->previous_particles[m];
          }

          virtual void set_initial_energy()
          {
            VLOG_FUNC_START("BorisSweeper");
            this->initial_energy = this->impl_solver->energy(this->particles.front(), this->get_controller()->get_time());
            VLOG_FUNC_END("BorisSweeper");
          }
          //! @}

          //! @{
          virtual void exact(shared_ptr<Encapsulation<time>> q, time t)
          {
            shared_ptr<encap_type> q_cast = dynamic_pointer_cast<encap_type>(q);
            assert(q_cast);
            this->exact(q_cast, t);
          }

          virtual void exact(shared_ptr<encap_type> q, const time t)
          {
            this->exact(*(q.get()), t);
          }

          virtual void exact(encap_type& q, const time t)
          {
            VLOG_FUNC_START("BorisSweeper") << "n=" << q.size() << ", t=" << t;
            UNUSED(t);
            if (!this->exact_updated) {
              typedef complex<scalar> C;
              C i(0.0, 1.0);
              auto initial = this->particles[0];
              scalar x0 = initial->positions()[0][0],
                     y0 = initial->positions()[0][1],
                     z0 = initial->positions()[0][2],
                     u0 = initial->velocities()[0][0],
                     v0 = initial->velocities()[0][1],
                     w0 = initial->velocities()[0][2],
                     omega_e = this->impl_solver->omega_e(),
                     omega_b = this->impl_solver->omega_b(),
                     epsilon = this->impl_solver->epsilon();
              time dt = this->get_controller()->get_time_step();

              C omega_tilde = sqrt(-2.0 * epsilon) * omega_e;
              q.positions()[0][2] = (z0 * cos(omega_tilde * (scalar)(dt)) 
                                    + w0 / omega_tilde * sin(omega_tilde * (scalar)(dt))).real();

              C sqrt_in_omega = sqrt(pow(omega_b, 2) + 4.0 * epsilon * pow(omega_e, 2));
              C omega_minus = 0.5 * (omega_b - sqrt_in_omega);
              C omega_plus = 0.5 * (omega_b + sqrt_in_omega);

              C r_minus = (omega_plus * x0 + v0) / (omega_plus - omega_minus);
              C r_plus = x0 - r_minus;

              C i_minus = (omega_plus * y0 - u0) / (omega_plus - omega_minus);
              C i_plus = y0 - i_minus;

              C x_y_move = (r_plus + i * i_plus) * exp(- i * omega_plus * (scalar)(dt))
                                                       + (r_minus + i * i_minus) * exp(- i * omega_minus * (scalar)(dt));
              q.positions()[0][0] = x_y_move.real();
              q.positions()[0][1] = x_y_move.imag();

              q.velocities()[0][2] = (- z0 * omega_tilde * sin(omega_tilde * (scalar)(dt)) + w0 * cos(omega_tilde * (scalar)(dt))).real();
              C u_v_move = (- i * omega_plus * (r_plus + i * i_plus)) * exp(-i * omega_plus * (scalar)(dt))
                                         - (i * omega_minus * (r_minus + i * i_minus)) * exp(-i * omega_minus * (scalar)(dt));
              q.velocities()[0][0] = u_v_move.real();
              q.velocities()[0][1] = u_v_move.imag();

              VLOG_INDENT(1) << "exact solution at time t=" << t << q;

              this->exact_cache = make_shared<encap_type>(q);
              this->exact_updated = true;
            } else {
              LOG(DEBUG) << "exact solution has been computed previously.";
              q = *(this->exact_cache.get());
              VLOG_INDENT(1) << "exact solution at time t=" << t << q;
            }
            VLOG_FUNC_END("BorisSweeper");
          }

          virtual void echo_error(const time t, bool predict = false)
          {
            if (this->particles[0]->size() == 1) {
              auto end = this->particles.back();
              shared_ptr<encap_type> ex = dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution));
              this->exact(ex, t);
              scalar e_end = this->impl_solver->energy(this->particles.back(), t);

              ErrorTuple<scalar> e_tuple;
              e_tuple.e_drift = this->initial_energy - e_end;
              e_tuple.p_err.x = ex->positions()[0][0] - end->positions()[0][0];
              e_tuple.p_err.y = ex->positions()[0][1] - end->positions()[0][1];
              e_tuple.p_err.z = ex->positions()[0][2] - end->positions()[0][2];
              e_tuple.p_err.u = ex->velocities()[0][0] - end->velocities()[0][0];
              e_tuple.p_err.v = ex->velocities()[0][1] - end->velocities()[0][1];
              e_tuple.p_err.w = ex->velocities()[0][2] - end->velocities()[0][2];
//               e_tuple.res = this->compute_residual();

              size_t n = this->get_controller()->get_step();
              size_t k = this->get_controller()->get_iteration();
              error_index nk(n, k);

              LOG(INFO) << scientific << setprecision(5)
                        << "step" << n+1
                        << "\titer" << k
                        << "\tres" << e_tuple.res
                        << "\tdrift" << e_tuple.e_drift
                        << "\tenergy" << e_end;
              VLOG_INDENT(1) << "particle at t_end:" << *end;
              VLOG_INDENT(2) << "absolute error:" << e_tuple.p_err;

              this->errors.insert(pair<error_index, ErrorTuple<scalar>>(nk, e_tuple));
            } else {
              // no-op
              // exact anlytical solution only valid for 1-particle-system
            }
          }

          virtual error_map<scalar> get_errors() const
          {
            return this->errors;
          }
          //! @}

          //! @{
          virtual void setup(bool coarse = false) override
          {
            VLOG_FUNC_START("BorisSweeper") << "coarse=" << coarse;
            auto nodes = this->get_nodes();
            assert(nodes.size() >= 1);
            const size_t nnodes = nodes.size();

            // compute delta nodes
            this->delta_nodes = vector<time>(nnodes, time(0.0));
            for (size_t m = 1; m < nnodes; m++) {
              this->delta_nodes[m] = nodes[m] - nodes[m - 1];
            }

            this->energy_evals.resize(nnodes);
            for (size_t m = 0; m < nnodes; ++m) {
              this->particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
              this->previous_particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
              this->forces.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
              this->previous_forces.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
              this->b_vecs.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
              this->previous_b_vecs.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
              this->s_integrals.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
              this->ss_integrals.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
              if (coarse) {
                this->tau_corrections.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
              }
            }

            this->q_mat = this->get_quadrature()->get_q_mat();
            auto qq_mat = this->q_mat * this->q_mat;
            this->s_mat = Matrix<time>(nnodes, nnodes);
            this->s_mat.fill(time(0.0));
            this->qt_mat = Matrix<time>(nnodes, nnodes);
            this->qt_mat.fill(time(0.0));
            this->qx_mat = Matrix<time>(nnodes, nnodes);
            this->qx_mat.fill(time(0.0));
            this->ss_mat = Matrix<time>(nnodes, nnodes);
            this->ss_mat.fill(time(0.0));
            this->sx_mat = Matrix<time>(nnodes, nnodes);
            this->sx_mat.fill(time(0.0));
            this->st_mat = Matrix<time>(nnodes, nnodes);
            this->st_mat.fill(time(0.0));

            // building rules for Q_E and Q_I:
            //  Q_E is striclty lower diagonal matrix with delta nodes of column index
            //  Q_I is lower diagonal matrix with first row and column all zero and delta nodes of 
            //      column index minus one
            Matrix<time> qe_mat = Matrix<time>(nnodes, nnodes);
            qe_mat.fill(time(0.0));
            Matrix<time> qi_mat = Matrix<time>(nnodes, nnodes);
            qi_mat.fill(time(0.0));
            for (size_t i = 0; i < nnodes; ++i) {
              for (size_t j = 0; j < nnodes; ++j) {
                if (j < i)           { qe_mat(i, j) = delta_nodes[j + 1]; }
                if (j > 0 && j <= i) { qi_mat(i, j) = delta_nodes[j]; }
              }
            }

            // Q_T = 0.5 * (Q_E + Q_I)
            qt_mat = 0.5 * (qe_mat + qi_mat);

            // Q_x = Q_E * Q_T + 0.5 * (Q_E ∘ Q_E)
            //  (only first term, i.e. matrix product)
            this->qx_mat = qe_mat * qt_mat;

            // compute QQ, S, SS, Q_T, Q_x
            for (size_t i = 1; i < nnodes; i++) {
              for (size_t j = 0; j < nnodes; j++) {
                // continuation of Q_x
                //  (i.e. Hadamard product of Q_E)
                this->qx_mat(i, j) += 0.5 * qe_mat(i, j) * qe_mat(i, j);
              }

              this->s_mat.row(i) = this->q_mat.row(i) - this->q_mat.row(i - 1);
              this->ss_mat.row(i) = qq_mat.row(i) - qq_mat.row(i - 1);
              this->sx_mat.row(i) = this->qx_mat.row(i) - this->qx_mat.row(i - 1);
              this->st_mat.row(i) = qt_mat.row(i) - qt_mat.row(i - 1);
            }

            UNUSED(coarse);
            VLOG_FUNC_END("BorisSweeper");
          }

          virtual void advance() override
          {
            VLOG_FUNC_START("BorisSweeper");
            this->set_state(const_pointer_cast<const encap_type>(this->particles.back()), 0);
            this->energy_evals.front() = this->energy_evals.back();
            this->forces.front() = this->forces.back();
            this->b_vecs.front() = this->b_vecs.back();
            this->exact_updated = false;
            VLOG_FUNC_END("BorisSweeper");
          }

          virtual void evaluate(size_t m)
          {
            VLOG_FUNC_START("BorisSweeper") << "node=" << m;
            time t = this->get_controller()->get_time() + this->get_controller()->get_time_step() * this->delta_nodes[m];
            VLOG_INDENT(2) << "particles:" << this->particles[m];
  //           Vector3d<scalar> vel; vel.fill(scalar(0.0));
  //           Vector3d<scalar> B; B.fill(scalar(0.0));
  //           Vector3d<scalar> b; b.fill(scalar(0.0));
  //           Eigen::Matrix<scalar, 3, 6> mat;
  //           mat << - this->epsilon * pow(this->get_e_field().omega_e, 2), 0, 0, 0, this->get_b_field().omega_b, 0,
  //                  0, - this->epsilon * pow(this->get_e_field().omega_e, 2), 0, - this->get_b_field().omega_b, 0, 0,
  //                  0, 0, 2 * this->epsilon * pow(this->get_e_field().omega_e, 2), 0, 0, 0;
  //           Eigen::Matrix<scalar, 6, 1> u_vec;
  //           u_vec.block(0, 0, 3, 1) = this->particles[m]->pos().as_matrix().transpose();
  //           u_vec.block(3, 0, 3, 1) = this->particles[m]->vel().as_matrix().transpose();

  //           Eigen::Matrix<scalar, 3, 1> f; f.fill(scalar(0.0));
  //           f = mat * u_vec;
  //           this->particles[m]->accel() = acceleration_type(f);
  //           auto e = this->get_e_field().evaluate(this->particles, m, t);
  //           vel = this->particles[m]->vel().as_matrix().transpose().array();
  //           B = this->get_b_field().get_field_vector().transpose().array() / this->particles[m]->alpha();
  //           b = dot(vel, B);
  //           this->particles[m]->accel() = this->particles[m]->alpha() * (e + acceleration_type(b));

            this->forces[m] = this->impl_solver->e_field_evaluate(this->particles[m], t);
            VLOG_INDENT(3) << "e_forces:" << this->forces[m];
            this->b_vecs[m] = this->impl_solver->b_field_vecs(this->particles[m], t);
            VLOG_INDENT(3) << "b_vecs:" << this->b_vecs[m];
            this->f_evals++;
            VLOG_FUNC_END("BorisSweeper");
          }

          virtual void predict(bool initial) override
          {
            VLOG_FUNC_START("BorisSweeper") << "initial=" << initial;
            UNUSED(initial);
            this->spread();
            for (size_t m = 0; m < this->particles.size(); ++m) {
              this->evaluate(m);
            }
            this->save();
            VLOG_FUNC_END("BorisSweeper");
          }

          virtual void sweep() override
          {
            VLOG_FUNC_START("BorisSweeper");
            const auto   nodes  = this->get_nodes();
            const size_t nnodes = nodes.size();
            assert(nnodes >= 1);
            time t  = this->get_controller()->get_time();
            time dt = this->get_controller()->get_time_step();

            this->impl_solver->energy(this->particles.front(), t);

            velocity_type c_k_term = cloud_component_factory<scalar>(this->particles[0]->size(), this->particles[0]->dim());

            // compute integrals
            for(size_t i = 0; i < nnodes; ++i) {
              zero(this->s_integrals[i]);
              zero(this->ss_integrals[i]);
            }
            // starting at m=1 as m=0 will only add zeros
            for (size_t m = 1; m < nnodes; m++) {
              for (size_t l = 0; l < nnodes; l++) {
                this->s_integrals[m] += this->build_rhs(l, true) * dt * this->s_mat(m, l);
                this->ss_integrals[m] += this->build_rhs(l, true) * dt * dt * this->ss_mat(m, l);
              }
            }

            this->evaluate(0);

            for (size_t m = 0; m < nnodes - 1; m++) {
              time ds = dt * this->delta_nodes[m+1];

              //// Update Position (explicit)
              //
              // x_{m+1}^{k+1} = x_{m}^{k+1}
              this->particles[m+1]->positions() = this->particles[m]->positions();
              //               + delta_nodes_{m+1} * v_{0}
              this->particles[m+1]->positions() += this->particles[0]->velocities() * ds;
              //               + \sum_{l=1}^{m} sx_{m+1,l}^{x} (f_{l}^{k+1} - f_{l}^{k})
              for (size_t l = 0; l <= m; l++) {
                this->particles[m+1]->positions() += this->build_rhs(l) * dt * dt * this->sx_mat(m+1, l);
                this->particles[m+1]->positions() -= this->build_rhs(l, true) * dt * dt * this->sx_mat(m+1, l);
              }

              this->particles[m+1]->positions() += this->ss_integrals[m+1];
              VLOG_INDENT(3) << "new positions:" << this->particles[m+1]->positions();

              // evaluate electric field with new position
              this->forces[m+1] = this->impl_solver->e_field_evaluate(this->particles[m+1], t + nodes[m]);

              //// Update Velocity (semi-implicit)
              zero(c_k_term);  // reset
              VLOG_INDENT(3) << "c_k_term:" << c_k_term;
              //                 - delta_nodes_{m} / 2 * f_{m+1}^{k}
              auto t1 = this->build_rhs(m+1, true);
              c_k_term -= (0.5 * t1 * ds);
//               VLOG_INDENT(3) << "  -= 0.5 *" << t1 << "*" << ds << "  =>" << c_k_term;
              //                 - delta_nodes_{m} / 2 * f_{m}^{k}
              auto t2 = this->build_rhs(m, true);
              c_k_term -= (0.5 * t2 * ds);
//               VLOG_INDENT(3) << "  -= 0.5 *" << t2 << "*" << ds << "  =>" << c_k_term;
              //                 + s_integral[m]
              c_k_term += this->s_integrals[m+1];
//               VLOG_INDENT(3) << "  +=" << this->s_integrals[m+1];
              VLOG_INDENT(3) << " ==>" << c_k_term;

              // doing Boris' magic
              this->boris_solve(nodes[m], nodes[m+1], ds, m, c_k_term);
              VLOG_INDENT(3) << "new velocities:" << this->particles[m+1]->velocities();

              // TODO XXX BUG
//               this->forces[m+1] += this->impl_solver->b_field_evaluate(this->particles[m+1], t);
            }
            this->save();
            this->echo_error(t + dt);
            VLOG_FUNC_END("BorisSweeper");
          }

          virtual void save(bool initial_only=false) override
          {
            VLOG_FUNC_START("BorisSweeper") << "initial_only=" << initial_only;
            if (initial_only) {
              this->previous_particles[0] = make_shared<encap_type>(*(this->particles[0].get()));
              this->previous_forces[0] = this->forces[0];
              this->previous_b_vecs[0] = this->b_vecs[0];
            } else {
              for (size_t m = 0; m < this->previous_particles.size(); m++) {
                this->previous_particles[m] = make_shared<encap_type>(*(this->particles[m].get()));
                this->previous_forces[m] = this->forces[m];
                this->previous_b_vecs[m] = this->b_vecs[m];
              }
            }
            VLOG_FUNC_END("BorisSweeper");
          }

          virtual void spread() override
          {
            VLOG_FUNC_START("BorisSweeper");
            for (size_t m = 1; m < this->particles.size(); ++m) {
              this->set_state(const_pointer_cast<const encap_type>(this->particles[0]), m);
            }
            VLOG_FUNC_END("BorisSweeper");
          }
          //! @}

          //! @{
          virtual void post(ICommunicator* comm, int tag) override
          {
            UNUSED(comm); UNUSED(tag);
            // TODO: implement BorisSweeper::post
          };

          virtual void send(ICommunicator* comm, int tag, bool blocking) override
          {
            UNUSED(comm); UNUSED(tag); UNUSED(blocking);
            // TODO: implement BorisSweeper::send
          }

          virtual void recv(ICommunicator* comm, int tag, bool blocking) override
          {
            UNUSED(comm); UNUSED(tag); UNUSED(blocking);
            // TODO: implement BorisSweeper::recv
          }

          virtual void broadcast(ICommunicator* comm) override
          {
            UNUSED(comm);
            // TODO: implement BorisSweeper::broadcast
          }
          //! @}

          //! @{
          virtual void boris_solve(const time tm, const time t_next, const time ds, const size_t m,
                                   const velocity_type& c_k_term)
          {
            VLOG_FUNC_START("BorisSweeper") << "tm=" << tm << ", t_next=" << t_next << ", node=" << m;
            UNUSED(t_next);
            velocity_type c_k_term_half = c_k_term / scalar(2.0);
            AttributeValues<scalar> beta = cmp_wise_div(this->particles[m]->charges(), this->particles[m]->masses()) / scalar(2.0) * ds;
            acceleration_type e_forces_mean = (this->forces[m] + this->forces[m+1]) / scalar(2.0);

            // first Boris' drift
            //   v^{-} = v^{k}
            velocity_type v_minus = this->particles[m]->velocities();
            //           + \beta * E_{mean} + c^{k} / 2
            v_minus += e_forces_mean * beta + c_k_term_half;

            // Boris' kick
            velocity_type boris_t(beta.size());
            auto b_field_vector = this->impl_solver->get_b_field_vector();
            for (size_t p = 0; p < beta.size(); ++p) {
              boris_t[p] = b_field_vector * beta[p];
            }
            velocity_type v_prime = v_minus + cross_prod(v_minus, boris_t);

            // final Boris' drift
            vector<scalar> boris_t_sqr(boris_t.size());
            for (size_t p = 0; p < boris_t.size(); ++p) {
              boris_t_sqr[p] = pow(norm0(boris_t[p]), 2);
            }
            velocity_type boris_s = (scalar(2.0) * boris_t) / (scalar(1.0) + boris_t_sqr);
            velocity_type v_plus = v_minus + cross_prod(v_prime, boris_s);

  //           assert(abs(v_minus.norm0() - v_plus.norm0()) <= 10e-8);

            this->particles[m+1]->velocities() = v_plus + e_forces_mean * beta + c_k_term_half;
            VLOG_FUNC_END("BorisSweeper");
          }
          //! @}
      };
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
