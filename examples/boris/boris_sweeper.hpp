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
#include "particle_3d.hpp"
#include "physics.hpp"

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
      struct ParticleError
      {
        scalar x;
        scalar y;
        scalar z;
        scalar u;
        scalar v;
        scalar w;
      };

      template<typename scalar>
      struct ErrorTuple
      {
        ParticleError<scalar> p_err;
        scalar err;
        scalar res;
      };

      template<typename scalar>
      using error_map = map<error_index, ErrorTuple<scalar>>;


      template<typename scalar>
      static Vector3d<scalar> dot(const Vector3d<scalar>& first, const Vector3d<scalar>& second)
      {
        Vector3d<scalar> vec;
        vec[0] = first[1] * second[2] - first[2] * second[1];
        vec[1] = first[2] * second[0] - first[0] * second[2];
        vec[2] = first[0] * second[1] - first[1] * second[0];
        return vec;
      }

      template<typename scalar>
      static scalar norm(const Vector3d<scalar>& vec)
      {
        return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
      }


      template<typename scalar>
      static void init_config_options(po::options_description& opts)
      {
        opts.add_options()
  //         ("num_particles", po::value<size_t>(), "number of particles in the cloud")
          ("epsilon", po::value<scalar>(), "Boris' epsilon")
          ("omega_e", po::value<scalar>(), "E-field constant")
          ("omega_b", po::value<scalar>(), "B-field constant")
          ;
      }

      template<typename scalar>
      static void enable_config_options(size_t index = -1)
      {
        pfasst::config::Options::get_instance()
          .register_init_function("Boris-SDC",
                                  std::function<void(po::options_description&)>(pfasst::examples::boris::init_config_options<scalar>),
                                  index);
      }


      template<
        typename scalar,
        typename time,
        typename EnergyOperatorT
      >
      class BorisSweeper
        : public EncapSweeper<time>
      {
        public:
          typedef Particle3DEncapsulation<scalar, time> encap_type;
          typedef typename encap_type::position_type position_type;
          typedef typename encap_type::velocity_type velocity_type;
          typedef typename encap_type::acceleration_type acceleration_type;
          typedef Particle3DFactory<scalar, time> factory_type;
          typedef EnergyOperatorT energy_operator_type;
          typedef typename energy_operator_type::e_field_type e_field_type;
          typedef typename energy_operator_type::b_field_type b_field_type;

        private:
          energy_operator_type energy_op;
          scalar epsilon;
          error_map<scalar> errors;

        protected:
          vector<shared_ptr<encap_type>> particles;
          vector<shared_ptr<encap_type>> previous_particles;
          vector<shared_ptr<encap_type>> tau_corrections;

          scalar initial_energy;
          vector<acceleration_type> energy_evals;
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
                pos += this->q_mat(m, j) * this->particles[j]->vel().convert(boris::dt<time>(dt));
                vel += this->q_mat(m, j) * this->particles[j]->accel().convert(boris::dt<time>(dt));
              }
              pos += this->particles[0]->pos() - this->particles[m]->pos();
              vel += this->particles[0]->vel() - this->particles[m]->vel();

              for (size_t j = 0; j < pos.DIM; ++j) {
                max_residual = max(max_residual, abs(pos[j]));
                max_residual = max(max_residual, abs(vel[j]));
              }
            }

            return max_residual;
          }


        public:
          //! @{
          BorisSweeper(const energy_operator_type& energy_operator, scalar epsilon = -1.0)
            :   energy_op(energy_operator)
              , epsilon(epsilon)
              , errors()
              , f_evals(0)
          {}

          BorisSweeper(const BorisSweeper<scalar, time, EnergyOperatorT>& other) = delete;
          BorisSweeper(BorisSweeper<scalar, time, EnergyOperatorT>&& other) = delete;

          virtual ~BorisSweeper()
          {
            LOG(INFO) << "f evals:" << this->f_evals;
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
            this->initial_energy = this->energy_op.evaluate({this->particles.front()}, this->get_controller()->get_time());
          }
          //! @}

          //! @{
          virtual const e_field_type& get_e_field() const
          {
            return this->energy_op.get_e_field();
          }

          virtual const b_field_type& get_b_field() const
          {
            return this->energy_op.get_b_field();
          }

          virtual void set_energy_operator(const energy_operator_type& e_operator)
          {
            this->energy_op = e_operator;
          }

          virtual const energy_operator_type& get_energy_operator() const
          {
            return this->energy_op;
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

          //! FIXME: there is probably still a bug in the calculation of the exact solution here
          virtual void exact(encap_type& q, const time t)
          {
            typedef complex<scalar> C;
            C i(0.0, 1.0);
            auto initial = this->particles[0];
            scalar x0 = initial->pos().x,
                   y0 = initial->pos().y,
                   z0 = initial->pos().z,
                   u0 = initial->vel().u,
                   v0 = initial->vel().v,
                   w0 = initial->vel().w,
                   omega_e = this->get_e_field().omega_e,
                   omega_b = this->get_b_field().omega_b;
            C c_epsilon = C(this->epsilon, 0.0);

            C omega_tilde = sqrt(-2.0 * c_epsilon) * omega_e;
            q.pos().z = (z0 * cos(omega_tilde * (scalar)(t)) + w0 / omega_tilde * sin(omega_tilde * (scalar)(t))).real();

            C sqrt_in_omega = sqrt(pow(omega_b, 2) + 4.0 * epsilon * pow(omega_e, 2));
            C omega_minus = 0.5 * (omega_b - sqrt_in_omega);
            C omega_plus = 0.5 * (omega_b + sqrt_in_omega);

            C r_minus = (omega_plus * x0 + v0) / (omega_plus - omega_minus);
            C r_plus = x0 - r_minus;

            C i_minus = (omega_plus * y0 - u0) / (omega_plus - omega_minus);
            C i_plus = y0 - i_minus;

            C x_y_move = (r_plus + i * i_plus) * exp(- i * omega_plus * (scalar)(t))
                                                     + (r_minus + i * i_minus) * exp(- i * omega_minus * (scalar)(t));
            q.pos().x = x_y_move.real();
            q.pos().y = x_y_move.imag();

            q.vel().w = (- z0 * omega_tilde * sin(omega_tilde * (scalar)(t)) + w0 * cos(omega_tilde * (scalar)(t))).real();
            C u_v_move = (- i * omega_plus * (r_plus + i * i_plus)) * exp(-i * omega_plus * (scalar)(t))
                                       - (i * omega_minus * (r_minus + i * i_minus)) * exp(-i * omega_minus * (scalar)(t));
            q.vel().u = u_v_move.real();
            q.vel().v = u_v_move.imag();
          }

          virtual void echo_error(const time t, bool predict = false)
          {
            auto end = this->particles.back();
            shared_ptr<encap_type> ex = make_shared<encap_type>();
            this->exact(ex, t);
            scalar e_end = this->energy_op.evaluate({this->particles.back()}, t);

            ErrorTuple<scalar> e_tuple;
            e_tuple.err = this->initial_energy - e_end;
            e_tuple.p_err.x = ex->pos().x - end->pos().x;
            e_tuple.p_err.y = ex->pos().y - end->pos().y;
            e_tuple.p_err.z = ex->pos().z - end->pos().z;
            e_tuple.p_err.u = ex->vel().u - end->vel().u;
            e_tuple.p_err.v = ex->vel().v - end->vel().v;
            e_tuple.p_err.w = ex->vel().w - end->vel().w;
            e_tuple.res = this->compute_residual();

            size_t n = this->get_controller()->get_step();
            size_t k = this->get_controller()->get_iteration();
            error_index nk(n, k);

            LOG(INFO) << "err:" << n << k
                      << "\tdrift:" << e_tuple.err
                      << "\tres:" << e_tuple.res
                      << "\t(energy:" << e_end
                      << "\tpos:" << end->pos()
                      << "\tvel:" << end->vel()
                      << "\t" << predict << ")";

            this->errors.insert(pair<error_index, ErrorTuple<scalar>>(nk, e_tuple));
          }

          virtual error_map<scalar> get_errors() const
          {
            return this->errors;
          }
          //! @}

          //! @{
          virtual void setup(bool coarse = false) override
          {
            auto nodes = this->get_nodes();
            assert(nodes.size() >= 1);
            const size_t nnodes = nodes.size();

            // compute delta nodes
            this->delta_nodes = vector<time>(nnodes, time(0.0));
            for (size_t m = 1; m < nnodes; m++) {
              this->delta_nodes[m] = nodes[m] - nodes[m - 1];
            }

            for (size_t m = 0; m < nnodes; ++m) {
              this->particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
              this->previous_particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
              this->s_integrals.push_back(velocity_type());
              this->ss_integrals.push_back(position_type());
              this->energy_evals.push_back(acceleration_type());
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
          }

          virtual void advance() override
          {
            this->set_state(const_pointer_cast<const encap_type>(this->particles.back()), 0);
            this->energy_evals.front() = this->energy_evals.back();
          }

          void evaluate(size_t m)
          {
      //       Vector3d<scalar> vel; vel.fill(scalar(0.0));
      //       Vector3d<scalar> B; B.fill(scalar(0.0));
      //       Vector3d<scalar> b; b.fill(scalar(0.0));
            Eigen::Matrix<scalar, 3, 6> mat;
            mat << - this->epsilon * pow(this->get_e_field().omega_e, 2), 0, 0, 0, this->get_b_field().omega_b, 0,
                   0, - this->epsilon * pow(this->get_e_field().omega_e, 2), 0, - this->get_b_field().omega_b, 0, 0,
                   0, 0, 2 * this->epsilon * pow(this->get_e_field().omega_e, 2), 0, 0, 0;
            Eigen::Matrix<scalar, 6, 1> u_vec;
            u_vec.block(0, 0, 3, 1) = this->particles[m]->pos().as_matrix().transpose();
            u_vec.block(3, 0, 3, 1) = this->particles[m]->vel().as_matrix().transpose();

            Eigen::Matrix<scalar, 3, 1> f; f.fill(scalar(0.0));
            f = mat * u_vec;
            this->particles[m]->accel() = acceleration_type(f);
            time t = this->get_controller()->get_time() + this->get_controller()->get_time_step() * this->delta_nodes[m];
      //       auto e = this->get_e_field().evaluate(this->particles, m, t);
      //       vel = this->particles[m]->vel().as_matrix().transpose().array();
      //       B = this->get_b_field().get_field_vector().transpose().array() / this->particles[m]->alpha();
      //       b = dot(vel, B);
      //       this->particles[m]->accel() = this->particles[m]->alpha() * (e + acceleration_type(b));
            this->f_evals++;
          }

          virtual void reevaluate(bool initial) override
          {
            for (size_t m = 0; m < this->get_quadrature()->get_num_nodes(); m++) {
              this->evaluate(m);
            }
          }

          virtual void predict(bool initial) override
          {
            UNUSED(initial);
            this->spread();
            for (size_t m = 0; m < this->particles.size(); ++m) {
              this->evaluate(m);
            }
            this->save();
          }

          virtual void sweep() override
          {
            const auto   nodes  = this->get_nodes();
            const size_t nnodes = nodes.size();
            assert(nnodes >= 1);
            time t  = this->get_controller()->get_time();
            time dt = this->get_controller()->get_time_step();
            {
              scalar e = this->energy_op.evaluate({this->particles.front()}, t);
            }
            velocity_type c_k_term;

            // compute integrals
            for(size_t i = 0; i < nnodes; ++i) {
              this->s_integrals[i].zero();
              this->ss_integrals[i].zero();
            }
            // starting at m=1 as m=0 will only add zeros
            for (size_t m = 1; m < nnodes; m++) {
              for (size_t l = 0; l < nnodes; l++) {
                this->s_integrals[m] += this->previous_particles[l]
                                            ->accel().convert(boris::dt<time>(dt * this->s_mat(m, l)));
                this->ss_integrals[m] += this->previous_particles[l]
                                             ->accel().convert(boris::dtdt<time>(dt * dt * this->ss_mat(m, l)));
              }
            }

            this->evaluate(0);

            for (size_t m = 0; m < nnodes - 1; m++) {
              boris::dt<time> ds(dt * this->delta_nodes[m+1]);

              //// Update Position (explicit)
              //
              // x_{m+1}^{k+1} = x_{m}^{k}
              this->particles[m+1]->pos() = position_type(this->particles[m]->pos());
              //               + delta_nodes_{m+1} * v_{0}
              this->particles[m+1]->pos() += this->particles[0]->vel().convert(ds);
              //               + \sum_{l=1}^{m} sx_{m+1,l}^{x} (f_{l}^{k+1} - f_{l}^{k})
              for (size_t l = 0; l <= m; l++) {
                this->particles[m+1]->pos() += this->particles[l]
                                                   ->accel().convert(boris::dtdt<time>(dt * dt * this->sx_mat(m+1, l)));
                this->particles[m+1]->pos() -= this->previous_particles[l]
                                                   ->accel().convert(boris::dtdt<time>(dt * dt * this->sx_mat(m+1, l)));
              }

              this->particles[m+1]->pos() += this->ss_integrals[m+1];

              // evaluate electric field with new position
              this->particles[m+1]->accel() = this->get_e_field().evaluate(this->particles, m+1, t + nodes[m]);

              //// Update Velocity (semi-implicit)
              c_k_term.zero();  // reset
              //                 - delta_nodes_{m} / 2 * f_{m+1}^{k}
              c_k_term -= 0.5 * this->previous_particles[m+1]->accel().convert(ds);
              //                 - delta_nodes_{m} / 2 * f_{m}^{k}
              c_k_term -= 0.5 * this->previous_particles[m]->accel().convert(ds);
              //                 + s_integral[m]
              c_k_term += this->s_integrals[m+1];

              // doing Boris' magic
              this->boris_solve(nodes[m], nodes[m+1], ds, m, c_k_term);

              this->particles[m+1]->accel() += this->get_b_field().evaluate(this->particles, m+1, t);
            }
            this->save();
            this->echo_error(t + dt);
          }

          virtual void save(bool initial_only=false) override
          {
            if (initial_only) {
              this->previous_particles[0] = make_shared<encap_type>(*(this->particles[0].get()));
            } else {
              for (size_t m = 0; m < this->previous_particles.size(); m++) {
                this->previous_particles[m] = make_shared<encap_type>(*(this->particles[m].get()));
              }
            }
          }

          virtual void spread() override
          {
            for (size_t m = 1; m < this->particles.size(); ++m) {
              this->set_state(const_pointer_cast<const encap_type>(this->particles[0]), m);
            }
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
          virtual void boris_solve(const time tm, const time t_next, const boris::dt<time> ds, const size_t m,
                                   const velocity_type& c_k_term)
          {
            UNUSED(t_next);
            velocity_type c_k_term_half = c_k_term / scalar(2.0);
            boris::dt<time> beta(this->particles[m]->alpha() / time(2.0) * ds.v);
            acceleration_type e_mean = (this->get_e_field().evaluate(this->particles, m, tm) + this->particles[m+1]->accel()) / scalar(2.0);

            // first Boris' drift
            //   v^{-} = v^{k}
            velocity_type v_minus = this->particles[m]->vel();
            //           + \beta * E_{mean} + c^{k} / 2
            v_minus += e_mean.convert(beta) + c_k_term_half;

              // Boris' kick
            velocity_type boris_t(beta.v * this->get_b_field().get_field_vector());
            velocity_type v_prime = v_minus + v_minus.cross(boris_t);

              // final Boris' drift
            scalar boris_t_sqr = pow(boris_t.norm0(), 2);
            velocity_type boris_s = (scalar(2.0) * boris_t) / (scalar(1.0) + boris_t_sqr);
            velocity_type v_plus = v_minus + v_prime.cross(boris_s);

      //       assert(abs(v_minus.norm0() - v_plus.norm0()) <= 10e-8);

            this->particles[m+1]->vel() = v_plus + e_mean.convert(beta) + c_k_term_half;
          }
          //! @}
      };

    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
