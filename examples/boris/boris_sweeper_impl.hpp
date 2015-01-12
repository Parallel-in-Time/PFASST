#include "boris_sweeper.hpp"

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
    namespace boris
    {
      template<typename scalar>
      void ParticleError<scalar>::log(el::base::type::ostream_t& os) const
      {
        os << "pos: [" << x << " " << y << " " << z << "]\tvel: [" << y << " " << v << " " << w << "]";
      }


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
      static void enable_config_options(size_t index)
      {
        config::Options::get_instance()
          .register_init_function("Boris-SDC",
                                  std::function<void(po::options_description&)>(pfasst::examples::boris::init_config_options<scalar>),
                                  index);
      }


      template<typename scalar, typename time>
      typename BorisSweeper<scalar, time>::acceleration_type
      BorisSweeper<scalar, time>::build_rhs(const size_t m, const bool previous) const
      {
        VLOG_FUNC_START("BorisSweeper") << " m=" << m << ", previous=" << boolalpha << previous;
        VLOG(6) << LOG_INDENT << "building rhs for node " << m
                              << ((previous) ? " of previous sweep" : " of current sweep");
        acceleration_type rhs = (previous) ? this->saved_forces[m] : this->forces[m];
        if (previous) {
          rhs += cross_prod(this->saved_particles[m]->velocities(), this->saved_b_vecs[m]);
        } else {
          rhs += cross_prod(this->particles[m]->velocities(), this->b_vecs[m]);
        }
        VLOG(6) << LOG_INDENT << "  => " << rhs;
        VLOG_FUNC_END("BorisSweeper");
        return rhs;
      }

      template<typename scalar, typename time>
      scalar BorisSweeper<scalar, time>::compute_residual()
      {
        VLOG_FUNC_START("BorisSweeper");
        VLOG(6) << LOG_INDENT << "computing residual";
        const auto   nodes  = this->get_nodes();
        const size_t nnodes = nodes.size();
        time dt = this->get_controller()->get_time_step();

        scalar max_residual = 0.0;

        for (size_t m = 1; m < nnodes; ++m) {
          position_type pos = cloud_component_factory<scalar>(this->particles[0]->size(), this->particles[0]->dim());
          velocity_type vel = cloud_component_factory<scalar>(this->particles[0]->size(), this->particles[0]->dim());
          for (size_t j = 0; j < nnodes; ++j) {
            pos += this->q_mat(m, j) * dt * this->particles[j]->velocities();
            vel += this->q_mat(m, j) * dt * this->build_rhs(j);
          }
          pos += this->start_particles->positions() - this->particles[m]->positions();
          vel += this->start_particles->velocities() - this->particles[m]->velocities();

          for (size_t j = 0; j < pos.size(); ++j) {
            for (size_t d = 0; d < pos[j].size(); ++d) {
              max_residual = std::max(max_residual, abs(pos[j][d]));
              max_residual = std::max(max_residual, abs(vel[j][d]));
            }
          }
        }

        VLOG(6) << LOG_INDENT << "residual: " << max_residual;
        VLOG_FUNC_END("BorisSweeper");
        return max_residual;
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::write_center_to_file(const size_t iter, const size_t sweep,
                                                            const ParticleComponent<scalar>& center,
                                                            const scalar energy, const scalar drift, const scalar residual)
      {
        VLOG(5) << "Formatting string: '" << this->DATA_STREAM_FORMAT_STR << "'";
        this->data_stream_fmt % (iter+1) % sweep % -1
                              % center[0] % center[1] % center[2]
                              % 0 % 0 % 0
                              % energy % drift % residual;
        this->data_stream << this->data_stream_fmt << endl;
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::write_particle_cloud_to_file(const size_t iter, const size_t sweep,
                                                                    const shared_ptr<encap_type>& cloud,
                                                                    const scalar energy, const scalar drift,
                                                                    const scalar residual,
                                                                    const bool with_center)
      {
        VLOG(5) << "Formatting string: '" << this->DATA_STREAM_FORMAT_STR << "'";
        for (size_t p = 0; p < cloud->size(); ++p) {
          this->data_stream_fmt % (iter+1) % sweep % p
                                % cloud->positions()[p][0] % cloud->positions()[p][1] % cloud->positions()[p][2]
                                % cloud->velocities()[p][0] % cloud->velocities()[p][1] % cloud->velocities()[p][2]
                                % energy % drift % residual;
          this->data_stream << this->data_stream_fmt << endl;
        }
        if (with_center) {
          this->write_center_to_file(iter, sweep, cloud->center_of_mass(), energy, drift, residual);
        }
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::update_position(const size_t m, const time dt, const time ds)
      {
        this->particles[m+1]->positions() = this->particles[m]->positions();
        VLOG(4) << LOG_INDENT << "pos = " << this->particles[m+1]->positions();
        //               + delta_nodes_{m+1} * v_{0}
        this->particles[m+1]->positions() += this->start_particles->velocities() * ds;
        VLOG(5) << LOG_INDENT << "   += " << this->start_particles->velocities() << " * " << ds;
        //               + \sum_{l=1}^{m} sx_{m+1,l}^{x} (f_{l}^{k+1} - f_{l}^{k})
        for (size_t l = 0; l <= m; l++) {
          auto rhs_this = this->build_rhs(l);
          auto rhs_prev = this->build_rhs(l, true);
          this->particles[m+1]->positions() += rhs_this * dt * dt * this->sx_mat(m+1, l);
          VLOG(5) << LOG_INDENT << "   += " << rhs_this << " * " << dt * dt << " * " << this->sx_mat(m+1, l);
          this->particles[m+1]->positions() -= rhs_prev * dt * dt * this->sx_mat(m+1, l);
          VLOG(5) << LOG_INDENT << "   -= " << rhs_prev << " * " << dt * dt << " * " << this->sx_mat(m+1, l);
        }

        this->particles[m+1]->positions() += this->ss_integrals[m+1];
        VLOG(5) << LOG_INDENT << "   += " << this->ss_integrals[m+1];
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::update_velocity(const size_t m, const time ds, const vector<time> nodes)
      {
        velocity_type c_k_term = cloud_component_factory<scalar>(this->particles[0]->size(), this->particles[0]->dim());
        zero(c_k_term);  // TODO: check if this is required

        VLOG(4) << LOG_INDENT << "c_k_term: " << c_k_term;
        //                 - delta_nodes_{m} / 2 * f_{m+1}^{k}
        auto t1 = this->build_rhs(m+1, true);
        c_k_term -= (0.5 * t1 * ds);
        VLOG(5) << LOG_INDENT << "  -= 0.5 * " << t1 << " * " << ds << "  => " << c_k_term;
        //                 - delta_nodes_{m} / 2 * f_{m}^{k}
        auto t2 = this->build_rhs(m, true);
        c_k_term -= (0.5 * t2 * ds);
        VLOG(5) << LOG_INDENT << "  -= 0.5 * " << t2 << " * " << ds << "  => " << c_k_term;
        //                 + s_integral[m]
        c_k_term += this->s_integrals[m+1];
        VLOG(5) << LOG_INDENT << "  += " << this->s_integrals[m+1];
        VLOG(4) << LOG_INDENT << " ==> " << c_k_term;

        // doing Boris' magic
        this->boris_solve(nodes[m], nodes[m+1], ds, m, c_k_term);
      }


      template<typename scalar, typename time>
      BorisSweeper<scalar, time>::BorisSweeper(shared_ptr<bindings::WrapperInterface<scalar, time>>& impl_solver,
                                               const string& data_file)
        :   impl_solver(impl_solver)
          , errors()
          , exact_updated(false)
          , f_evals(0)
          , data_stream(data_file, ios_base::out | ios_base::trunc)
      {
        VLOG_FUNC_START("BorisSweeper");
        assert(data_stream.is_open() && data_stream.good());
        LOG(INFO) << "writing particle data to: " << data_file;
        // CSV format specification:  [step],[iter],[particle],[x],[y],[z],[u],[v],[w],[energy],[drift],[residual]
        this->DATA_STREAM_FORMAT_STR = "%d,%d,%d,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f";
        //                              i  sw p  x    y    z    u    v    w    engy edr  res
        this->data_stream_fmt = boost::format(this->DATA_STREAM_FORMAT_STR);
      }


      template<typename scalar, typename time>
      BorisSweeper<scalar, time>::~BorisSweeper()
      {
        LOG(INFO) << "number force computations:" << this->f_evals;
        VLOG_FUNC_END("BorisSweeper");
      }


      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::set_state(shared_ptr<const Encapsulation<time>> u0, size_t m)
      {
        shared_ptr<const encap_type> u0_cast = dynamic_pointer_cast<const encap_type>(u0);
        assert(u0_cast);
        this->set_state(u0_cast, m);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::set_state(shared_ptr<const encap_type> u0, size_t m)
      {
        this->particles[m]->operator=(*(u0.get()));
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::set_start_state(shared_ptr<const encap_type> u0)
      {
        this->start_particles->operator=(*(u0.get()));
      }

      template<typename scalar, typename time>
      shared_ptr<Encapsulation<time>> BorisSweeper<scalar, time>::get_state(size_t m) const
      {
        return this->particles[m];
      }

      template<typename scalar, typename time>
      shared_ptr<typename BorisSweeper<scalar, time>::encap_type> BorisSweeper<scalar, time>::get_start_state() const
      {
        return this->start_particles;
      }

      template<typename scalar, typename time>
      shared_ptr<Encapsulation<time>> BorisSweeper<scalar, time>::get_tau(size_t m) const
      {
        return this->tau_corrections[m];
      }

      template<typename scalar, typename time>
      shared_ptr<Encapsulation<time>> BorisSweeper<scalar, time>::get_saved_state(size_t m) const
      {
        return this->saved_particles[m];
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::set_initial_energy()
      {
        VLOG_FUNC_START("BorisSweeper");
        VLOG(6) << LOG_INDENT << "computing and setting initial energy";

        auto p0 = this->start_particles;
        VLOG(7) << LOG_INDENT << "initial particles: " << p0;
        this->initial_energy = this->impl_solver->energy(p0, this->get_controller()->get_time());
        LOG(INFO) << OUT::green << "initial total energy: " << this->initial_energy;

        this->write_particle_cloud_to_file(0, 0, p0, this->initial_energy, 0, 0);
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::exact(shared_ptr<Encapsulation<time>> q, time t)
      {
        shared_ptr<encap_type> q_cast = dynamic_pointer_cast<encap_type>(q);
        assert(q_cast);
        this->exact(q_cast, t);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::exact(shared_ptr<encap_type> q, const time t)
      {
        this->exact(*(q.get()), t);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::exact(encap_type& q, const time t)
      {
        VLOG_FUNC_START("BorisSweeper") << " n=" << q.size() << ", t=" << t;
        VLOG(6) << LOG_INDENT << "computing exact solution at t=" << t;
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

          this->exact_cache = make_shared<encap_type>(q);
          this->exact_updated = true;
        } else {
          LOG(DEBUG) << "exact solution has been computed previously.";
          q = *(this->exact_cache.get());
        }
        VLOG(5) << LOG_INDENT << "exact solution at t=" << t << ": " << q;
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::echo_error(const time t, bool predict)
      {
        VLOG_FUNC_START("BorisSweeper");
        auto end = this->end_particles;
        ErrorTuple<scalar> e_tuple;
        scalar e_end = this->impl_solver->energy(end, t);
        e_tuple.e_drift = abs(this->initial_energy - e_end);
        e_tuple.res = this->compute_residual();

        size_t n = this->get_controller()->get_step();
        size_t k = this->get_controller()->get_iteration();
        error_index nk(n, k);

        LOG(INFO) << boost::format(this->FORMAT_STR) % (n+1) % k % e_tuple.res % e_tuple.e_drift % e_end;
        VLOG(3) << LOG_INDENT << "particle at t_end: " << end;

        // exact anlytical solution only valid for 1-particle-system
        if (this->particles[0]->size() == 1) {
          shared_ptr<encap_type> ex = dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution));
          this->exact(ex, t);

          e_tuple.p_err.x = ex->positions()[0][0] - end->positions()[0][0];
          e_tuple.p_err.y = ex->positions()[0][1] - end->positions()[0][1];
          e_tuple.p_err.z = ex->positions()[0][2] - end->positions()[0][2];
          e_tuple.p_err.u = ex->velocities()[0][0] - end->velocities()[0][0];
          e_tuple.p_err.v = ex->velocities()[0][1] - end->velocities()[0][1];
          e_tuple.p_err.w = ex->velocities()[0][2] - end->velocities()[0][2];

          VLOG(2) << LOG_INDENT << "absolute error at end point: " << e_tuple.p_err;
        }
        this->errors.insert(pair<error_index, ErrorTuple<scalar>>(nk, e_tuple));

        this->write_particle_cloud_to_file(n, k, end, e_end, e_tuple.e_drift, e_tuple.res);
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      error_map<scalar> BorisSweeper<scalar, time>::get_errors() const
      {
        return this->errors;
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::setup(bool coarse)
      {
        EncapSweeper<time>::setup(coarse);
        VLOG_FUNC_START("BorisSweeper") << " coarse=" << boolalpha << coarse;
        auto const nodes = this->get_nodes();
        assert(nodes.size() >= 1);
        const size_t nnodes = nodes.size();
        const size_t num_s_integrals = this->quadrature->left_is_node() ? nnodes : nnodes - 1;
        VLOG(4) << "there will be " << num_s_integrals << " integrals for " << nnodes << " nodes";

        // compute delta nodes
        this->delta_nodes = vector<time>(nnodes, time(0.0));
        for (size_t m = 1; m < nnodes; m++) {
          this->delta_nodes[m] = nodes[m] - nodes[m - 1];
        }

        this->start_particles = dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution));
        this->end_particles = dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution));

        this->energy_evals.resize(nnodes);
        for (size_t m = 0; m < nnodes; ++m) {
          this->particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
          this->saved_particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
          this->forces.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->saved_forces.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->b_vecs.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->saved_b_vecs.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          if (coarse) {
            this->tau_corrections.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
          }
        }

        for (size_t m = 0; m < num_s_integrals; ++m) {
          this->s_integrals.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->ss_integrals.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
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

        size_t nsteps = this->get_controller()->get_end_time() / this->get_controller()->get_time_step();
        size_t digit_step = to_string(nsteps + 1).length();
        size_t digit_iter = to_string(this->get_controller()->get_max_iterations() - 1).length();
        this->FORMAT_STR = "step: %|" + to_string(digit_step) + "|      iter: %|" + to_string(digit_iter) + "|"
                           + "      residual: %10.4e      energy drift: %10.4e      total energy: %10.2f";

        UNUSED(coarse);
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::advance()
      {
        VLOG_FUNC_START("BorisSweeper");
        VLOG(6) << LOG_INDENT << "advancing to next step";
        this->start_particles->copy(this->end_particles);
        this->energy_evals.front() = this->energy_evals.back();
        this->forces.front() = this->forces.back();
        this->b_vecs.front() = this->b_vecs.back();
        this->exact_updated = false;
        VLOG(7) << LOG_INDENT << "particles: " << this->particles;
        VLOG(7) << LOG_INDENT << "energies:  " << this->energy_evals;
        VLOG(8) << LOG_INDENT << "forces:    " << this->forces;
        VLOG(8) << LOG_INDENT << "b_vecs:    " << this->b_vecs;
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::evaluate(size_t m)
      {
        VLOG_FUNC_START("BorisSweeper") << " node=" << m;
        time t = this->get_controller()->get_time() + this->get_controller()->get_time_step() * this->delta_nodes[m];
        VLOG(6) << LOG_INDENT << "computing forces at node " << m << " (t=" << t << ")";

        this->forces[m] = this->impl_solver->e_field_evaluate(this->particles[m], t);
        this->b_vecs[m] = this->impl_solver->b_field_vecs(this->particles[m], t);

        VLOG(8) << LOG_INDENT << "particles:" << this->particles[m];
        VLOG(7) << LOG_INDENT << "e_forces:" << this->forces[m];
        VLOG(7) << LOG_INDENT << "b_vecs:" << this->b_vecs[m];
        this->f_evals++;
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::predict(bool initial)
      {
        VLOG_FUNC_START("BorisSweeper") << " initial=" << boolalpha << initial;

        // in any case copy as we have a simple spread as predict
        this->particles[0]->copy(this->start_particles);

        this->spread();
        for (size_t m = 0; m < this->particles.size(); ++m) {
          this->evaluate(m);
        }

        // in any case copy as we have a simple spread as predict
        this->end_particles->copy(this->particles.back());

        this->save();
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::sweep()
      {
        VLOG_FUNC_START("BorisSweeper");
        const auto   nodes  = this->get_nodes();
        const size_t nnodes = nodes.size();
        assert(nnodes >= 1);
        time t  = this->get_controller()->get_time();
        time dt = this->get_controller()->get_time_step();
        VLOG(2) << LOG_INDENT << "sweeping for t=" << t << " and dt=" << dt;
        VLOG(6) << LOG_INDENT << "nodes: " << nodes;

        this->impl_solver->energy(this->particles.front(), t);

        // compute integrals
        VLOG(7) << LOG_INDENT << "computing integrals";
        zero(this->s_integrals);
        zero(this->ss_integrals);
        if (this->get_quadrature()->left_is_node()) {
          // starting at m=1 as m=0 will only add zeros
          for (size_t m = 1; m < nnodes; m++) {
            for (size_t l = 0; l < nnodes; l++) {
              auto rhs = this->build_rhs(l);
              this->s_integrals[m] += rhs * dt * this->s_mat(m, l);
              this->ss_integrals[m] += rhs * dt * dt * this->ss_mat(m, l);
            }
          }
        } else {
          throw NotImplementedYet("left-is-NOT-node");
        }
        VLOG(7) << LOG_INDENT << "s_int:  " << this->s_integrals;
        VLOG(7) << LOG_INDENT << "ss_int: " << this->ss_integrals;

        this->evaluate(0);

        for (size_t m = 0; m < nnodes - 1; m++) {
          time ds = dt * this->delta_nodes[m+1];
          VLOG(3) << LOG_INDENT << "node " << m << ", ds=" << ds;
          pfasst::log::stack_position++;

          //// Update Position (explicit)
          //
          // x_{m+1}^{k+1} = x_{m}^{k+1}
          this->update_position(m, dt, ds);
          VLOG(3) << LOG_INDENT << "new positions: " << this->particles[m+1]->positions();

          // evaluate electric field with new position
          this->forces[m+1] = this->impl_solver->e_field_evaluate(this->particles[m+1], t + nodes[m]);

          //// Update Velocity (semi-implicit)
          this->update_velocity(m, ds, nodes);
          VLOG(3) << LOG_INDENT << "new velocities: " << this->particles[m+1]->velocities();

          pfasst::log::stack_position--;
        }

        if (this->get_quadrature()->right_is_node()) {
          this->end_particles->copy(this->particles.back());
        } else {
          throw NotImplementedYet("right-is-NOT-node");
        }

        this->save();
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::save(bool initial_only)
      {
        VLOG_FUNC_START("BorisSweeper") << " initial_only=" << boolalpha << initial_only;
        if (initial_only) {
          this->saved_particles[0] = make_shared<encap_type>(*(this->particles[0].get()));
          this->saved_forces[0] = this->forces[0];
          this->saved_b_vecs[0] = this->b_vecs[0];
        } else {
          for (size_t m = 0; m < this->saved_particles.size(); m++) {
            VLOG(8) << LOG_INDENT << "node" << m;
            VLOG(8) << LOG_INDENT << "  particle:" << this->particles[m];
            this->saved_particles[m] = make_shared<encap_type>(*(this->particles[m].get()));
            VLOG(8) << LOG_INDENT << "  previous_particle:" << this->saved_particles[m];
          }
          VLOG(8) << LOG_INDENT << "forces:" << this->forces;
          this->saved_forces = this->forces;
          VLOG(8) << LOG_INDENT << "saved_forces:" << this->saved_forces;
          VLOG(8) << LOG_INDENT << "b_vecs:" << this->b_vecs;
          this->saved_b_vecs = this->b_vecs;
          VLOG(8) << LOG_INDENT << "saved_b_vecs:" << this->saved_b_vecs;
        }
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::spread()
      {
        VLOG_FUNC_START("BorisSweeper");
        for (size_t m = 1; m < this->particles.size(); ++m) {
          this->set_state(const_pointer_cast<const encap_type>(this->particles[0]), m);
        }
        VLOG_FUNC_END("BorisSweeper");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post_sweep()
      {
        time t  = this->get_controller()->get_time();
        time dt = this->get_controller()->get_time_step();
        this->echo_error(t + dt);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post_predict()
      {
        time t  = this->get_controller()->get_time();
        time dt = this->get_controller()->get_time_step();
        this->echo_error(t + dt);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post_step()
      {}

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post(ICommunicator* comm, int tag)
      {
        UNUSED(comm); UNUSED(tag);
        // TODO: implement BorisSweeper::post
      };

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::send(ICommunicator* comm, int tag, bool blocking)
      {
        UNUSED(comm); UNUSED(tag); UNUSED(blocking);
        // TODO: implement BorisSweeper::send
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::recv(ICommunicator* comm, int tag, bool blocking)
      {
        UNUSED(comm); UNUSED(tag); UNUSED(blocking);
        // TODO: implement BorisSweeper::recv
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::broadcast(ICommunicator* comm)
      {
        UNUSED(comm);
        // TODO: implement BorisSweeper::broadcast
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::boris_solve(const time tm, const time t_next, const time ds, const size_t m,
                                                   const velocity_type& c_k_term)
      {
        VLOG_FUNC_START("BorisSweeper") << " tm=" << tm << ", t_next=" << t_next << ", ds=" << ds << ", node=" << m;
        VLOG(4) << LOG_INDENT << "solving with Boris' method";
        UNUSED(t_next);
        velocity_type c_k_term_half = c_k_term / scalar(2.0);
        VLOG(6) << LOG_INDENT << "c_k_term/2: " << c_k_term_half;
        AttributeValues<scalar> beta = cmp_wise_div(this->particles[m]->charges(), this->particles[m]->masses()) / scalar(2.0) * ds;
        VLOG(6) << LOG_INDENT << "beta: " << beta;
        acceleration_type e_forces_mean = (this->forces[m] + this->forces[m+1]) / scalar(2.0);
        VLOG(5) << LOG_INDENT << "e_mean: " << e_forces_mean << " (<=" << this->forces[m] << " +" << this->forces[m+1] << " / 2)";

        // first Boris' drift
        //   v^{-} = v^{k}
        velocity_type v_minus = this->particles[m]->velocities();
        //           + \beta * E_{mean} + c^{k} / 2
        v_minus += e_forces_mean * beta + c_k_term_half;
        VLOG(5) << LOG_INDENT << "v-: " << v_minus;

        // Boris' kick
        velocity_type boris_t(beta.size());
        auto b_field_vector = this->impl_solver->get_b_field_vector();
        for (size_t p = 0; p < beta.size(); ++p) {
          boris_t[p] = b_field_vector * beta[p];
        }
        velocity_type v_prime = v_minus + cross_prod(v_minus, boris_t);
        VLOG(5) << LOG_INDENT << "v': " << v_prime;

        // final Boris' drift
        vector<scalar> boris_t_sqr(boris_t.size());
        for (size_t p = 0; p < boris_t.size(); ++p) {
          boris_t_sqr[p] = pow(norm0(boris_t[p]), 2);
        }
        velocity_type boris_s = (scalar(2.0) * boris_t) / (scalar(1.0) + boris_t_sqr);
        velocity_type v_plus = v_minus + cross_prod(v_prime, boris_s);
        VLOG(5) << LOG_INDENT << "v+: " << v_plus;

//           assert(abs(v_minus.norm0() - v_plus.norm0()) <= 10e-8);

        this->particles[m+1]->velocities() = v_plus + e_forces_mean * beta + c_k_term_half;
        VLOG_FUNC_END("BorisSweeper");
      }


      template<typename scalar, typename time>
      BorisSweeper<scalar, time>& as_boris_sweeper(shared_ptr<ISweeper<time>> x)
      {
        shared_ptr<BorisSweeper<scalar, time>> y = dynamic_pointer_cast<BorisSweeper<scalar, time>>(x);
        assert(y);
        return *y.get();
      }


      template<typename scalar, typename time>
      const BorisSweeper<scalar, time>& as_boris_sweeper(shared_ptr<const ISweeper<time>> x)
      {
        shared_ptr<const BorisSweeper<scalar, time>> y = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(x);
        assert(y);
        return *y.get();
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst
