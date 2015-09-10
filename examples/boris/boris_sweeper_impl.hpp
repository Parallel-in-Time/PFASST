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


#define BCVLOG(level) CVLOG(level, "Boris") << this->log_indent->indent(level)


/*
 * VLOG Levels for 'Boris' logger:
 *  1: predict, sweep, setup
 *  2: same as 1 but more verose (will include basic transfer notes)
 *  3: boris_solve
 *  4: update_position, update_velocity
 *  5: same as 3+4 but more verbose (will include verbose transfer notes)
 *  6: integrate, evaluate
 *  7: build_rhs
 *  8: compute_residual, exact solution
 *  9: printing, verbose save/advance/spread
 */

namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      template<typename precision>
      void ParticleError<precision>::log(el::base::type::ostream_t& os) const
      {
        os << "pos: [" << x << " " << y << " " << z << "]\tvel: [" << y << " " << v << " " << w << "]";
      }


      template<typename precision>
      static void init_opts()
      {
        pfasst::config::options::add_option<size_t>("Boris-SDC", "num_particles", "number of particles in the cloud");
        pfasst::config::options::add_option<precision>("Boris-SDC", "epsilon", "Boris' epsilon");
        pfasst::config::options::add_option<precision>("Boris-SDC", "omega_e", "E-field constant");
        pfasst::config::options::add_option<precision>("Boris-SDC", "omega_b", "B-field constant");
      }

      template<typename precision>
      static void init_logs()
      {
        pfasst::log::add_custom_logger("Boris");
        pfasst::log::add_custom_logger("SolverBinding");
        pfasst::log::add_custom_logger("Solver");
      }


      LogIndent::LogIndent()
      {
        this->vlog_levels.fill(0);
      }
      LogIndent::~LogIndent()
      {}

      void LogIndent::increment(const size_t vlevel)
      {
        this->vlog_levels.at(vlevel - 1)++;
      }
      void LogIndent::decrement(const size_t vlevel)
      {
        this->vlog_levels.at(vlevel - 1)--;
      }
      const string LogIndent::indent(const size_t vlevel) const
      {
        size_t count = 0;
        for(size_t vl = 0; vl < vlevel; ++vl) {
          count += this->vlog_levels.at(vl);
        }
        return string(count * 2, ' ');
      }


      template<typename scalar, typename time>
      typename BorisSweeper<scalar, time>::acceleration_type
      BorisSweeper<scalar, time>::build_rhs(const size_t m, const bool previous) const
      {
        BCVLOG(7) << "building rhs for node " << m << " of "
                           << ((previous) ? "previous" : "current") << " sweep";
        this->log_indent->increment(7);
        acceleration_type rhs = (previous) ? this->saved_forces[m] : this->forces[m];
        BCVLOG(7) << "e-forces: " << rhs;

        if (previous) {
          rhs += cross_prod(this->saved_particles[m]->velocities(), this->saved_b_vecs[m]);
        } else {
          rhs += cross_prod(this->particles[m]->velocities(), this->b_vecs[m]);
        }

        BCVLOG(7) << "=> rhs: " << rhs;
        this->log_indent->decrement(7);
        return rhs;
      }

      template<typename scalar, typename time>
      scalar BorisSweeper<scalar, time>::compute_residual_max()
      {
        BCVLOG(8) << "computing max residual";
        this->log_indent->increment(8);

        this->residual(this->get_controller()->get_step_size(), this->residuals);

        scalar max_residual = scalar(0.0);

        for (size_t m = 1; m < this->residuals.size(); ++m) {
          auto residual_m = dynamic_pointer_cast<encap_type>(this->residuals[m]);
          assert(residual_m);
          max_residual = std::max(max_residual, residual_m->norm0());
        }

        BCVLOG(8) << "=> max residual: " << max_residual;
        return max_residual;
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::write_center_to_file(const size_t iter, const size_t sweep,
                                                            const ParticleComponent<scalar>& center,
                                                            const scalar energy, const scalar drift, const scalar residual)
      {
        BCVLOG(9) << "writing center particle to file";
        this->data_stream_fmt % (iter+1) % sweep % -1
                              % center[0] % center[1] % center[2]
                              // cppcheck-suppress zerodiv
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
        this->log_indent->increment(9);
        for (size_t p = 0; p < cloud->size(); ++p) {
          BCVLOG(9) << "writing cloud particle " << p << " to file";
          // cppcheck-suppress zerodiv
          this->data_stream_fmt % (iter+1) % sweep % p
                                % cloud->positions()[p * cloud->dim()] % cloud->positions()[p * cloud->dim() + 1] % cloud->positions()[p * cloud->dim() + 2]
                                % cloud->velocities()[p * cloud->dim()] % cloud->velocities()[p * cloud->dim() + 1] % cloud->velocities()[p * cloud->dim() + 2]
                                % energy % drift % residual;
          this->data_stream << this->data_stream_fmt << endl;
        }
        if (with_center) {
          this->write_center_to_file(iter, sweep, cloud->center_of_mass(), energy, drift, residual);
        }
        this->log_indent->decrement(9);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::update_position(const size_t m, const time dt, const time ds)
      {
        BCVLOG(4) << "updating position (" << m << "->" << m+1 << ") with dt=" << dt << ", ds=" << ds;
        this->log_indent->increment(4);
        this->particles[m+1]->positions() = this->particles[m]->positions();
        BCVLOG(4) << "old: " << this->particles[m+1]->positions();
        this->log_indent->increment(5);
        //               + delta_nodes_{m+1} * v_{0}
        this->particles[m+1]->positions() += this->start_particles->velocities() * ds;
        BCVLOG(5) << "+= " << this->start_particles->velocities() << " * " << ds;
        //               + \sum_{l=1}^{m} sx_{m+1,l}^{x} (f_{l}^{k+1} - f_{l}^{k})
        for (size_t l = 0; l <= m; l++) {
          auto rhs_this = this->build_rhs(l);
          auto rhs_prev = this->build_rhs(l, true);
          this->particles[m+1]->positions() += rhs_this * dt * dt * this->sx_mat(m+1, l);
          BCVLOG(5) << "+= " << rhs_this << " * " << dt * dt << " * " << this->sx_mat(m+1, l);
          this->particles[m+1]->positions() -= rhs_prev * dt * dt * this->sx_mat(m+1, l);
          BCVLOG(5) << "-= " << rhs_prev << " * " << dt * dt << " * " << this->sx_mat(m+1, l);
        }

        this->particles[m+1]->positions() += this->ss_integrals[m+1];
        BCVLOG(5) << "+= " << this->ss_integrals[m+1];
        this->log_indent->decrement(5);
        this->log_indent->decrement(4);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::update_velocity(const size_t m, const time ds, const vector<time>& nodes)
      {
        BCVLOG(4) << "updating velocity (" << m << "->" << m+1 << ") with ds=" << ds;
        this->log_indent->increment(4);
        velocity_type c_k_term = cloud_component_factory<scalar>(this->particles[0]->size(), this->particles[0]->dim());
        zero(c_k_term);  // TODO: check if this is required

        BCVLOG(5) << "c_k: " << c_k_term;
        this->log_indent->increment(5);
        //                 - delta_nodes_{m} / 2 * f_{m+1}^{k}
        auto t1 = this->build_rhs(m+1, true);
        c_k_term -= t1 * 0.5 * ds;
        BCVLOG(5) << "-= 0.5 * " << t1 << " * " << ds << "  => " << c_k_term;
        //                 - delta_nodes_{m} / 2 * f_{m}^{k}
        auto t2 = this->build_rhs(m, true);
        c_k_term -= t2 * 0.5 * ds;
        BCVLOG(5) << "-= 0.5 * " << t2 << " * " << ds << "  => " << c_k_term;
        //                 + s_integral[m]
        c_k_term += this->s_integrals[m+1];
        BCVLOG(5) << "+= " << this->s_integrals[m+1];
        this->log_indent->decrement(5);
        BCVLOG(4) << "=> c_k: " << c_k_term;

        // doing Boris' magic
        this->boris_solve(nodes[m], nodes[m+1], ds, m, c_k_term);
        this->log_indent->decrement(4);
      }


      template<typename scalar, typename time>
      BorisSweeper<scalar, time>::BorisSweeper(shared_ptr<bindings::WrapperInterface<scalar, time>>& impl_solver,
                                               const string& data_file)
        :   impl_solver(impl_solver)
          , errors()
          , exact_updated(false)
          , log_indent(make_shared<LogIndent>())
          , f_evals(0)
          , data_stream(data_file, ios_base::out | ios_base::trunc)
      {
        assert(data_stream.is_open() && data_stream.good());
        CLOG(INFO, "Boris") << "writing particle data to: " << data_file;
        // CSV format specification:  [step],[iter],[particle],[x],[y],[z],[u],[v],[w],[energy],[drift],[residual]
        this->DATA_STREAM_FORMAT_STR = "%d,%d,%d,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f,%.16f";
        //                              i  sw p  x    y    z    u    v    w    engy edr  res
        BCVLOG(2) << "formatting string: '" << this->DATA_STREAM_FORMAT_STR << "'";
        this->data_stream_fmt = boost::format(this->DATA_STREAM_FORMAT_STR);
      }


      template<typename scalar, typename time>
      BorisSweeper<scalar, time>::~BorisSweeper()
      {
        CLOG(INFO, "Boris") << "number force computations: " << this->f_evals;
      }


      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::set_state(shared_ptr<const encap_type> u0, size_t m)
      {
        this->particles[m]->copy(u0);
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
      shared_ptr<Encapsulation<time>> BorisSweeper<scalar, time>::get_start_state() const
      {
        return this->start_particles;
      }

      template<typename scalar, typename time>
      shared_ptr<typename BorisSweeper<scalar, time>::acceleration_type>
      BorisSweeper<scalar, time>::get_tau_q_as_force(size_t m) const
      {
        return this->tau_q_corrections[m];
      }

      template<typename scalar, typename time>
      shared_ptr<typename BorisSweeper<scalar, time>::acceleration_type>
      BorisSweeper<scalar, time>::get_tau_qq_as_force(size_t m) const
      {
        return this->tau_qq_corrections[m];
      }

      template<typename scalar, typename time>
      shared_ptr<Encapsulation<time>> BorisSweeper<scalar, time>::get_saved_state(size_t m) const
      {
        return this->saved_particles[m];
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::set_initial_energy()
      {
        BCVLOG(1) << "computing and setting initial energy";
        this->log_indent->increment(1);

        auto p0 = this->start_particles;
        BCVLOG(2) << "initial particles: " << p0;
        this->initial_energy = this->impl_solver->energy(p0, this->get_controller()->get_time());
        CLOG(INFO, "Boris") << "initial total energy of system: " << this->initial_energy;

        this->write_particle_cloud_to_file(0, 0, p0, this->initial_energy, 0, 0);
        this->log_indent->decrement(1);
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
        BCVLOG(8) << "computing exact solution at t=" << t;
        this->log_indent->increment(8);
        UNUSED(t);
        if (!this->exact_updated) {
          typedef complex<scalar> C;
          C i(0.0, 1.0);
          auto initial = this->particles[0];
          scalar x0 = initial->positions()[0],
                 y0 = initial->positions()[1],
                 z0 = initial->positions()[2],
                 u0 = initial->velocities()[0],
                 v0 = initial->velocities()[1],
                 w0 = initial->velocities()[2],
                 omega_e = this->impl_solver->omega_e(),
                 omega_b = this->impl_solver->omega_b(),
                 epsilon = this->impl_solver->epsilon();
          time dt = this->get_controller()->get_step_size();

          C omega_tilde = sqrt(-2.0 * epsilon) * omega_e;
          q.positions()[2] = (z0 * cos(omega_tilde * (scalar)(dt)) 
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
          q.positions()[0] = x_y_move.real();
          q.positions()[1] = x_y_move.imag();

          q.velocities()[2] = (- z0 * omega_tilde * sin(omega_tilde * (scalar)(dt)) + w0 * cos(omega_tilde * (scalar)(dt))).real();
          C u_v_move = (- i * omega_plus * (r_plus + i * i_plus)) * exp(-i * omega_plus * (scalar)(dt))
                                     - (i * omega_minus * (r_minus + i * i_minus)) * exp(-i * omega_minus * (scalar)(dt));
          q.velocities()[0] = u_v_move.real();
          q.velocities()[1] = u_v_move.imag();

          this->exact_cache = make_shared<encap_type>(q);
          this->exact_updated = true;
        } else {
          BCVLOG(8) << "exact solution has been computed previously.";
          q = *(this->exact_cache.get());
        }
        BCVLOG(8) << "exact solution at t=" << t << ": " << q;
        this->log_indent->decrement(8);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::echo_error(const time t, bool predict)
      {
        auto end = this->end_particles;
        ErrorTuple<scalar> e_tuple;
        scalar e_end = this->impl_solver->energy(end, t);
        e_tuple.e_drift = abs(this->initial_energy - e_end);
        e_tuple.res = this->compute_residual_max();

        size_t n = this->get_controller()->get_step();
        size_t k = this->get_controller()->get_iteration();
        error_index nk(n, k);

        CLOG(INFO, "Boris") << boost::format(this->FORMAT_STR) % (n+1) % k % (this->coarse ? "coarse" : "fine") % e_tuple.res % e_tuple.e_drift % e_end;
        BCVLOG(9) << "particle at t_end: " << end;

        // exact anlytical solution only valid for 1-particle-system
        if (this->particles[0]->size() == 1) {
          shared_ptr<encap_type> ex = dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution));
          this->exact(ex, t);

          e_tuple.p_err.x = ex->positions()[0] - end->positions()[0];
          e_tuple.p_err.y = ex->positions()[1] - end->positions()[1];
          e_tuple.p_err.z = ex->positions()[2] - end->positions()[2];
          e_tuple.p_err.u = ex->velocities()[0] - end->velocities()[0];
          e_tuple.p_err.v = ex->velocities()[1] - end->velocities()[1];
          e_tuple.p_err.w = ex->velocities()[2] - end->velocities()[2];

          BCVLOG(9) << "absolute error at end point: " << e_tuple.p_err;
        }
        this->errors.insert(pair<error_index, ErrorTuple<scalar>>(nk, e_tuple));

        this->write_particle_cloud_to_file(n, k, end, e_end, e_tuple.e_drift, e_tuple.res);
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
        BCVLOG(1) << "setting up Boris Sweeper for " << ((coarse) ? "coarse" : "fine") << " level";
        this->log_indent->increment(1);
        this->coarse = coarse;
        auto const nodes = this->get_nodes();
        assert(nodes.size() >= 1);
        const size_t nnodes = nodes.size();
        const size_t num_s_integrals = this->quadrature->left_is_node() ? nnodes : nnodes - 1;
        BCVLOG(2) << "there will be " << num_s_integrals << " integrals for " << nnodes << " nodes";

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
          this->residuals.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
          this->saved_particles.push_back(dynamic_pointer_cast<encap_type>(this->get_factory()->create(pfasst::encap::solution)));
          this->forces.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->saved_forces.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->b_vecs.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->saved_b_vecs.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          if (coarse) {
            this->tau_q_corrections.push_back(make_shared<acceleration_type>(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim())));
            this->tau_qq_corrections.push_back(make_shared<acceleration_type>(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim())));
          }
        }

        for (size_t m = 0; m < num_s_integrals; ++m) {
          this->s_integrals.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
          this->ss_integrals.push_back(cloud_component_factory<scalar>(this->particles[m]->size(), this->particles[m]->dim()));
        }

        this->q_mat = this->get_quadrature()->get_q_mat();
        this->qq_mat = this->q_mat * this->q_mat;
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
        this->qt_mat = 0.5 * (qe_mat + qi_mat);

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
          this->ss_mat.row(i) = this->qq_mat.row(i) - this->qq_mat.row(i - 1);
          this->sx_mat.row(i) = this->qx_mat.row(i) - this->qx_mat.row(i - 1);
          this->st_mat.row(i) = this->qt_mat.row(i) - this->qt_mat.row(i - 1);
        }

        size_t nsteps = this->get_controller()->get_end_time() / this->get_controller()->get_step_size();
        size_t digit_step = to_string(nsteps + 1).length();
        size_t digit_iter = to_string(this->get_controller()->get_max_iterations() - 1).length();
        this->FORMAT_STR = "step: %|" + to_string(digit_step) + "|      iter: %|" + to_string(digit_iter) + "| (%-6s)"
                           + "      residual: %10.4e      energy drift: %10.4e      total energy: %10.2f";

        UNUSED(coarse);
        this->log_indent->decrement(1);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
      {
        throw NotImplementedYet("Boris::integrate for basic Encap type");
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::integrate(time dt, vector<shared_ptr<acceleration_type>> dst_q,
                                                 vector<shared_ptr<acceleration_type>> dst_qq) const
      {
        BCVLOG(6) << "integrating over dt=" << dt;
        this->log_indent->increment(6);
        const size_t nnodes = this->get_nodes().size();

        vector<acceleration_type> rhs(nnodes);
        for (size_t m = 0; m < nnodes; ++m) {
          rhs[m] = this->build_rhs(m);
        }

        for (size_t m = 0; m < nnodes; ++m) {
          zero(dst_q[m]);
          zero(dst_qq[m]);
          for (size_t n = 0; n < nnodes; ++n) {
            // for positions
            *(dst_qq[m].get()) += rhs[n] * dt * dt * this->qq_mat(m, n)
                                  + this->start_particles->velocities() * dt * this->q_mat(m, n);
            // for velocities
            *(dst_q[m].get()) += rhs[n] * dt * this->q_mat(m, n);
          }
          BCVLOG(6) << "integral(QQ)[" << m << "]: " << *(dst_qq[m].get());
          BCVLOG(6) << "integral(Q)[" << m << "]:  " << *(dst_q[m].get());
        }
        this->log_indent->decrement(6);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
      {
        BCVLOG(8) << "computing residual";
        const auto   nodes  = this->get_nodes();
        const size_t nnodes = nodes.size();
        assert(dst.size() == nnodes);

        vector<shared_ptr<encap_type>> dst_cast(nnodes);
        vector<shared_ptr<acceleration_type>> q_int(nnodes), qq_int(nnodes);
        for (size_t m = 0; m < nnodes; ++m) {
          dst_cast[m] = dynamic_pointer_cast<encap_type>(dst[m]);
          assert(dst_cast[m]);

          qq_int[m] = make_shared<acceleration_type>(this->start_particles->size() * this->start_particles->dim());
          q_int[m] = make_shared<acceleration_type>(this->start_particles->size() * this->start_particles->dim());
        }

        // cf. pySDC: QF(u)
        this->integrate(dt, q_int, qq_int);

        for (size_t m = 1; m < nnodes; ++m) {
          BCVLOG(8) << "for node " << m;
          zero(dst_cast[m]->positions());
          zero(dst_cast[m]->velocities());

          dst_cast[m]->positions() = *(qq_int[m].get());
          dst_cast[m]->velocities() = *(q_int[m].get());

          // cf. pySDC: L.u[0] - L.u[m+1] (in boris_2nd_order#compute_residual())
          dst_cast[m]->positions() += this->start_particles->positions() - this->particles[m]->positions();
          dst_cast[m]->velocities() += this->start_particles->velocities() - this->particles[m]->velocities();

          // add tau correction (if available)
          if (this->tau_q_corrections.size() > 0 && this->tau_qq_corrections.size()) {
            dst_cast[m]->positions() += *(this->tau_qq_corrections[m].get());
            dst_cast[m]->velocities() += *(this->tau_q_corrections[m].get());
          }
        }

        this->log_indent->decrement(8);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::advance()
      {
        BCVLOG(2) << "advancing to next step";
        this->log_indent->increment(2);
        this->start_particles->copy(this->end_particles);
        this->energy_evals.front() = this->energy_evals.back();
        this->forces.front() = this->forces.back();
        this->b_vecs.front() = this->b_vecs.back();
        this->exact_updated = false;
        BCVLOG(9) << "new starting values:";
        BCVLOG(9) << "  => start_particles: " << this->start_particles;
        BCVLOG(9) << "  => energies:        " << this->energy_evals.front();
        BCVLOG(9) << "  => forces:          " << this->forces.front();
        BCVLOG(9) << "  => b_vecs:          " << this->b_vecs.front();
        this->log_indent->decrement(2);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::evaluate(size_t m)
      {
        time t = this->get_controller()->get_time() + this->get_controller()->get_step_size() * this->delta_nodes[m];
        BCVLOG(2) << "computing forces at node " << m << " (t=" << t << ")";
        this->log_indent->increment(2);

#ifndef BORIS_SAME_LEVELS
        if (this->coarse) {
          BCVLOG(2) << "only external electric field (because on coarse level)";
          // don't compute internal electric forces when on coarse level
          this->forces[m] = this->impl_solver->external_e_field_evaluate(this->particles[m], t);
        } else {
          BCVLOG(2) << "internal and external electric field (because not on coarse level)";
          this->forces[m] = this->impl_solver->e_field_evaluate(this->particles[m], t);
        }
#else
        BCVLOG(2) << "internal and external electric field";
        this->forces[m] = this->impl_solver->e_field_evaluate(this->particles[m], t);
#endif
        this->b_vecs[m] = this->impl_solver->b_field_vecs(this->particles[m], t);

        BCVLOG(9) << "for particles:" << this->particles[m];
        BCVLOG(9) << "  => e_forces:" << this->forces[m];
        BCVLOG(9) << "  => b_vecs:  " << this->b_vecs[m];
        this->log_indent->decrement(2);
        this->f_evals++;
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::predict(bool initial)
      {
        BCVLOG(1) << "predicting with initial particle cloud: " << *(this->start_particles.get());
        this->log_indent->increment(1);
        // in any case copy as we have a simple spread as predict
        this->particles[0]->copy(this->start_particles);

        this->spread();
        for (size_t m = 0; m < this->particles.size(); ++m) {
          this->evaluate(m);
        }

        // in any case copy as we have a simple spread as predict
        this->end_particles->copy(this->particles.back());

        this->save();
        this->log_indent->decrement(1);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::sweep()
      {
        const auto   nodes  = this->get_nodes();
        const size_t nnodes = nodes.size();
        assert(nnodes >= 1);
        time t  = this->get_controller()->get_time();
        time dt = this->get_controller()->get_step_size();
        BCVLOG(1) << "sweeping for t=" << t << " and dt=" << dt;
        this->log_indent->increment(1);
        BCVLOG(2) << "with nodes: " << nodes;
        BCVLOG(2) << "initial: " << *(this->start_particles.get());
        BCVLOG(2) << "previous particles:";
        for (size_t m = 0; m < nnodes; ++m) {
          BCVLOG(2) << "  [" << m << "]: " << *(this->saved_particles[m].get());
        }
        BCVLOG(2) << "current particles:";
        for (size_t m = 0; m < nnodes; ++m) {
          BCVLOG(2) << "  [" << m << "]: " << *(this->particles[m].get());
        }

        this->impl_solver->energy(this->particles.front(), t);

        // compute integrals
        BCVLOG(1) << "computing integrals";
        if (this->get_quadrature()->left_is_node()) {
          // starting at m=1 as m=0 will only add zeros
          for (size_t m = 1; m < nnodes; m++) {
            zero(this->s_integrals[m]);
            zero(this->ss_integrals[m]);
            for (size_t l = 0; l < nnodes; l++) {
              auto rhs = this->build_rhs(l);
              this->s_integrals[m] += rhs * dt * this->s_mat(m, l);
              this->ss_integrals[m] += rhs * dt * dt * this->ss_mat(m, l);
            }
          }
          if (this->tau_q_corrections.size() > 0 && this->tau_qq_corrections.size()) {
            BCVLOG(2) << "adding FAS correction to integrals";
            this->log_indent->increment(2);
            for (size_t m = 0; m < nnodes; ++m) {
              BCVLOG(2) << "+= tau_q[" << m << "]  (<" << this->tau_q_corrections[m]  << ">" << *(this->tau_q_corrections[m].get()) << ")";
              BCVLOG(2) << "+= tau_qq[" << m << "] (<" << this->tau_qq_corrections[m] << ">" << *(this->tau_qq_corrections[m].get()) << ")";
              this->s_integrals[m] += *(this->tau_q_corrections[m].get());
              this->ss_integrals[m] += *(this->tau_qq_corrections[m].get());
              if (m > 0) {
                BCVLOG(2) << "-= tau_q[" << m-1 << "]  (<" << this->tau_q_corrections[m-1]  << ">" << *(this->tau_q_corrections[m-1].get()) << ")";
                BCVLOG(2) << "-= tau_qq[" << m-1 << "] (<" << this->tau_qq_corrections[m-1] << ">" << *(this->tau_qq_corrections[m-1].get()) << ")";
                this->s_integrals[m] -= *(this->tau_q_corrections[m-1].get());
                this->ss_integrals[m] -= *(this->tau_qq_corrections[m-1].get());
              }
            }
            this->log_indent->decrement(2);
          }
        } else {
          throw NotImplementedYet("left-is-NOT-node");
        }
        BCVLOG(2) << "s_int:  " << this->s_integrals;
        BCVLOG(2) << "ss_int: " << this->ss_integrals;

        this->evaluate(0);

        for (size_t m = 0; m < nnodes - 1; m++) {
          time ds = dt * this->delta_nodes[m+1];
          BCVLOG(1) << "node " << m << " (ds=" << ds << ")";
          this->log_indent->increment(1);
          BCVLOG(2) << "old m+1 particle: " << this->particles[m+1];

          //// Update Position (explicit)
          //
          // x_{m+1}^{k+1} = x_{m}^{k+1}
          this->update_position(m, dt, ds);
          BCVLOG(1) << "new positions: " << this->particles[m+1]->positions();

          // evaluate electric field with new position
#ifndef BORIS_SAME_LEVELS
          if (this->coarse) {
            BCVLOG(2) << "only external electric field (because on coarse level)";
            // don't compute internal electric forces when on coarse level
            this->forces[m+1] = this->impl_solver->external_e_field_evaluate(this->particles[m+1], t + nodes[m]);
          } else {
            BCVLOG(2) << "internal and external electric field (because not on coarse level)";
            this->forces[m+1] = this->impl_solver->e_field_evaluate(this->particles[m+1], t + nodes[m]);
          }
#else
          BCVLOG(2) << "internal and external electric field";
          this->forces[m+1] = this->impl_solver->e_field_evaluate(this->particles[m+1], t + nodes[m]);
#endif

          //// Update Velocity (semi-implicit)
          this->update_velocity(m, ds, nodes);
          BCVLOG(1) << "new velocities: " << this->particles[m+1]->velocities();

          this->log_indent->decrement(1);
        }

        if (this->get_quadrature()->right_is_node()) {
          this->end_particles->copy(this->particles.back());
        } else {
          throw NotImplementedYet("right-is-NOT-node");
        }

        this->save();
        this->log_indent->decrement(1);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::save(bool initial_only)
      {
        BCVLOG(2) << "saving current state" << ((initial_only) ? " (only initial)" : "");
        this->log_indent->increment(2);
        if (initial_only) {
          this->saved_particles[0] = make_shared<encap_type>(*(this->particles[0].get()));
          this->saved_forces[0] = this->forces[0];
          this->saved_b_vecs[0] = this->b_vecs[0];
        } else {
          for (size_t m = 0; m < this->saved_particles.size(); m++) {
            BCVLOG(9) << "node " << m;
            BCVLOG(9) << "  particle:          " << this->particles[m];
            this->saved_particles[m] = make_shared<encap_type>(*(this->particles[m].get()));
            BCVLOG(9) << "  previous_particle: " << this->saved_particles[m];
          }
          BCVLOG(9) << "forces:       " << this->forces;
          this->saved_forces = this->forces;
          BCVLOG(9) << "saved_forces: " << this->saved_forces;
          BCVLOG(9) << "b_vecs:       " << this->b_vecs;
          this->saved_b_vecs = this->b_vecs;
          BCVLOG(9) << "saved_b_vecs: " << this->saved_b_vecs;
        }
        this->log_indent->decrement(2);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::spread()
      {
        for (size_t m = 1; m < this->particles.size(); ++m) {
          this->set_state(const_pointer_cast<const encap_type>(this->particles[0]), m);
        }
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post_sweep()
      {
        time t  = this->get_controller()->get_time();
        time dt = this->get_controller()->get_step_size();
        this->echo_error(t + dt);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post_predict()
      {
        time t  = this->get_controller()->get_time();
        time dt = this->get_controller()->get_step_size();
        this->echo_error(t + dt);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post_step()
      {}

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::post(ICommunicator* comm, int tag)
      {
        this->start_particles->post(comm, tag);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::send(ICommunicator* comm, int tag, bool blocking)
      {
        this->end_particles->send(comm, tag, blocking);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::recv(ICommunicator* comm, int tag, bool blocking)
      {
        this->start_particles->recv(comm, tag, blocking);
        if (this->quadrature->left_is_node()) {
          this->particles[0]->copy(this->start_particles);
        }
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::broadcast(ICommunicator* comm)
      {
        if (comm->rank() == comm->size() - 1) {
          this->start_particles->copy(this->end_particles);
        }
        this->start_particles->broadcast(comm);
      }

      template<typename scalar, typename time>
      void BorisSweeper<scalar, time>::boris_solve(const time tm, const time t_next, const time ds, const size_t m,
                                                   const velocity_type& c_k_term)
      {
        BCVLOG(3) << "solving with Boris' method";
        this->log_indent->increment(3);
        const size_t npart = this->start_particles->size();
        UNUSED(t_next);

        velocity_type c_k_term_half = c_k_term / scalar(2.0);
        BCVLOG(5) << "c_k_term/2: " << c_k_term_half;
        vector<scalar> beta = cmp_wise_div(this->particles[m]->charges(), this->particles[m]->masses()) / scalar(2.0) * ds;
        BCVLOG(5) << "beta: " << beta;
        acceleration_type e_forces_mean = (this->forces[m] + this->forces[m+1]) / scalar(2.0);
        BCVLOG(5) << "e_mean: " << e_forces_mean << " (<=" << this->forces[m] << " +" << this->forces[m+1] << " / 2)";

        // first Boris' drift
        //   v^{-} = v^{k}
        velocity_type v_minus = this->particles[m]->velocities();
        //           + \beta * E_{mean} + c^{k} / 2
        v_minus += e_forces_mean * beta + c_k_term_half;
        BCVLOG(3) << "v-: " << v_minus;

        // Boris' kick
        vector<scalar> b_field_vector = this->impl_solver->get_b_field_vector();
        velocity_type boris_t = kronecker(beta, b_field_vector);
        velocity_type v_prime = v_minus + cross_prod(v_minus, boris_t);
        BCVLOG(3) << "v': " << v_prime;

        // final Boris' drift
        vector<scalar> boris_t_sqr = norm_sq_npart(boris_t, npart);  // particle-wise scalar product of boris_t
//         for (size_t p = 0; p < boris_t.size(); ++p) {
//           boris_t_sqr[p] = pow(norm0(boris_t[p]), 2);
//         }
        velocity_type boris_s = (boris_t * scalar(2.0)) / (boris_t_sqr + scalar(1.0));
        velocity_type v_plus = v_minus + cross_prod(v_prime, boris_s);
        BCVLOG(3) << "v+: " << v_plus;

//           assert(abs(v_minus.norm0() - v_plus.norm0()) <= 10e-8);

        this->particles[m+1]->velocities() = v_plus + e_forces_mean * beta + c_k_term_half;
        this->log_indent->decrement(3);
      }
    }  // ::pfasst::examples::boris
  }  // ::pfasst::examples
}  // ::pfasst
