#include "pfasst/encap/imex_sweeper.hpp"

#include <algorithm>
#include <cassert>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace encap
  {
    template<typename time>
    void IMEXSweeper<time>::integrate_end_state(time dt)
    {
      vector<shared_ptr<Encapsulation<time>>> dst = { this->end_state };
      dst[0]->copy(this->start_state);
      dst[0]->mat_apply(dst, dt, this->quadrature->get_b_mat(), this->fs_expl, false);
      dst[0]->mat_apply(dst, dt, this->quadrature->get_b_mat(), this->fs_impl, false);
    }

    template<typename time>
    void IMEXSweeper<time>::setup(bool coarse)
    {
      pfasst::encap::EncapSweeper<time>::setup(coarse);

      auto const num_nodes = this->quadrature->get_num_nodes();
      auto const num_s_integrals = this->quadrature->left_is_node() ? num_nodes - 1 : num_nodes;

      for (size_t m = 0; m < num_nodes; m++) {
        this->fs_expl.push_back(this->get_factory()->create(pfasst::encap::function));
        this->fs_impl.push_back(this->get_factory()->create(pfasst::encap::function));
      }
      for (size_t m = 0; m < num_s_integrals; m++) {
        this->s_integrals.push_back(this->get_factory()->create(pfasst::encap::solution));
      }

      if (! this->quadrature->left_is_node()) {
        this->fs_expl_start = this->get_factory()->create(pfasst::encap::function);
      }

      size_t nsteps = this->get_controller()->get_end_time() / this->get_controller()->get_step_size();
      size_t digit_step = (this->get_controller()->get_step_size() > 0) ? 
                            to_string(nsteps + 1).length() : 3;
      size_t digit_iter = (this->get_controller()->get_max_iterations() > 0) ? 
                            to_string(this->get_controller()->get_max_iterations() - 1).length() : 3;
      this->FORMAT_STR = "step: %|" + to_string(digit_step) + "|      iter: %|" + to_string(digit_iter) + "|"
                         + "      n1: %|2|      n2: %|3|"
                         + "      residual: %10.4e" + "      err: %10.4e";
    }

    template<typename time>
    void IMEXSweeper<time>::predict(bool initial)
    {
      if (this->quadrature->left_is_node()) {
        this->predict_with_left(initial);
      } else {
        this->predict_without_left(initial);
      }

      if (this->quadrature->right_is_node()) {
        this->end_state->copy(this->state.back());
      } else {
        this->integrate_end_state(this->get_controller()->get_step_size());
      }
    }

    template<typename time>
    void IMEXSweeper<time>::sweep()
    {
      if (this->quadrature->left_is_node()) {
        this->sweep_with_left();
      } else {
        this->sweep_without_left();
      }

      if (this->quadrature->right_is_node()) {
        this->end_state->copy(this->state.back());
      } else {
        this->integrate_end_state(this->get_controller()->get_step_size());
      }
    }

    template<typename time>
    void IMEXSweeper<time>::advance()
    {
      this->start_state->copy(this->end_state);
      if (this->quadrature->left_is_node() && this->quadrature->right_is_node()) {
        this->state[0]->copy(this->start_state);
        this->fs_expl.front()->copy(this->fs_expl.back());
        this->fs_impl.front()->copy(this->fs_impl.back());
      }
    }

    template<typename time>
    void IMEXSweeper<time>::reevaluate(bool initial_only)
    {
      time t0 = this->get_controller()->get_time();
      time dt = this->get_controller()->get_step_size();
      if (initial_only) {
        if (this->quadrature->left_is_node()) {
          this->f_expl_eval(this->fs_expl[0], this->state[0], t0);
          this->f_impl_eval(this->fs_impl[0], this->state[0], t0);
        } else {
          throw NotImplementedYet("reevaluate");
        }
      } else {
        for (size_t m = 0; m < this->quadrature->get_num_nodes(); m++) {
          time t =  t0 + dt * this->quadrature->get_nodes()[m];
          this->f_expl_eval(this->fs_expl[m], this->state[m], t);
          this->f_impl_eval(this->fs_impl[m], this->state[m], t);
        }
      }
    }

    template<typename time>
    void IMEXSweeper<time>::integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
    {
      auto const q_mat = this->quadrature->get_q_mat();
      dst[0]->mat_apply(dst, dt, q_mat, this->fs_expl, true);
      dst[0]->mat_apply(dst, dt, q_mat, this->fs_impl, false);
    }

    template<typename time>
    void IMEXSweeper<time>::residual(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
    {
      for (size_t m = 0; m < this->quadrature->get_num_nodes(); m++) {
        dst[m]->copy(this->start_state);
        dst[m]->saxpy(-1.0, this->state[m]);
      }
      if (this->fas_corrections.size() > 0) {
        // XXX: build a matrix and use mat_apply to do this
        for (size_t m = 0; m < this->quadrature->get_num_nodes(); m++) {
          for (size_t n = 0; n <= m; n++) {
            dst[m]->saxpy(1.0, this->fas_corrections[n]);
          }
        }
      }
      dst[0]->mat_apply(dst, dt, this->quadrature->get_q_mat(), this->fs_expl, false);
      dst[0]->mat_apply(dst, dt, this->quadrature->get_q_mat(), this->fs_impl, false);
    }

    template<typename time>
    void IMEXSweeper<time>::f_expl_eval(shared_ptr<Encapsulation<time>> f_expl_encap,
                                        shared_ptr<Encapsulation<time>> u_encap,
                                        time t)
    {
      UNUSED(f_expl_encap); UNUSED(u_encap); UNUSED(t);
      throw NotImplementedYet("imex (f_expl_eval)");
    }

    template<typename time>
    void IMEXSweeper<time>::f_impl_eval(shared_ptr<Encapsulation<time>> f_impl_encap,
                                        shared_ptr<Encapsulation<time>> u_encap,
                                        time t)
    {
      UNUSED(f_impl_encap); UNUSED(u_encap); UNUSED(t);
      throw NotImplementedYet("imex (f_impl_eval)");
    }

    template<typename time>
    void IMEXSweeper<time>::impl_solve(shared_ptr<Encapsulation<time>> f_encap,
                                       shared_ptr<Encapsulation<time>> u_encap,
                                       time t, time dt,
                                       shared_ptr<Encapsulation<time>> rhs_encap)
    {
      UNUSED(f_encap); UNUSED(u_encap); UNUSED(t); UNUSED(dt); UNUSED(rhs_encap);
      throw NotImplementedYet("imex (impl_solve)");
    }

    template<typename time>
    void IMEXSweeper<time>::predict_with_left(bool initial)
    {
      time dt = this->get_controller()->get_step_size();
      time t  = this->get_controller()->get_time();
      ML_CVLOG(1, "Sweeper", "predicting step " << this->get_controller()->get_step() + 1
                             << " (t=" << t << ", dt=" << dt << ")");

      if (initial) {
        this->state[0]->copy(this->start_state);
        this->f_expl_eval(this->fs_expl[0], this->state[0], t);
        this->f_impl_eval(this->fs_impl[0], this->state[0], t);
      }

      shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

      auto const nodes = this->quadrature->get_nodes();

      // step across all nodes
      for (size_t m = 0; m < nodes.size() - 1; ++m) {
        time ds = dt * (nodes[m+1] - nodes[m]);
        rhs->copy(this->state[m]);
        rhs->saxpy(ds, this->fs_expl[m]);
        this->impl_solve(this->fs_impl[m + 1], this->state[m + 1], t, ds, rhs);
        this->f_expl_eval(this->fs_expl[m + 1], this->state[m + 1], t + ds);
        t += ds;
      }
    }

    template<typename time>
    void IMEXSweeper<time>::predict_without_left(bool initial)
    {
      UNUSED(initial);
      time dt = this->get_controller()->get_step_size();
      time t  = this->get_controller()->get_time();
      ML_CVLOG(1, "Sweeper", "predicting step " << this->get_controller()->get_step() + 1
                             << " (t=" << t << ", dt=" << dt << ")");
      time ds;

      shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

      auto const nodes = this->quadrature->get_nodes();

      // step to first node
      ds = dt * nodes[0];
      this->f_expl_eval(this->fs_expl_start, this->start_state, t);
      rhs->copy(this->start_state);
      rhs->saxpy(ds, this->fs_expl_start);
      this->impl_solve(this->fs_impl[0], this->state[0], t, ds, rhs);
      this->f_expl_eval(this->fs_expl[0], this->state[0], t + ds);

      // step across all nodes
      for (size_t m = 0; m < nodes.size() - 1; ++m) {
        ds = dt * (nodes[m+1] - nodes[m]);
        rhs->copy(this->state[m]);
        rhs->saxpy(ds, this->fs_expl[m]);
        this->impl_solve(this->fs_impl[m+1], this->state[m+1], t, ds, rhs);
        this->f_expl_eval(this->fs_expl[m+1], this->state[m+1], t + ds);
        t += ds;
      }
    }

    template<typename time>
    void IMEXSweeper<time>::sweep_with_left()
    {
      auto const nodes = this->quadrature->get_nodes();
      auto const dt    = this->get_controller()->get_step_size();
      auto const s_mat = this->quadrature->get_s_mat().block(1, 0, nodes.size()-1, nodes.size());
      ML_CVLOG(1, "Sweeper", "sweeping on step " << this->get_controller()->get_step() + 1
                             << " in iteration " << this->get_controller()->get_iteration()
                             << " (dt=" << dt << ")");
      time ds;

      this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_expl, true);
      this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_impl, false);
      for (size_t m = 0; m < nodes.size() - 1; ++m) {
        ds = dt * (nodes[m+1] - nodes[m]);
        this->s_integrals[m]->saxpy(-ds, this->fs_expl[m]);
        this->s_integrals[m]->saxpy(-ds, this->fs_impl[m+1]);
      }
      if (this->fas_corrections.size() > 0) {
        for (size_t m = 0; m < this->s_integrals.size(); m++) {
          this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m+1]);
        }
      }

      shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

      // step across all nodes
      auto t = this->get_controller()->get_time();
      for (size_t m = 0; m < nodes.size() - 1; ++m) {
        ds = dt * (nodes[m+1] - nodes[m]);
        rhs->copy(this->state[m]);
        rhs->saxpy(ds, this->fs_expl[m]);
        rhs->saxpy(1.0, this->s_integrals[m]);
        this->impl_solve(this->fs_impl[m + 1], this->state[m + 1], t, ds, rhs);
        this->f_expl_eval(this->fs_expl[m + 1], this->state[m + 1], t + ds);
        t += ds;
      }
    }

    template<typename time>
    void IMEXSweeper<time>::sweep_without_left()
    {
      auto const nodes = this->quadrature->get_nodes();
      auto const dt    = this->get_controller()->get_step_size();
      auto const s_mat = this->quadrature->get_s_mat();
      ML_CVLOG(1, "Sweeper", "sweeping on step " << this->get_controller()->get_step() + 1
                             << " in iteration " << this->get_controller()->get_iteration()
                             << " (dt=" << dt << ")");
      time ds;

      this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_expl, true);
      this->s_integrals[0]->mat_apply(this->s_integrals, dt, s_mat, this->fs_impl, false);
      ds = dt * nodes[0];
      this->s_integrals[0]->saxpy(-ds, this->fs_expl_start);
      this->s_integrals[0]->saxpy(-ds, this->fs_impl[0]);
      for (size_t m = 0; m < nodes.size() - 1; ++m) {
        ds = dt * (nodes[m+1] - nodes[m]);
        this->s_integrals[m+1]->saxpy(-ds, this->fs_expl[m]);
        this->s_integrals[m+1]->saxpy(-ds, this->fs_impl[m+1]);
      }
      if (this->fas_corrections.size() > 0) {
        for (size_t m = 0; m < this->s_integrals.size(); m++) {
          this->s_integrals[m]->saxpy(1.0, this->fas_corrections[m]);
        }
      }

      shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

      // step to first node
      auto t = this->get_controller()->get_time();
      ds = dt * nodes[0];
      this->f_expl_eval(this->fs_expl_start, this->start_state, t);
      rhs->copy(this->start_state);
      rhs->saxpy(ds, this->fs_expl_start);
      rhs->saxpy(1.0, this->s_integrals[0]);
      this->impl_solve(this->fs_impl[0], this->state[0], t, ds, rhs);
      this->f_expl_eval(this->fs_expl[0], this->state[0], t + ds);

      // step across all nodes
      for (size_t m = 0; m < nodes.size() - 1; ++m) {
        ds = dt * (nodes[m+1] - nodes[m]);
        rhs->copy(this->state[m]);
        rhs->saxpy(ds, this->fs_expl[m]);
        rhs->saxpy(1.0, this->s_integrals[m+1]);
        this->impl_solve(this->fs_impl[m+1], this->state[m+1], t, ds, rhs);
        this->f_expl_eval(this->fs_expl[m+1], this->state[m+1], t + ds);
        t += ds;
      }
    }
  }  // ::pfasst::encap
}  // ::pfasst
