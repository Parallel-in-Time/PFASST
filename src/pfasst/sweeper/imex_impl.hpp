#include "pfasst/sweeper/imex.hpp"

#include <stdexcept>
using namespace std;


namespace pfasst
{
  template<class SweeperTrait, typename Enabled>
  IMEX<SweeperTrait, Enabled>::IMEX()
    :   Sweeper<SweeperTrait, Enabled>()
      , _q_integrals(0)
      , _expl_rhs(0)
      , _impl_rhs(0)
      , _num_expl_f_evals(0)
      , _num_impl_f_evals(0)
      , _num_impl_solves(0)
  {}

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::setup()
  {
    pfasst::Sweeper<SweeperTrait, Enabled>::setup();

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id())
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const auto num_nodes = this->get_quadrature()->get_num_nodes();
    assert(this->get_states().size() == num_nodes + 1);

    this->_q_integrals.resize(num_nodes + 1);
    generate(this->_q_integrals.begin(), this->_q_integrals.end(),
             bind(&encap_type::factory_type::create, this->encap_factory()));

    this->_expl_rhs.resize(num_nodes + 1);
    generate(this->_expl_rhs.begin(), this->_expl_rhs.end(),
             bind(&encap_type::factory_type::create, this->encap_factory()));

    this->_impl_rhs.resize(num_nodes + 1);
    generate(this->_impl_rhs.begin(), this->_impl_rhs.end(),
             bind(&encap_type::factory_type::create, this->encap_factory()));

    this->compute_delta_matrices();
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::pre_predict()
  {
    Sweeper<SweeperTrait, Enabled>::pre_predict();

    assert(this->get_quadrature() != nullptr);
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), time_type(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    CVLOG(2, this->get_logger_id()) << "initial values for prediction";
    for (size_t m = 0; m <= num_nodes; ++m) {
      CVLOG(2, this->get_logger_id()) << LOG_FIXED << "  t["<<m<<"]=" << nodes[m];
      CVLOG(2, this->get_logger_id()) << LOG_FLOAT << "       u: " << to_string(this->get_states()[m]);
      CVLOG(2, this->get_logger_id()) << LOG_FLOAT << "    f_ex: " << to_string(this->_expl_rhs[m]);
      CVLOG(2, this->get_logger_id()) << LOG_FLOAT << "    f_im: " << to_string(this->_impl_rhs[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::predict()
  {
    Sweeper<SweeperTrait, Enabled>::predict();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id())
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const time_type t = this->get_status()->get_time();
    const time_type dt = this->get_status()->get_dt();

    assert(this->get_quadrature() != nullptr);
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), time_type(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    this->_expl_rhs.front() = this->evaluate_rhs_expl(t, this->get_states().front());

    CLOG(INFO, this->get_logger_id()) << LOG_FIXED << "Predicting from t=" << t << " over " << num_nodes << " nodes"
                          << " to t=" << (t + dt) << " (dt=" << dt << ")";

    time_type tm = t;
    for (size_t m = 0; m < num_nodes; ++m) {
      CVLOG(1, this->get_logger_id()) << LOG_FIXED << "propagating from t["<<m<<"]=" << dt << " * " << nodes[m]
                          << " to t["<<(m+1)<<"]=" << dt << " * " << nodes[m+1];
      CVLOG(2, this->get_logger_id()) << LOG_FLOAT << "  u["<<m<<"] = " << to_string(this->get_states()[m]);

      // compute right hand side for implicit solve (i.e. the explicit part of the propagation)
      shared_ptr<encap_type> rhs = this->get_encap_factory()->create();
      rhs->data() = this->get_states()[m]->get_data();
      CVLOG(2, this->get_logger_id()) << "  rhs = u["<<m<<"]                    = " << to_string(rhs);
      rhs->scaled_add(dt * this->_q_delta_expl(m + 1, m), this->_expl_rhs[m]);
      CVLOG(2, this->get_logger_id()) << "     += dt * QE_{"<<(m+1)<<","<<m<<"} * f_ex["<<m<<"] = "
                          << LOG_FIXED << dt << " * " << this->_q_delta_expl(m + 1, m) << " * "
                          << LOG_FLOAT << to_string(this->_expl_rhs[m]);
      CVLOG(2, this->get_logger_id()) << "                                = " << to_string(rhs);

      // solve the implicit part
      CVLOG(2, this->get_logger_id()) << "  solve(u["<<(m+1)<<"] - dt * QI_{"<<(m+1)<<","<<(m+1)<<"} * f_im["<<(m+1)<<"] = rhs)";
      this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, dt * this->_q_delta_impl(m + 1, m + 1), rhs);
      CVLOG(2, this->get_logger_id()) << "  u["<<(m+1)<<"] = " << to_string(this->get_states()[m + 1]);

      // reevaluate the explicit part with the new solution value
      tm += dt * this->_q_delta_expl(m + 1, m);
      this->_expl_rhs[m + 1] = this->evaluate_rhs_expl(tm, this->get_states()[m + 1]);

      CVLOG(1, this->get_logger_id()) << LOG_FIXED << "  ==> values at t["<<(m+1)<<"]=" << (dt * nodes[m+1]);
      CVLOG(1, this->get_logger_id()) << LOG_FLOAT << "         u["<<m+1<<"]: " << to_string(this->get_states()[m + 1]);
      CVLOG(2, this->get_logger_id()) << LOG_FLOAT << "      f_ex["<<m+1<<"]: " << to_string(this->_expl_rhs[m + 1]);
      CVLOG(2, this->get_logger_id()) << LOG_FLOAT << "      f_im["<<m+1<<"]: " << to_string(this->_impl_rhs[m + 1]);
      CVLOG(1, this->get_logger_id()) << "";
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::post_predict()
  {
    Sweeper<SweeperTrait, Enabled>::post_predict();
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::pre_sweep()
  {
    Sweeper<SweeperTrait, Enabled>::pre_sweep();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id())
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const time_type dt = this->get_status()->get_dt();
    const auto q_mat = this->get_quadrature()->get_q_mat();
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), time_type(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    CVLOG(4, this->get_logger_id()) << "initial values for sweeping";
    for (size_t m = 0; m <= num_nodes; ++m) {
      CVLOG(5, this->get_logger_id()) << "  t["<<m<<"]=" << LOG_FIXED << dt << " * " << nodes[m];
      CVLOG(6, this->get_logger_id()) << LOG_FLOAT << "       u: " << to_string(this->get_states()[m]);
      CVLOG(6, this->get_logger_id()) << LOG_FLOAT << "    f_ex: " << to_string(this->_expl_rhs[m]);
      CVLOG(6, this->get_logger_id()) << LOG_FLOAT << "    f_im: " << to_string(this->_impl_rhs[m]);
    }

    CVLOG(4, this->get_logger_id()) << "computing integrals";
    CVLOG(6, this->get_logger_id()) << "  q_int     = dt * Q * f_ex";
    this->_q_integrals = encap::mat_mul_vec(dt, q_mat, this->_expl_rhs);
    CVLOG(6, this->get_logger_id()) << "           += dt * Q * f_im";
    encap::mat_apply(this->_q_integrals, dt, q_mat, this->_impl_rhs, false);

    CVLOG(4, this->get_logger_id()) << "  subtracting function evaluations of previous iteration and adding FAS correction";

    // XXX do we need to do that ?! isn't it always zero ?!
    CVLOG(6, this->get_logger_id()) << LOG_FLOAT << "  q_int[0] += tau[0]                  = " << to_string(this->get_tau()[0]);

    for (size_t m = 0; m < num_nodes; ++m) {
      for (size_t n = 0; n < m + 1; ++n) {
        CVLOG(6, this->get_logger_id()) << LOG_FIXED << "  q_int["<<m<<"] -= dt * QE_{"<<(m+1)<<","<<n<<"} * f_ex["<<n<<"] = "
                                         << -dt << " * " << this->_q_delta_expl(m + 1, n) << " * "
                                         << LOG_FLOAT << to_string(this->_expl_rhs[n]);
        this->_q_integrals[m + 1]->scaled_add(-dt * this->_q_delta_expl(m + 1, n), this->_expl_rhs[n]);

        CVLOG(6, this->get_logger_id()) << LOG_FIXED << "  q_int["<<(m+1)<<"] -= dt * QI_{"<<(m+1)<<","<<(n+1)<<"} * f_im["<<(n+1)<<"] = "
                                         << -dt << " * " << this->_q_delta_impl(m + 1, n + 1) << " * "
                                         << LOG_FLOAT << to_string(this->_impl_rhs[n+1]);
        this->_q_integrals[m + 1]->scaled_add(-dt * this->_q_delta_impl(m + 1, n + 1), this->_impl_rhs[n + 1]);
      }

      CVLOG(6, this->get_logger_id()) << LOG_FLOAT << "  q_int["<<(m+1)<<"] += tau["<<(m+1)<<"]                  = "
                                       << to_string(this->get_tau()[m + 1]);
      this->_q_integrals[m + 1]->scaled_add(1.0, this->get_tau()[m + 1]);
    }

    for (size_t m = 0; m < num_nodes + 1; ++m) {
      CVLOG(5, this->get_logger_id()) << LOG_FLOAT << "  q_int["<<m<<"] = " << to_string(this->_q_integrals[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::sweep()
  {
    Sweeper<SweeperTrait, Enabled>::sweep();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, this->get_logger_id())
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const time_type t = this->get_status()->get_time();
    const time_type dt = this->get_status()->get_dt();
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), time_type(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    CLOG(INFO, this->get_logger_id()) << LOG_FIXED << "Sweeping from t=" << t << " over " << num_nodes
                          << " nodes to t=" << (t + dt) << " (dt=" << dt << ")";

    this->_expl_rhs.front() = this->evaluate_rhs_expl(t, this->get_states().front());

    time_type tm = t;
    // note: m=0 is initial value and not a quadrature node
    for (size_t m = 0; m < num_nodes; ++m) {
      CVLOG(4, this->get_logger_id()) << LOG_FIXED << "propagating from t["<<m<<"]=" << dt << " * " << nodes[m]
                          << " to t["<<(m+1)<<"]=" << dt << " * " << nodes[m+1];
      CVLOG(5, this->get_logger_id()) << LOG_FLOAT << "  u["<<m<<"] = " << to_string(this->get_states()[m]);

      // compute right hand side for implicit solve (i.e. the explicit part of the propagation)
      shared_ptr<encap_type> rhs = this->get_encap_factory()->create();
      // rhs = u_0
      rhs->data() = this->get_states().front()->get_data();
      CVLOG(6, this->get_logger_id()) << "  rhs = u[0]                    = " << to_string(rhs);

      // rhs += dt * \sum_{i=0}^m (QI_{m+1,i} fI(u_i^{k+1}) + QE_{m+1,i-1} fE(u_{i-1}^{k+1}) ) + QE_{m+1,m} fE(u_{m}^{k+1})
      for (size_t n = 0; n <= m; ++n) {
        rhs->scaled_add(dt * this->_q_delta_impl(m + 1, n), this->_impl_rhs[n]);
        CVLOG(6, this->get_logger_id()) << "     += dt * QI_{"<<(m+1)<<","<<(n+1)<<"} * f_im["<<(n+1)<<"] = "
                            << LOG_FIXED << dt << " * " << this->_q_delta_impl(m + 1, n + 1) << " * "
                            << LOG_FLOAT << to_string(this->_impl_rhs[n]);

        rhs->scaled_add(dt * this->_q_delta_expl(m + 1, n), this->_expl_rhs[n]);
        CVLOG(6, this->get_logger_id()) << "     += dt * QE_{"<<(m+1)<<","<<n<<"} * f_ex["<<n<<"] = "
                            << LOG_FIXED << dt << " * " << this->_q_delta_expl(m + 1, n) << " * "
                            << LOG_FLOAT << to_string(this->_expl_rhs[n]);
      }

      rhs->scaled_add(1.0, this->_q_integrals[m + 1]);
      CVLOG(6, this->get_logger_id()) << "     += 1.0 * q_int["<<(m+1)<<"]          = " << to_string(this->_q_integrals[m + 1]);
      CVLOG(6, this->get_logger_id()) << "      = " << to_string(rhs);

      // solve the implicit part
      CVLOG(4, this->get_logger_id()) << "  solve(u["<<(m+1)<<"] - dt * QI_{"<<(m+1)<<","<<(m+1)<<"} * f_im["<<(m+1)<<"] = rhs)";
      this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, dt * this->_q_delta_impl(m+1, m+1), rhs);
      CVLOG(5, this->get_logger_id()) << "  u["<<(m+1)<<"] = " << to_string(this->get_states()[m + 1]);

      // reevaluate the explicit part with the new solution value
      tm += dt * this->_q_delta_impl(m+1, m+1);
      this->_expl_rhs[m + 1] = this->evaluate_rhs_expl(tm, this->get_states()[m + 1]);

      CVLOG(4, this->get_logger_id()) << LOG_FIXED << "  ==> values at t["<<(m+1)<<"]=" << (dt * nodes[m+1]);
      CVLOG(5, this->get_logger_id()) << LOG_FLOAT << "         u["<<m+1<<"]: " << to_string(this->get_states()[m + 1]);
      CVLOG(6, this->get_logger_id()) << LOG_FLOAT << "      f_ex["<<m+1<<"]: " << to_string(this->_expl_rhs[m + 1]);
      CVLOG(6, this->get_logger_id()) << LOG_FLOAT << "      f_im["<<m+1<<"]: " << to_string(this->_impl_rhs[m + 1]);
      CVLOG(4, this->get_logger_id()) << "";
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::post_sweep()
  {
    Sweeper<SweeperTrait, Enabled>::post_sweep();
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::post_step()
  {
    Sweeper<SweeperTrait, Enabled>::post_step();

    const time_type t_end = this->get_status()->get_time() + this->get_status()->get_dt();
    CLOG(INFO, this->get_logger_id()) << LOG_FIXED << "Solution at t_end=" << t_end << ": "
                          << LOG_FLOAT << to_string(this->get_end_state());
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::advance(const size_t& num_steps)
  {
    UNUSED(num_steps);

    assert(this->get_end_state() != nullptr);
    CVLOG(4, this->get_logger_id()) << "advancing";

    this->initial_state()->data() = this->get_end_state()->get_data();
    assert(this->get_quadrature() != nullptr);

    if (this->get_quadrature()->left_is_node() && this->get_quadrature()->right_is_node()) {
      assert(this->_expl_rhs.front() != nullptr && this->_expl_rhs.back() != nullptr);
      assert(this->_impl_rhs.front() != nullptr && this->_impl_rhs.back() != nullptr);

      this->_expl_rhs.front()->data() = this->_expl_rhs.back()->get_data();
      this->_impl_rhs.front()->data() = this->_impl_rhs.back()->get_data();
    } else {
      // TODO: this might not be necessary as it is probably dealt with in pre_predict and pre_sweep
//       throw NotImplementedYet("advancing IMEX for nodes not containing left and right time interval borders");
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::reevaluate(const bool initial_only)
  {
    assert(this->get_status() != nullptr);
    assert(this->get_quadrature() != nullptr);

    const time_type t0 = this->get_status()->get_time();

    if (initial_only) {
      assert(this->_expl_rhs.front() != nullptr && this->_impl_rhs.front() != nullptr);

      this->_expl_rhs.front() = this->evaluate_rhs_expl(t0, this->get_initial_state());
      this->_impl_rhs.front() = this->evaluate_rhs_impl(t0, this->get_initial_state());

    } else {
      const time_type dt = this->get_status()->get_dt();
      auto nodes = this->get_quadrature()->get_nodes();
      nodes.insert(nodes.begin(), 0.0);

      for (size_t m = 0; m < this->get_quadrature()->get_num_nodes() + 1; ++m) {
        const time_type t = t0 + dt * nodes[m];
        assert(this->_expl_rhs[m] != nullptr && this->_impl_rhs[m] != nullptr);

        this->_expl_rhs[m] = this->evaluate_rhs_expl(t, this->get_states()[m]);
        this->_impl_rhs[m] = this->evaluate_rhs_impl(t, this->get_states()[m]);
      }
    }
  }

  template<class SweeperTrait, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_type>>
  IMEX<SweeperTrait, Enabled>::integrate(const typename SweeperTrait::time_type& dt)
  {
    auto const q_mat = this->get_quadrature()->get_q_mat();

    auto result = encap::mat_mul_vec(dt, q_mat, this->_expl_rhs);
    encap::mat_apply(result, dt, q_mat, this->_impl_rhs, false);

    return result;
  }


  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::integrate_end_state(const typename SweeperTrait::time_type& dt)
  {
    try {
      Sweeper<SweeperTrait, Enabled>::integrate_end_state(dt);
    } catch (runtime_error err) {
      assert(this->get_quadrature() != nullptr);
      assert(this->get_initial_state() != nullptr);

      this->end_state()->data() = this->get_initial_state()->get_data();
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_expl_rhs)[0]);
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_impl_rhs)[0]);
      CVLOG(1, this->get_logger_id()) << "end state: " << to_string(this->get_end_state());
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::compute_residuals()
  {
    CVLOG(4, this->get_logger_id()) << "computing residuals";

    assert(this->get_status() != nullptr);
    assert(this->get_quadrature() != nullptr);
    assert(this->get_initial_state() != nullptr);

    const time_type dt = this->get_status()->get_dt();
    const size_t num_nodes = this->get_quadrature()->get_num_nodes() + 1;

    for (size_t m = 0; m < num_nodes; ++m) {
      assert(this->get_states()[m] != nullptr);
      assert(this->residuals()[m] != nullptr);

      CVLOG(5, this->get_logger_id()) << "  res["<<m<<"] = u[0]   = " << to_string(this->get_initial_state());
      this->residuals()[m]->data() = this->get_initial_state()->get_data();
      CVLOG(5, this->get_logger_id()) << "        -= u["<<m<<"]   = " << to_string(this->get_states()[m]);
      this->residuals()[m]->scaled_add(-1.0, this->get_states()[m]);

      assert(this->get_tau()[m] != nullptr);
      CVLOG(5, this->get_logger_id()) << "        += tau["<<m<<"] = " << to_string(this->get_tau()[m]);
      this->residuals()[m]->scaled_add(1.0, this->get_tau()[m]);
    }

    CVLOG(5, this->get_logger_id()) << "  res += dt * Q * F_ex";
    encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_expl_rhs, false);

    CVLOG(5, this->get_logger_id()) << "  res += dt * Q * F_im";
    encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_impl_rhs, false);

    CVLOG(5, this->get_logger_id()) << "  ==>";
    for (size_t m = 0; m < num_nodes; ++m) {
      CVLOG(5, this->get_logger_id()) << "    |res["<<m<<"]| = " << LOG_FLOAT << this->get_residuals()[m]->norm0()
                                      << "    res["<<m<<"] = " << to_string(this->get_residuals()[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_type>
  IMEX<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                 const shared_ptr<typename SweeperTrait::encap_type> u)
  {
    UNUSED(t); UNUSED(u);
    throw runtime_error("evaluation of explicit part of right-hand-side");
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_type>
  IMEX<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                 const shared_ptr<typename SweeperTrait::encap_type> u)
  {
    UNUSED(t); UNUSED(u);
    throw runtime_error("evaluation of implicit part of right-hand-side");
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                              shared_ptr<typename SweeperTrait::encap_type> u,
                                              const typename SweeperTrait::time_type& t,
                                              const typename SweeperTrait::time_type& dt,
                                              const shared_ptr<typename SweeperTrait::encap_type> rhs)
  {
    UNUSED(f); UNUSED(u); UNUSED(t); UNUSED(dt); UNUSED(dt); UNUSED(rhs);
    throw runtime_error("spacial solver");
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::compute_delta_matrices()
  {
    assert(this->get_quadrature() != nullptr);
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();
    auto nodes = this->get_quadrature()->get_nodes();
    if (this->get_quadrature()->left_is_node()) {
      CLOG(ERROR, this->get_logger_id()) << "Don't know how to compute delta matrices for quadrature containing left time point.";
      throw runtime_error("IMEX with quadrature containing left time point");
    } else {
      nodes.insert(nodes.begin(), time_type(0.0));
    }
    
    CVLOG(4, this->get_logger_id()) << "computing Q_delta matrices for IMEX scheme";

    this->_q_delta_expl = Matrix<time_type>::Zero(num_nodes + 1, num_nodes + 1);
    this->_q_delta_impl = Matrix<time_type>::Zero(num_nodes + 1, num_nodes + 1);

    for (size_t m = 1; m < num_nodes + 1; ++m) {
      for (size_t n = m; n < num_nodes + 1; ++n) {
        this->_q_delta_expl(n, m - 1) = nodes[m] - nodes[m - 1];
        this->_q_delta_impl(n, m) = nodes[m] - nodes[m - 1];
      }
    }

    CVLOG(5, this->get_logger_id()) << "QE:";
    for (int row = 0; row < this->_q_delta_expl.rows(); ++row) {
      CVLOG(5, this->get_logger_id()) << "  " << LOG_FIXED << this->_q_delta_expl.block(row, 0, 1, this->_q_delta_expl.cols());
    }

    CVLOG(5, this->get_logger_id()) << "QI:";
    for (int row = 0; row < this->_q_delta_impl.rows(); ++row) {
      CVLOG(5, this->get_logger_id()) << "  " << LOG_FIXED << this->_q_delta_impl.block(row, 0, 1, this->_q_delta_impl.cols());
    }
  }
}  // ::pfasst
