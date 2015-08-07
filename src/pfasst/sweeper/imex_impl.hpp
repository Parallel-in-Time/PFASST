#include "pfasst/sweeper/imex.hpp"


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

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, "SWEEPER")
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

    CVLOG(2, "SWEEPER") << "initial values for prediction";
    for (size_t m = 0; m <= num_nodes; ++m) {
      CVLOG(2, "SWEEPER") << LOG_FIXED << "  t["<<m<<"]=" << nodes[m];
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "       u: " << to_string(this->get_states()[m]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "    f_ex: " << to_string(this->_expl_rhs[m]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "    f_im: " << to_string(this->_impl_rhs[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::predict()
  {
    Sweeper<SweeperTrait, Enabled>::predict();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, "SWEEPER")
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const time_type t = this->get_status()->get_time();
    const time_type dt = this->get_status()->get_dt();

    assert(this->get_quadrature() != nullptr);
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), time_type(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    this->_expl_rhs.front() = this->evaluate_rhs_expl(t, this->get_states().front());

    CLOG(INFO, "SWEEPER") << LOG_FIXED << "Predicting from t=" << t << " over " << num_nodes << " nodes"
                          << " to t=" << (t + dt) << " (dt=" << dt << ")";

    time_type tm = t;
    for (size_t m = 0; m < num_nodes; ++m) {
      const time_type ds = dt * (nodes[m + 1] - nodes[m]);

      CVLOG(1, "SWEEPER") << LOG_FIXED << "propagating from t["<<m<<"]=" << (dt * nodes[m])
                          << " to t["<<(m+1)<<"]=" << (dt * nodes[m+1])
                          << " with ds = dt * (t["<<(m+1)<<"] - t["<<m<<"]) = "
                          << ds << " = " << dt << " * (" << nodes[m+1] << " - " << nodes[m] << ")";
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "  u["<<m<<"] = " << to_string(this->get_states()[m]);

      // compute right hand side for implicit solve (i.e. the explicit part of the propagation)
      shared_ptr<encap_type> rhs = this->get_encap_factory()->create();
      rhs->data() = this->get_states()[m]->get_data();
      CVLOG(2, "SWEEPER") << "  rhs = u["<<m<<"]         = " << to_string(rhs);
      rhs->scaled_add(ds, this->_expl_rhs[m]);
      CVLOG(2, "SWEEPER") << "     += ds * f_ex["<<m<<"] = " << LOG_FIXED << ds << " * "
                                                             << LOG_FLOAT << to_string(this->_expl_rhs[m]);
      CVLOG(2, "SWEEPER") << "                     = " << to_string(rhs);

      // solve the implicit part
      CVLOG(2, "SWEEPER") << "  solve(u["<<(m+1)<<"] - ds * f_im["<<(m+1)<<"] = rhs)";
      this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, ds, rhs);
      CVLOG(2, "SWEEPER") << "  u["<<(m+1)<<"] = " << to_string(this->get_states()[m + 1]);

      // reevaluate the explicit part with the new solution value
      tm += ds;
      this->_expl_rhs[m + 1] = this->evaluate_rhs_expl(tm, this->get_states()[m + 1]);

      CVLOG(1, "SWEEPER") << LOG_FIXED << "  ==> values at t["<<(m+1)<<"]=" << (dt * nodes[m+1]);
      CVLOG(1, "SWEEPER") << LOG_FLOAT << "         u["<<m+1<<"]: " << to_string(this->get_states()[m + 1]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "      f_ex["<<m+1<<"]: " << to_string(this->_expl_rhs[m + 1]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "      f_im["<<m+1<<"]: " << to_string(this->_impl_rhs[m + 1]);
      CVLOG(1, "SWEEPER") << "";
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

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, "SWEEPER")
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const time_type dt = this->get_status()->get_dt();
    const auto q_mat = this->get_quadrature()->get_q_mat();
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), time_type(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    CVLOG(2, "SWEEPER") << "initial values for sweeping";
    for (size_t m = 0; m <= num_nodes; ++m) {
      CVLOG(2, "SWEEPER") << "  t["<<m<<"]=" << LOG_FIXED << (dt * nodes[m]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "       u: " << to_string(this->get_states()[m]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "    f_ex: " << to_string(this->_expl_rhs[m]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "    f_im: " << to_string(this->_impl_rhs[m]);
    }

    CVLOG(1, "SWEEPER") << "computing integrals";
    CVLOG(2, "SWEEPER") << "  Q_int     = dt * Q * f_ex";
    this->_q_integrals = encap::mat_mul_vec(dt, q_mat, this->_expl_rhs);
    CVLOG(2, "SWEEPER") << "           += dt * Q * f_im";
    encap::mat_apply(this->_q_integrals, dt, q_mat, this->_impl_rhs, false);

    // notes[0] == 0 thus q[0] -= ds * f_ex[0] can be omitted

    CVLOG(2, "SWEEPER") << "  subtracting function evaluations of previous iteration";
    for (size_t m = 0; m < num_nodes; ++m) {
      // XXX: nodes[m] - nodes[m] ?!
      const time_type ds = dt * (nodes[m + 1] - nodes[m]);

      CVLOG(1, "SWEEPER") << LOG_FIXED << "  Q_int["<<(m+1)<<"] -= ds * f_ex["<<m<<"] = "
                                       << ds << " * " << LOG_FLOAT << to_string(this->_expl_rhs[m]);
      this->_q_integrals[m + 1]->scaled_add(-ds, this->_expl_rhs[m]);
      CVLOG(1, "SWEEPER") << LOG_FIXED << "  Q_int["<<(m+1)<<"] -= ds * f_im["<<(m+1)<<"] = "
                                       << ds << " * " << LOG_FLOAT << to_string(this->_impl_rhs[m+1]);
      this->_q_integrals[m + 1]->scaled_add(-ds, this->_impl_rhs[m + 1]);
    }

    for (size_t m = 0; m < num_nodes + 1; ++m) {
      CVLOG(1, "SWEEPER") << LOG_FLOAT << "  Q_int["<<m<<"] += tau["<<m<<"] = " << to_string(this->get_tau()[m]);
      this->_q_integrals[m]->scaled_add(1.0, this->get_tau()[m]);
      CVLOG(1, "SWEEPER") << LOG_FLOAT << "            = " << to_string(this->_q_integrals[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::sweep()
  {
    Sweeper<SweeperTrait, Enabled>::sweep();

    assert(this->get_quadrature() != nullptr);
    assert(this->get_status() != nullptr);

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, "SWEEPER")
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const time_type t = this->get_status()->get_time();
    const time_type dt = this->get_status()->get_dt();
    auto nodes = this->get_quadrature()->get_nodes();
    nodes.insert(nodes.begin(), time_type(0.0));
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    this->_expl_rhs.front() = this->evaluate_rhs_expl(t, this->get_states().front());

    CLOG(INFO, "SWEEPER") << LOG_FIXED << "Sweeping from t=" << t << " over " << num_nodes
                          << " nodes to t=" << (t + dt) << " (dt=" << dt << ")";

    time_type tm = t;
    // note: m=0 is initial value and not a quadrature node
    for (size_t m = 0; m < num_nodes; ++m) {
      // XXX: nodes[m] - nodes[0] ?!
      const time_type ds = dt * (nodes[m + 1] - nodes.front());

      CVLOG(1, "SWEEPER") << LOG_FIXED << "propagating from t["<<m<<"]=" << (dt * nodes[m])
                          << " to t["<<(m+1)<<"]=" << (dt * nodes[m+1])
                          << " with ds = dt * (t["<<(m+1)<<"] - t["<<m<<"]) = "
                          << ds << " = " << dt << " * (" << nodes[m+1] << " - " << nodes[m] << ")";
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "  u["<<m<<"] = " << to_string(this->get_states()[m]);

      // compute right hand side for implicit solve (i.e. the explicit part of the propagation)
      shared_ptr<encap_type> rhs = this->get_encap_factory()->create();
      rhs->data() = this->get_states()[m]->get_data();
      CVLOG(2, "SWEEPER") << "  rhs = u["<<m<<"]           = " << to_string(rhs);
      rhs->scaled_add(ds, this->_expl_rhs[m]);
      CVLOG(2, "SWEEPER") << "     += ds * f_ex["<<m<<"]   = " << LOG_FIXED << ds << " * "
                                                               << LOG_FLOAT << to_string(this->_expl_rhs[m]);
      rhs->scaled_add(1.0, this->_q_integrals[m + 1]);
      CVLOG(2, "SWEEPER") << "     += 1.0 * q_int["<<(m+1)<<"] = " << to_string(this->_q_integrals[m + 1]);
      CVLOG(2, "SWEEPER") << "      = " << to_string(rhs);

      // solve the implicit part
      CVLOG(2, "SWEEPER") << "  solve(u["<<(m+1)<<"] - ds * f_im["<<(m+1)<<"] = rhs)";
      this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, ds, rhs);
      CVLOG(2, "SWEEPER") << "  u["<<(m+1)<<"] = " << to_string(this->get_states()[m + 1]);

      // reevaluate the explicit part with the new solution value
      tm += ds;
      this->_expl_rhs[m + 1] = this->evaluate_rhs_expl(tm, this->get_states()[m + 1]);

      CVLOG(1, "SWEEPER") << LOG_FIXED << "  ==> values at t["<<(m+1)<<"]=" << (dt * nodes[m+1]);
      CVLOG(1, "SWEEPER") << LOG_FLOAT << "         u["<<m+1<<"]: " << to_string(this->get_states()[m + 1]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "      f_ex["<<m+1<<"]: " << to_string(this->_expl_rhs[m + 1]);
      CVLOG(2, "SWEEPER") << LOG_FLOAT << "      f_im["<<m+1<<"]: " << to_string(this->_impl_rhs[m + 1]);
      CVLOG(1, "SWEEPER") << "";
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
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::advance()
  {
    assert(this->get_end_state() != nullptr);
    CLOG(DEBUG, "SWEEPER") << "advancing";

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
      if (this->get_quadrature()->left_is_node()) {
        assert(this->_expl_rhs.front() != nullptr && this->_impl_rhs.front() != nullptr);

        this->_expl_rhs.front() = this->evaluate_rhs_expl(t0, this->_expl_rhs.front());
        this->_impl_rhs.front() = this->evaluate_rhs_impl(t0, this->_impl_rhs.front());

      } else {
        throw NotImplementedYet("reevaluation for quadrature not containing left interval boundary");
      }

    } else {
      const time_type dt = this->get_status()->get_dt();
      const auto nodes = this->get_quadrature()->get_nodes();

      for (size_t m = 0; this->get_quadrature()->get_num_nodes(); ++m) {
        const time_type t = t0 + dt * nodes[m];
        assert(this->_expl_rhs[m] != nullptr && this->_impl_rhs[m] != nullptr);

        this->_expl_rhs[m] = this->evaluate_rhs_expl(t, this->_expl_rhs[m]);
        this->_impl_rhs[m] = this->evaluate_rhs_impl(t, this->_impl_rhs[m]);
      }
    }
  }


  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::integrate_end_state(const typename SweeperTrait::time_type& dt)
  {
    try {
      Sweeper<SweeperTrait, Enabled>::integrate_end_state(dt);
    } catch (NotImplementedYet niy) {
      assert(this->get_quadrature() != nullptr);
      assert(this->get_initial_state() != nullptr);

      this->end_state()->data() = this->get_initial_state()->get_data();
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_expl_rhs)[0]);
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_impl_rhs)[0]);
      CVLOG(1, "SWEEPER") << "end state: " << to_string(this->get_end_state());
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::compute_residuals()
  {
    CLOG(DEBUG, "SWEEPER") << "computing residuals";

    assert(this->get_status() != nullptr);
    assert(this->get_quadrature() != nullptr);
    assert(this->get_initial_state() != nullptr);

    const time_type dt = this->get_status()->get_dt();
    const size_t num_nodes = this->get_quadrature()->get_num_nodes() + 1;

    for (size_t m = 0; m < num_nodes; ++m) {
      assert(this->get_states()[m] != nullptr);
      assert(this->residuals()[m] != nullptr);

      CVLOG(5, "SWEEPER") << "  res["<<m<<"] = u[0]   = " << to_string(this->get_initial_state());
      this->residuals()[m]->data() = this->get_initial_state()->get_data();
      CVLOG(5, "SWEEPER") << "        -= u["<<m<<"]   = " << to_string(this->get_states()[m]);
      this->residuals()[m]->scaled_add(-1.0, this->get_states()[m]);

      for (size_t n = 0; n <= m; ++n) {
        assert(this->get_tau()[n] != nullptr);
        CVLOG(5, "SWEEPER") << "        += tau["<<n<<"] = " << to_string(this->get_tau()[n]);
        this->residuals()[m]->scaled_add(1.0, this->get_tau()[n]);
      }
    }

    CVLOG(5, "SWEEPER") << "  res += dt * Q * F_ex";
    encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_expl_rhs, false);

    CVLOG(5, "SWEEPER") << "  res += dt * Q * F_im";
    encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_impl_rhs, false);

    CVLOG(2, "SWEEPER") << "  ==>";
    for (size_t m = 0; m < num_nodes; ++m) {
      CVLOG(2, "SWEEPER") << "    |res["<<m<<"]| = " << LOG_FLOAT << this->get_residuals()[m]->norm0()
                          << "    res["<<m<<"] = " << to_string(this->get_residuals()[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_type>
  IMEX<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                 const shared_ptr<typename SweeperTrait::encap_type> u)
  {
    UNUSED(t); UNUSED(u);
    throw NotImplementedYet("evaluation of explicit part of right-hand-side");
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_type>
  IMEX<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                 const shared_ptr<typename SweeperTrait::encap_type> u)
  {
    UNUSED(t); UNUSED(u);
    throw NotImplementedYet("evaluation of implicit part of right-hand-side");
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
    throw NotImplementedYet("spacial solver");
  }
}  // ::pfasst
