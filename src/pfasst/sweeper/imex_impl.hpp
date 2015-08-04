#include "pfasst/sweeper/imex.hpp"

namespace pfasst
{
  template<class SweeperTrait, typename Enabled>
  IMEX<SweeperTrait, Enabled>::IMEX()
    :   Sweeper<SweeperTrait, Enabled>()
      , _q_integrals(0)
      , _expl_rhs(0)
      , _impl_rhs(0)
  {}

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::setup()
  {
    pfasst::Sweeper<SweeperTrait, Enabled>::setup();

    CLOG_IF(this->get_quadrature()->left_is_node(), WARNING, "SWEEPER")
      << "IMEX Sweeper for quadrature nodes containing t_0 not implemented and tested.";

    const auto num_nodes = this->quadrature()->get_num_nodes();
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
  {}

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
    const auto nodes = this->get_quadrature()->get_nodes();
    const size_t num_nodes = nodes.size();

    time_type tm = t;
    for (size_t m = 0; m < num_nodes; ++m) {
      const time_type ds = dt * (nodes[m + 1] - nodes.front());

      shared_ptr<encap_type> rhs = this->get_states()[m];
      rhs->scaled_add(ds, this->_expl_rhs[m]);
      this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, ds, rhs);
      tm += ds;
      this->_expl_rhs[m + 1] = this->evaluate_rhs_expl(tm, this->get_states()[m + 1]);
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
    const auto nodes = this->get_quadrature()->get_nodes();
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    this->_q_integrals = encap::mat_mul_vec(dt, q_mat, this->_expl_rhs);
    encap::mat_apply(this->_q_integrals, dt, q_mat, this->_impl_rhs, false);

    for (size_t m = 1; m < num_nodes; ++m) {
      // XXX: nodes[m] - nodes[0] ?!
      const size_t ds = dt * (nodes[m] - nodes.front());

      this->_q_integrals[m + 1]->scaled_add(-ds, this->_expl_rhs[m]);
      this->_q_integrals[m + 1]->scaled_add(-ds, this->_impl_rhs[m + 1]);
    }
    for (size_t m = 0; m < this->_q_integrals.size(); ++m) {
      this->_q_integrals[m]->scaled_add(1.0, this->get_tau()[m]);
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
    const auto nodes = this->get_quadrature()->get_nodes();
    const size_t num_nodes = this->get_quadrature()->get_num_nodes();

    this->_expl_rhs.front() = this->evaluate_rhs_expl(t, this->get_states().front());

    size_t tm = t;
    // note: m=0 is initial value and not a quadrature node
    for (size_t m = 0; m < num_nodes; ++m) {
      // XXX: nodes[m] - nodes[0] ?!
      const size_t ds = dt * (nodes[m + 1] - nodes.front());

      // we don't reevaluate F at first note as it should have been done after setting it
      // TODO: may need to add some checks for that later...

      shared_ptr<encap_type> rhs = this->get_initial_state();
      rhs->data() = this->get_states()[m]->data();
      rhs->scaled_add(ds, this->_expl_rhs[m]);
      rhs->scaled_add(1.0, this->_q_integrals[m + 1]);
      this->implicit_solve(this->_impl_rhs[m + 1], this->states()[m + 1], tm, ds, rhs);
      tm += ds;
      this->_expl_rhs[m + 1] = this->evaluate_rhs_expl(tm, this->get_states()[m + 1]);
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
  {}

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::advance()
  {
    assert(this->get_end_state() != nullptr);

    this->initial_state()->data() = this->get_end_state()->get_data();
    assert(this->get_quadrature() != nullptr);

    if (this->get_quadrature()->left_is_node() && this->get_quadrature()->right_is_node()) {
      assert(this->_expl_rhs.front() != nullptr && this->_expl_rhs.back() != nullptr);
      assert(this->_impl_rhs.front() != nullptr && this->_impl_rhs.back() != nullptr);

      this->_expl_rhs.front()->data() = this->_expl_rhs.back()->get_data();
      this->_impl_rhs.front()->data() = this->_impl_rhs.back()->get_data();
    } else {
      // TODO: this might not be necessary as it is probably dealt with in pre_predict and pre_sweep
      throw NotImplementedYet("advancing IMEX for nodes not containing left and right time interval borders");
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

      this->end_state() = this->get_initial_state();
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_expl_rhs)[0]);
      this->end_state()->scaled_add(1.0, encap::mat_mul_vec(dt, this->get_quadrature()->get_b_mat(), this->_impl_rhs)[0]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  IMEX<SweeperTrait, Enabled>::compute_residuals()
  {
    assert(this->get_status() != nullptr);
    assert(this->get_quadrature() != nullptr);
    const time_type& dt = this->get_status()->get_dt();

    for (size_t m = 0; this->get_quadrature()->get_num_nodes(); ++m) {
      assert(this->get_initial_state() != nullptr);
      assert(this->get_states()[m] != nullptr);

      this->residuals()[m] = this->get_initial_state();
      this->residuals()[m]->scaled_add(-1.0, this->get_states()[m]);

      for (size_t n = 0; n <= m; ++n) {
        assert(this->get_tau()[n] != nullptr);
        this->residuals()[m]->scaled_add(1.0, this->get_tau()[n]);
      }
    }

    encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_expl_rhs, false);
    encap::mat_apply(this->residuals(), dt, this->get_quadrature()->get_q_mat(), this->_impl_rhs, false);
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
