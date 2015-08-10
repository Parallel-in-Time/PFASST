#include "pfasst/sweeper/sweeper.hpp"

#include <algorithm>
#include <cassert>
using namespace std;

#include "pfasst/exceptions.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"
using pfasst::quadrature::IQuadrature;


namespace pfasst
{
  template<class SweeperTrait, typename Enabled>
  Sweeper<SweeperTrait, Enabled>::Sweeper()
    :   _quadrature(nullptr)
      , _factory(make_shared<typename SweeperTrait::encap_type::factory_type>())
      , _states(0)
      , _previous_states(0)
      , _end_state(nullptr)
      , _tau(0)
      , _residuals(0)
      , _status(nullptr)
      , _abs_residual_tol(0.0)
      , _rel_residual_tol(0.0)
  {}

  template<class SweeperTrait, typename Enabled>
  shared_ptr<IQuadrature<typename SweeperTrait::time_type>>&
  Sweeper<SweeperTrait, Enabled>::quadrature()
  {
    return this->_quadrature;
  }

  template<class SweeperTrait, typename Enabled>
  const shared_ptr<IQuadrature<typename SweeperTrait::time_type>>
  Sweeper<SweeperTrait, Enabled>::get_quadrature() const
  {
    return this->_quadrature;
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<Status<typename SweeperTrait::time_type>>&
  Sweeper<SweeperTrait, Enabled>::status()
  {
    return this->_status;
  }

  template<class SweeperTrait, typename Enabled>
  const shared_ptr<Status<typename SweeperTrait::time_type>>
  Sweeper<SweeperTrait, Enabled>::get_status() const
  {
    return this->_status;
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_type::factory_type>&
  Sweeper<SweeperTrait, Enabled>::encap_factory()
  {
    return this->_factory;
  }

  template<class SweeperTrait, typename Enabled>
  const shared_ptr<typename SweeperTrait::encap_type::factory_type>
  Sweeper<SweeperTrait, Enabled>::get_encap_factory() const
  {
    return this->_factory;
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_type>&
  Sweeper<SweeperTrait, Enabled>::initial_state()
  {
    if (this->get_states().size() == 0) {
      CLOG(ERROR, "SWEEPER") << "Sweeper need to be setup first before querying initial state.";
      throw runtime_error("sweeper not setup before querying initial state");
    }
    return this->states().front();
  }

  template<class SweeperTrait, typename Enabled>
  const shared_ptr<typename SweeperTrait::encap_type>
  Sweeper<SweeperTrait, Enabled>::get_initial_state() const
  {
    if (this->get_states().size() == 0) {
      CLOG(ERROR, "SWEEPER") << "Sweeper need to be setup first before querying initial state.";
      throw runtime_error("sweeper not setup before querying initial state");
    }
    return this->get_states().front();
  }

  template<class SweeperTrait, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_type>>&
  Sweeper<SweeperTrait, Enabled>::states()
  {
    return this->_states;
  }

  template<class SweeperTrait, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_type>>
  Sweeper<SweeperTrait, Enabled>::get_states() const
  {
    return this->_states;
  }

  template<class SweeperTrait, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_type>>&
  Sweeper<SweeperTrait, Enabled>::previous_states()
  {
    return this->_previous_states;
  }

  template<class SweeperTrait, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_type>>
  Sweeper<SweeperTrait, Enabled>::get_previous_states() const
  {
    return this->_previous_states;
  }

  template<class SweeperTrait, typename Enabled>
  shared_ptr<typename SweeperTrait::encap_type>&
  Sweeper<SweeperTrait, Enabled>::end_state()
  {
    return this->_end_state;
  }

  template<class SweeperTrait, typename Enabled>
  const shared_ptr<typename SweeperTrait::encap_type>
  Sweeper<SweeperTrait, Enabled>::get_end_state() const
  {
    return this->_end_state;
  }

  template<class SweeperTrait, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_type>>&
  Sweeper<SweeperTrait, Enabled>::tau()
  {
    return this->_tau;
  }

  template<class SweeperTrait, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_type>>
  Sweeper<SweeperTrait, Enabled>::get_tau() const
  {
    return this->_tau;
  }

  template<class SweeperTrait, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_type>>&
  Sweeper<SweeperTrait, Enabled>::residuals()
  {
    return this->_residuals;
  }

  template<class SweeperTrait, typename Enabled>
  const vector<shared_ptr<typename SweeperTrait::encap_type>>
  Sweeper<SweeperTrait, Enabled>::get_residuals() const
  {
    return this->_residuals;
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::set_options()
  {
    CLOG(DEBUG, "SWEEPER") << "setting options";
    this->_abs_residual_tol = config::get_value<typename traits::spacial_type>("abs_res_tol", this->_abs_residual_tol);
    this->_rel_residual_tol = config::get_value<typename traits::spacial_type>("rel_res_tol", this->_rel_residual_tol);
    CVLOG(1, "SWEEPER") << "absolut residual tolerance:  " << this->_abs_residual_tol;
    CVLOG(1, "SWEEPER") << "relative residual tolerance: " << this->_rel_residual_tol;
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::set_abs_residual_tol(const typename SweeperTrait::spacial_type& abs_res_tol)
  {
    this->_abs_residual_tol = abs_res_tol;
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::set_rel_residual_tol(const typename SweeperTrait::spacial_type& rel_res_tol)
  {
    this->_rel_residual_tol = rel_res_tol;
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::setup()
  {
    if (this->get_status() == nullptr) {
      throw runtime_error("Status not yet set.");
    }
    CLOG(DEBUG, "SWEEPER") << "setting up sweeper with " << to_string(this->get_status());

    if (this->get_quadrature() == nullptr) {
      throw runtime_error("Quadrature not yet set.");
    }
    assert(this->get_encap_factory() != nullptr);

    const auto nodes = this->get_quadrature()->get_nodes();
    const auto num_nodes = this->get_quadrature()->get_num_nodes();

    this->states().resize(num_nodes + 1);
    generate(this->states().begin(), this->states().end(),
             bind(&traits::encap_type::factory_type::create, this->get_encap_factory()));

    this->previous_states().resize(num_nodes + 1);
    generate(this->previous_states().begin(), this->previous_states().end(),
             bind(&traits::encap_type::factory_type::create, this->get_encap_factory()));

    this->end_state() = this->get_encap_factory()->create();

    this->tau().resize(num_nodes + 1);
    generate(this->tau().begin(), this->tau().end(),
             bind(&traits::encap_type::factory_type::create, this->get_encap_factory()));

    this->residuals().resize(num_nodes + 1);
    generate(this->residuals().begin(), this->residuals().end(),
             bind(&traits::encap_type::factory_type::create, this->get_encap_factory()));
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::pre_predict()
  {
    CLOG(DEBUG, "SWEEPER") << "pre-predicting";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::predict()
  {
    CLOG(DEBUG, "SWEEPER") << "predicting";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::post_predict()
  {
    CLOG(DEBUG, "SWEEPER") << "post-predicting";

    assert(this->get_status() != nullptr);
    this->integrate_end_state(this->get_status()->get_dt());

    assert(this->get_quadrature() != nullptr);
    CVLOG(3, "SWEEPER") << "solution at nodes:";
    for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
      CVLOG(3, "SWEEPER") << "  " << m << ": " << to_string(this->get_states()[m]);
    }
    CVLOG(3, "SWEEPER") << "solution at t_end: " << to_string(this->get_end_state());
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::pre_sweep()
  {
    CLOG(DEBUG, "SWEEPER") << "pre-sweeping";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::sweep()
  {
    CLOG(DEBUG, "SWEEPER") << "sweeping";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::post_sweep()
  {
    CLOG(DEBUG, "SWEEPER") << "post-sweeping";

    assert(this->get_status() != nullptr);
    this->integrate_end_state(this->get_status()->get_dt());

    assert(this->get_quadrature() != nullptr);
    CVLOG(3, "SWEEPER") << "solution at nodes:";
    for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
      CVLOG(3, "SWEEPER") << "\t" << m << ": " << to_string(this->get_states()[m]);
    }
    CVLOG(3, "SWEEPER") << "solution at t_end:" << to_string(this->get_end_state());
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::post_step()
  {
    CLOG(DEBUG, "SWEEPER") << "post step";

    assert(this->get_quadrature() != nullptr);
    CVLOG(3, "SWEEPER") << "initial value: " << to_string(this->get_initial_state());
    CVLOG(3, "SWEEPER") << "solution at nodes:";
    for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
      CVLOG(3, "SWEEPER") << "\t" << m << ": " << to_string(this->get_states()[m]);
    }
    CVLOG(3, "SWEEPER") << "solution at t_end: " << to_string(this->get_end_state());
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::advance()
  {
    CLOG(DEBUG, "SWEEPER") << "advancing to next time step";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::spread()
  {
    CLOG(DEBUG, "SWEEPER") << "spreading initial value to all states";

    assert(this->get_initial_state() != nullptr);

    for(size_t m = 1; m < this->get_states().size(); ++m) {
      assert(this->states()[m] != nullptr);
      this->states()[m]->data() = this->get_initial_state()->get_data();
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::save()
  {
    CLOG(DEBUG, "SWEEPER") << "saving states to previous states";

    assert(this->get_quadrature() != nullptr);
    assert(this->get_states().size() == this->get_quadrature()->get_num_nodes() + 1);
    assert(this->get_previous_states().size() == this->get_states().size());

    for (size_t m = 0; m < this->get_states().size(); ++m) {
      CVLOG(2, "SWEEPER") << "\t" << m << ": " << to_string(this->get_states()[m]);
      this->previous_states()[m]->data() = this->get_states()[m]->get_data();
      CVLOG(2, "SWEEPER") << "\t\t-> " << to_string(this->get_previous_states()[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::reevaluate(const bool initial_only)
  {
    UNUSED(initial_only);
    throw NotImplementedYet("reevaluation of right-hand-side");
  }

  template<class SweeperTrait, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_type>>
  Sweeper<SweeperTrait, Enabled>::integrate(const typename SweeperTrait::time_type& dt)
  {
    UNUSED(dt);
    throw NotImplementedYet("integration over dt");
  }

  template<class SweeperTrait, typename Enabled>
  bool
  Sweeper<SweeperTrait, Enabled>::converged()
  {
    CLOG(DEBUG, "SWEEPER") << "convergence check";

    if (this->_abs_residual_tol > 0.0 || this->_rel_residual_tol > 0.0) {
      this->compute_residuals();
      const size_t num_residuals = this->get_residuals().size();
      vector<typename traits::spacial_type>   abs_norms(num_residuals)
                                            , rel_norms(num_residuals);

      assert(this->get_residuals().back() != nullptr);
      abs_norms.back() = this->get_residuals().back()->norm0();

      for (size_t m = 0; m < num_residuals - 1; ++m) {
        assert(this->get_residuals()[m] != nullptr);
        const auto norm = this->get_residuals()[m]->norm0();
        abs_norms[m] = norm;
        rel_norms[m] = abs_norms[m] / abs_norms.back();
      }
      rel_norms.back() = 1.0;

      CVLOG(1, "SWEEPER") << "absolute residuals: " << abs_norms;
      CVLOG(1, "SWEEPER") << "relative residuals: " << rel_norms;

      auto max_abs_norm = *(max_element(abs_norms.cbegin(), abs_norms.cend()));
      auto max_rel_norm = *(max_element(rel_norms.cbegin(), rel_norms.cend()));

      if (max_abs_norm < this->_abs_residual_tol) {
        CLOG(INFO, "SWEEPER") << "Sweeper has converged w.r.t. absolute residual tolerance: "
                              << LOG_FLOAT << max_abs_norm << " < " << this->_abs_residual_tol;
      } else if (max_rel_norm < this->_rel_residual_tol) {
        CLOG(INFO, "SWEEPER") << "Sweeper has converged w.r.t. relative residual tolerance: "
                              << LOG_FLOAT << max_rel_norm << " < " << this->_rel_residual_tol;
      } else {
        CLOG(INFO, "SWEEPER") << "Sweeper has not yet converged to neither residual tolerance.";
      }

      return (   *(max_element(abs_norms.cbegin(), abs_norms.cend())) < this->_abs_residual_tol
              || *(max_element(rel_norms.cbegin(), rel_norms.cend())) < this->_rel_residual_tol);
    } else {
      CLOG(WARNING, "SWEEPER") << "No residual tolerances set. Thus skipping convergence check.";
      return false;
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::integrate_end_state(const typename SweeperTrait::time_type& dt)
  {
    UNUSED(dt);
    assert(this->get_quadrature() != nullptr);
    CLOG(DEBUG, "SWEEPER") << "integrating end state";

    if (this->get_quadrature()->right_is_node()) {
      assert(this->get_end_state() != nullptr);
      assert(this->get_states().size() > 0);

      this->end_state()->data() = this->get_states().back()->get_data();
      CVLOG(1, "SWEEPER") << "end state: " << to_string(this->get_end_state());
    } else {
      throw NotImplementedYet("integration of end state for quadrature not including right time interval boundary");
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::compute_residuals()
  {
    throw NotImplementedYet("computation of residuals");
  }
}  // ::pfasst
