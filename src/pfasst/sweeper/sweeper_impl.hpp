#include "pfasst/sweeper/sweeper.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
using namespace std;

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
      , _logger_id("SWEEPER")
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
      CLOG(ERROR, this->get_logger_id()) << "Sweeper need to be setup first before querying initial state.";
      throw runtime_error("sweeper not setup before querying initial state");
    }
    return this->states().front();
  }

  template<class SweeperTrait, typename Enabled>
  const shared_ptr<typename SweeperTrait::encap_type>
  Sweeper<SweeperTrait, Enabled>::get_initial_state() const
  {
    if (this->get_states().size() == 0) {
      CLOG(ERROR, this->get_logger_id()) << "Sweeper need to be setup first before querying initial state.";
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
  Sweeper<SweeperTrait, Enabled>::set_logger_id(const string& logger_id)
  {
    this->_logger_id = logger_id;
  }

  template<class SweeperTrait, typename Enabled>
  const char*
  Sweeper<SweeperTrait, Enabled>::get_logger_id() const
  {
    return this->_logger_id.c_str();
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::set_options()
  {
    CLOG(DEBUG, this->get_logger_id()) << "setting options";
    this->_abs_residual_tol = config::get_value<typename traits::spacial_type>("abs_res_tol", this->_abs_residual_tol);
    this->_rel_residual_tol = config::get_value<typename traits::spacial_type>("rel_res_tol", this->_rel_residual_tol);
    CVLOG(1, this->get_logger_id()) << "absolut residual tolerance:  " << this->_abs_residual_tol;
    CVLOG(1, this->get_logger_id()) << "relative residual tolerance: " << this->_rel_residual_tol;
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
    CLOG(DEBUG, this->get_logger_id()) << "setting up sweeper with " << to_string(this->get_status());

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
    CLOG(DEBUG, this->get_logger_id()) << "pre-predicting";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::predict()
  {
    CLOG(DEBUG, this->get_logger_id()) << "predicting";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::post_predict()
  {
    CLOG(DEBUG, this->get_logger_id()) << "post-predicting";

    assert(this->get_status() != nullptr);
    this->integrate_end_state(this->get_status()->get_dt());

    assert(this->get_quadrature() != nullptr);
    CVLOG(3, this->get_logger_id()) << "solution at nodes:";
    for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
      CVLOG(3, this->get_logger_id()) << "  " << m << ": " << to_string(this->get_states()[m]);
    }
    CVLOG(3, this->get_logger_id()) << "solution at t_end: " << to_string(this->get_end_state());
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::pre_sweep()
  {
    CLOG(DEBUG, this->get_logger_id()) << "pre-sweeping";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::sweep()
  {
    CLOG(DEBUG, this->get_logger_id()) << "sweeping";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::post_sweep()
  {
    CLOG(DEBUG, this->get_logger_id()) << "post-sweeping";

    assert(this->get_status() != nullptr);
    this->integrate_end_state(this->get_status()->get_dt());

    assert(this->get_quadrature() != nullptr);
    CVLOG(3, this->get_logger_id()) << "solution at nodes:";
    for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
      CVLOG(3, this->get_logger_id()) << "\t" << m << ": " << to_string(this->get_states()[m]);
    }
    CVLOG(3, this->get_logger_id()) << "solution at t_end:" << to_string(this->get_end_state());
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::post_step()
  {
    CLOG(DEBUG, this->get_logger_id()) << "post step";

    assert(this->get_quadrature() != nullptr);
    CVLOG(3, this->get_logger_id()) << "initial value: " << to_string(this->get_initial_state());
    CVLOG(3, this->get_logger_id()) << "solution at nodes:";
    for (size_t m = 0; m <= this->get_quadrature()->get_num_nodes(); ++m) {
      CVLOG(3, this->get_logger_id()) << "\t" << m << ": " << to_string(this->get_states()[m]);
    }
    CVLOG(3, this->get_logger_id()) << "solution at t_end: " << to_string(this->get_end_state());
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::advance(const size_t& num_steps)
  {
    CLOG(DEBUG, this->get_logger_id()) << "advancing " << num_steps << " time steps";
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::spread()
  {
    CLOG(DEBUG, this->get_logger_id()) << "spreading initial value to all states";

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
    CLOG(DEBUG, this->get_logger_id()) << "saving states to previous states";

    assert(this->get_quadrature() != nullptr);
    assert(this->get_states().size() == this->get_quadrature()->get_num_nodes() + 1);
    assert(this->get_previous_states().size() == this->get_states().size());

    for (size_t m = 0; m < this->get_states().size(); ++m) {
      CVLOG(2, this->get_logger_id()) << "\t" << m << ": " << to_string(this->get_states()[m]);
      this->previous_states()[m]->data() = this->get_states()[m]->get_data();
      CVLOG(2, this->get_logger_id()) << "\t\t-> " << to_string(this->get_previous_states()[m]);
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::reevaluate(const bool initial_only)
  {
    UNUSED(initial_only);
    throw runtime_error("reevaluation of right-hand-side");
  }

  template<class SweeperTrait, typename Enabled>
  vector<shared_ptr<typename SweeperTrait::encap_type>>
  Sweeper<SweeperTrait, Enabled>::integrate(const typename SweeperTrait::time_type& dt)
  {
    UNUSED(dt);
    throw runtime_error("integration over dt");
  }

  template<class SweeperTrait, typename Enabled>
  bool
  Sweeper<SweeperTrait, Enabled>::converged()
  {
    CLOG(DEBUG, this->get_logger_id()) << "convergence check";

    this->compute_residuals();
    const size_t num_residuals = this->get_residuals().size();
    this->_abs_res_norms.resize(num_residuals);
    this->_rel_res_norms.resize(num_residuals);

    assert(this->get_residuals().back() != nullptr);
    this->_abs_res_norms.back() = this->get_residuals().back()->norm0();

    for (size_t m = 0; m < num_residuals; ++m) {
      assert(this->get_residuals()[m] != nullptr);
      const auto norm = this->get_residuals()[m]->norm0();
      this->_abs_res_norms[m] = norm;
      this->_rel_res_norms[m] = this->_abs_res_norms[m] / this->get_states()[m]->norm0();
    }

    CVLOG(1, this->get_logger_id()) << "absolute residuals: " << this->_abs_res_norms;
    CVLOG(1, this->get_logger_id()) << "relative residuals: " << this->_rel_res_norms;

    this->status()->abs_res_norm() = *(max_element(this->_abs_res_norms.cbegin(), this->_abs_res_norms.cend()));
    this->status()->rel_res_norm() = *(max_element(this->_rel_res_norms.cbegin(), this->_rel_res_norms.cend()));

    if (this->status()->abs_res_norm() < this->_abs_residual_tol) {
      CLOG(INFO, this->get_logger_id()) << "Sweeper has converged w.r.t. absolute residual tolerance: " << LOG_FLOAT
                                        << this->status()->abs_res_norm() << " < " << this->_abs_residual_tol;
    } else if (this->status()->rel_res_norm() < this->_rel_residual_tol) {
      CLOG(INFO, this->get_logger_id()) << "Sweeper has converged w.r.t. relative residual tolerance: " << LOG_FLOAT
                                        << this->status()->rel_res_norm() << " < " << this->_rel_residual_tol;
    } else {
      CLOG(INFO, this->get_logger_id()) << "Sweeper has not yet converged to neither residual tolerance.";
    }

    if (this->_abs_residual_tol > 0.0 || this->_rel_residual_tol > 0.0) {
      return (   this->status()->abs_res_norm() < this->_abs_residual_tol
              || this->status()->rel_res_norm() < this->_rel_residual_tol);
    } else {
      CLOG(WARNING, this->get_logger_id()) << "No residual tolerances set. Thus skipping convergence check.";
      return false;
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::integrate_end_state(const typename SweeperTrait::time_type& dt)
  {
    UNUSED(dt);
    assert(this->get_quadrature() != nullptr);
    CLOG(DEBUG, this->get_logger_id()) << "integrating end state";

    if (this->get_quadrature()->right_is_node()) {
      assert(this->get_end_state() != nullptr);
      assert(this->get_states().size() > 0);

      this->end_state()->data() = this->get_states().back()->get_data();
      CVLOG(1, this->get_logger_id()) << "end state: " << to_string(this->get_end_state());
    } else {
      throw runtime_error("integration of end state for quadrature not including right time interval boundary");
    }
  }

  template<class SweeperTrait, typename Enabled>
  void
  Sweeper<SweeperTrait, Enabled>::compute_residuals()
  {
    throw runtime_error("computation of residuals");
  }
}  // ::pfasst
