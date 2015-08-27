#ifndef _PFASST__SWEEPER__INTERFACE_HPP_
#define _PFASST__SWEEPER__INTERFACE_HPP_

#include <memory>
#include <type_traits>
#include <vector>
using namespace std;

#include "pfasst/sweeper/traits.hpp"
#include "pfasst/controller/status.hpp"
#include "pfasst/encap/encapsulation.hpp"
#include "pfasst/quadrature.hpp"
using pfasst::quadrature::IQuadrature;


namespace pfasst
{
  template<
    class SweeperTrait,
    typename Enabled = void
  >
  class Sweeper
    : public enable_shared_from_this<Sweeper<SweeperTrait, Enabled>>
  {
    public:
      typedef          SweeperTrait traits;
      // TODO: get rid of mirroring and query trait directly in code
      typedef typename traits::encap_type encap_type;

    static_assert(is_arithmetic<typename traits::time_type>::value,
                  "precision type must be arithmetic");
    static_assert(is_constructible<typename traits::encap_type>::value,
                  "Encapsulation type must be constructible");
    static_assert(is_destructible<typename traits::encap_type>::value,
                  "Encapsulation type must be destructible");

    protected:
      shared_ptr<IQuadrature<typename traits::time_type>>    _quadrature;
      shared_ptr<typename traits::encap_type::factory_type>  _factory;

      //! size = #nodes + 1
      vector<shared_ptr<typename traits::encap_type>>        _states;
      vector<shared_ptr<typename traits::encap_type>>        _previous_states;
      shared_ptr<typename traits::encap_type>                _end_state;

      //! size = #nodes + 1
      vector<shared_ptr<typename traits::encap_type>>        _tau;
      vector<shared_ptr<typename traits::encap_type>>        _residuals;
      vector<typename traits::spacial_type>                  _abs_res_norms;
      vector<typename traits::spacial_type>                  _rel_res_norms;

      shared_ptr<Status<typename traits::time_type>>         _status;
      typename traits::spacial_type                          _abs_residual_tol;
      typename traits::spacial_type                          _rel_residual_tol;

      string                                                 _logger_id;

      virtual void integrate_end_state(const typename SweeperTrait::time_type& dt);
      virtual void compute_residuals();

      virtual vector<shared_ptr<typename SweeperTrait::encap_type>>& previous_states();
      virtual shared_ptr<typename SweeperTrait::encap_type>&         end_state();
      virtual vector<shared_ptr<typename SweeperTrait::encap_type>>& residuals();

    public:
      explicit Sweeper();
      Sweeper(const Sweeper<SweeperTrait, Enabled>& other) = default;
      Sweeper(Sweeper<SweeperTrait, Enabled>&& other) = default;
      virtual ~Sweeper() = default;
      Sweeper<SweeperTrait, Enabled>& operator=(const Sweeper<SweeperTrait, Enabled>& other) = default;
      Sweeper<SweeperTrait, Enabled>& operator=(Sweeper<SweeperTrait, Enabled>&& other) = default;

      virtual       shared_ptr<IQuadrature<typename SweeperTrait::time_type>>& quadrature();
      virtual const shared_ptr<IQuadrature<typename SweeperTrait::time_type>> get_quadrature() const;

      virtual       shared_ptr<Status<typename SweeperTrait::time_type>>& status();
      virtual const shared_ptr<Status<typename SweeperTrait::time_type>>  get_status() const;

      virtual       shared_ptr<typename SweeperTrait::encap_type::factory_type>& encap_factory();
      virtual const shared_ptr<typename SweeperTrait::encap_type::factory_type>  get_encap_factory() const;

      virtual       shared_ptr<typename SweeperTrait::encap_type>&         initial_state();
      virtual       vector<shared_ptr<typename SweeperTrait::encap_type>>& states();
      virtual       vector<shared_ptr<typename SweeperTrait::encap_type>>& tau();

      virtual const shared_ptr<typename SweeperTrait::encap_type>          get_initial_state() const;
      virtual const vector<shared_ptr<typename SweeperTrait::encap_type>>& get_states() const;
      virtual const vector<shared_ptr<typename SweeperTrait::encap_type>>& get_previous_states() const;
      virtual const shared_ptr<typename SweeperTrait::encap_type>          get_end_state() const;
      virtual const vector<shared_ptr<typename SweeperTrait::encap_type>>& get_tau() const;
      virtual const vector<shared_ptr<typename SweeperTrait::encap_type>>& get_residuals() const;

      virtual       void  set_logger_id(const string& logger_id);
      virtual const char* get_logger_id() const;

      virtual void set_options();
      virtual void set_abs_residual_tol(const typename SweeperTrait::spacial_type& abs_res_tol);
      virtual void set_rel_residual_tol(const typename SweeperTrait::spacial_type& rel_res_tol);

      virtual void setup();

      virtual void pre_predict();
      virtual void predict();
      virtual void post_predict();

      virtual void pre_sweep();
      virtual void sweep();
      virtual void post_sweep();

      virtual void post_step();

      virtual void advance(const size_t& num_steps = 1);

      virtual void spread();
      virtual void save();
      virtual void reevaluate(const bool initial_only=false);
      virtual vector<shared_ptr<typename SweeperTrait::encap_type>> integrate(const typename SweeperTrait::time_type& dt);

      virtual bool converged();
  };
}

#include "pfasst/sweeper/sweeper_impl.hpp"

#endif  // _PFASST__SWEEPER__INTERFACE_HPP_
