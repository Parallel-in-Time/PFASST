#ifndef _PFASST__SWEEPER__IMEX_HPP_
#define _PFASST__SWEEPER__IMEX_HPP_

#include "pfasst/sweeper/sweeper.hpp"


namespace pfasst
{
  template<
    class SweeperTrait,
    typename Enabled = void
  >
  class IMEX
    : public Sweeper<SweeperTrait, Enabled>
  {
    public:
      typedef          SweeperTrait                 sweeper_traits;
      typedef typename sweeper_traits::encap_type   encap_type;
      typedef typename sweeper_traits::time_type    time_type;
      typedef typename sweeper_traits::spacial_type spacial_type;

    protected:
      Matrix<time_type> _q_delta_expl;
      Matrix<time_type> _q_delta_impl;

      //! size = #nodes + 1
      vector<shared_ptr<encap_type>> _q_integrals;
      vector<shared_ptr<encap_type>> _expl_rhs;
      vector<shared_ptr<encap_type>> _impl_rhs;

      size_t _num_expl_f_evals;
      size_t _num_impl_f_evals;
      size_t _num_impl_solves;

      virtual void integrate_end_state(const typename SweeperTrait::time_type& dt) override;
      virtual void compute_residuals() override;

      virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                                              const shared_ptr<typename SweeperTrait::encap_type> u);
      virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                                              const shared_ptr<typename SweeperTrait::encap_type> u);

      virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                  shared_ptr<typename SweeperTrait::encap_type> u,
                                  const typename SweeperTrait::time_type& t,
                                  const typename SweeperTrait::time_type& dt,
                                  const shared_ptr<typename SweeperTrait::encap_type> rhs);

      virtual void compute_delta_matrices();

    public:
      explicit IMEX();
      IMEX(const IMEX<SweeperTrait, Enabled>& other) = default;
      IMEX(IMEX<SweeperTrait, Enabled>&& other) = default;
      virtual ~IMEX() = default;
      IMEX<SweeperTrait, Enabled>& operator=(const IMEX<SweeperTrait, Enabled>& other) = default;
      IMEX<SweeperTrait, Enabled>& operator=(IMEX<SweeperTrait, Enabled>&& other) = default;

      virtual void setup() override;

      virtual void pre_predict() override;
      /**
       * first order Euler method
       */
      virtual void predict() override;
      virtual void post_predict() override;

      virtual void pre_sweep() override;
      virtual void sweep() override;
      virtual void post_sweep() override;

      virtual void post_step() override;
      virtual void advance(const size_t& num_steps = 1) override;
      virtual void reevaluate(const bool initial_only = false) override;
      virtual vector<shared_ptr<typename SweeperTrait::encap_type>> integrate(const typename SweeperTrait::time_type& dt) override;
  };
}  // ::pfasst

#include "pfasst/sweeper/imex_impl.hpp"

#endif  // _PFASST__SWEEPER__IMEX_HPP_
