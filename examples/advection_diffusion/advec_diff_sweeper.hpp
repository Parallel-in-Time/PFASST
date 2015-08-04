#ifndef _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_
#define _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/contrib/fft.hpp>


namespace pfasst
{
  namespace examples
  {
    template<
      class SweeperTrait,
      typename Enabled = void
    >
    class AdvecDiff
      : public IMEX<SweeperTrait, Enabled>
    {
      static_assert(is_same<
                      vector<typename SweeperTrait::time_type>,
                      typename SweeperTrait::encap_type::data_type
                    >::value,
                    "Advection Diffusion Sweeper requires encapsulated vectors");

      public:
        typedef          SweeperTrait         traits;
        typedef typename traits::encap_type   encap_type;
        typedef typename traits::time_type    time_type;
        typedef typename traits::spacial_type spacial_type;

      private:
        spacial_type _v;
        time_type    _t0;
        spacial_type _nu;

        pfasst::contrib::FFT<time_type> _fft;
        vector<complex<time_type>>      _ddx;
        vector<complex<time_type>>      _lap;

      protected:
        virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                                                const shared_ptr<typename SweeperTrait::encap_type> u);
        virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                                                const shared_ptr<typename SweeperTrait::encap_type> u);

        virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                    shared_ptr<typename SweeperTrait::encap_type> u,
                                    const typename SweeperTrait::time_type& t,
                                    const typename SweeperTrait::time_type& dt,
                                    const shared_ptr<typename SweeperTrait::encap_type> rhs);

        virtual shared_ptr<typename SweeperTrait::encap_type> compute_error(const shared_ptr<typename SweeperTrait::encap_type> q,
                                                                            const typename SweeperTrait::time_type& t);

      public:
        explicit AdvecDiff(const size_t& ndofs);
        AdvecDiff(const AdvecDiff<SweeperTrait, Enabled>& other) = default;
        AdvecDiff(AdvecDiff<SweeperTrait, Enabled>&& other) = default;
        virtual ~AdvecDiff() = default;
        AdvecDiff<SweeperTrait, Enabled>& operator=(const AdvecDiff<SweeperTrait, Enabled>& other) = default;
        AdvecDiff<SweeperTrait, Enabled>& operator=(AdvecDiff<SweeperTrait, Enabled>&& other) = default;

        virtual shared_ptr<typename SweeperTrait::encap_type> exact(const typename SweeperTrait::time_type& t);

        virtual size_t get_num_dofs() const;
    };
  }  // ::pfasst::examples
}  // ::pfasst

#include "advec_diff_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__ADVEC_DIFF__ADVEC_DIFF_SWEEPER_HPP_
