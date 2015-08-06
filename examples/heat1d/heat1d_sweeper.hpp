#ifndef _PFASST__EXAMPLES__HEAD1D__HEAD1D_SWEEPER_HPP_
#define _PFASST__EXAMPLES__HEAD1D__HEAD1D_SWEEPER_HPP_

#include <memory>
#include <vector>
using namespace std;

#include <pfasst/sweeper/imex.hpp>
#include <pfasst/contrib/fft.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace heat1d
    {
      template<
        class SweeperTrait,
        typename Enabled = void
      >
      class Heat1D
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
          time_type    _t0;
          spacial_type _nu;

          pfasst::contrib::FFT<spacial_type> _fft;
          vector<complex<spacial_type>>      _lap;

        protected:
          virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                                                  const shared_ptr<typename SweeperTrait::encap_type> u) override;
          virtual shared_ptr<typename SweeperTrait::encap_type> evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                                                  const shared_ptr<typename SweeperTrait::encap_type> u) override;

          virtual void implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                      shared_ptr<typename SweeperTrait::encap_type> u,
                                      const typename SweeperTrait::time_type& t,
                                      const typename SweeperTrait::time_type& dt,
                                      const shared_ptr<typename SweeperTrait::encap_type> rhs) override;

          virtual shared_ptr<typename SweeperTrait::encap_type> compute_error(const shared_ptr<typename SweeperTrait::encap_type> q,
                                                                              const typename SweeperTrait::time_type& t);

        public:
          explicit Heat1D(const size_t& ndofs);
          Heat1D(const Heat1D<SweeperTrait, Enabled>& other) = default;
          Heat1D(Heat1D<SweeperTrait, Enabled>&& other) = default;
          virtual ~Heat1D() = default;
          Heat1D<SweeperTrait, Enabled>& operator=(const Heat1D<SweeperTrait, Enabled>& other) = default;
          Heat1D<SweeperTrait, Enabled>& operator=(Heat1D<SweeperTrait, Enabled>&& other) = default;

          virtual shared_ptr<typename SweeperTrait::encap_type> exact(const typename SweeperTrait::time_type& t);

          virtual void post_predict() override;
          virtual void post_sweep() override;
          virtual void post_step() override;

          virtual size_t get_num_dofs() const;
      };
    }  // ::pfasst::examples::heat1d
  }  // ::pfasst::examples
}  // ::pfasst

#include "heat1d_sweeper_impl.hpp"

#endif  // _PFASST__EXAMPLES__HEAD1D__HEAD1D_SWEEPER_HPP_
