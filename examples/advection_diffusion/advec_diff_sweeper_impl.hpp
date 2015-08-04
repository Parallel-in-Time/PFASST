#include "advec_diff_sweeper.hpp"

#include <cmath>
#include <complex>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <boost/math/constants/constants.hpp>
#include <leathers/pop>

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>


namespace pfasst
{
  namespace examples
  {
    template<class SweeperTrait, typename Enabled>
    AdvecDiff<SweeperTrait, Enabled>::AdvecDiff(const size_t& ndofs)
      :   IMEX<SweeperTrait, Enabled>()
        , _v(1.0)
        , _t0(1.0)
        , _nu(0.02)
        , _ddx(ndofs)
        , _lap(ndofs)
    {
      this->encap_factory()->set_size(ndofs);

      for (size_t i = 0; i < ndofs; ++i) {
        time_type kx = 2 * boost::math::constants::pi<time_type>() * ((i <= ndofs / 2) ? int(i) : int(i) - int(ndofs));
        this->_ddx[i] = complex<time_type>(0.0, 1.0) * kx;
        this->_lap[i] = pfasst::almost_zero(kx * kx) ? 0.0 : -kx * kx;
      }
    }

    template<class SweeperTrait, typename Enabled>
    shared_ptr<typename SweeperTrait::encap_type>
    AdvecDiff<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_type& t)
    {
      spacial_type a = 1.0 / sqrt(4.0 * boost::math::constants::pi<spacial_type>() * this->_nu * (t + this->_t0));

      auto result = this->get_encap_factory()->create();

      for (int ii = -2; ii < 3; ++ii) {
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          spacial_type x = spacial_type(i) / this->get_num_dofs() - 0.5 + ii - t * this->_v;
          result->data()[i] += a * exp(-x * x / (4.0 * this->_nu * (t + this->_t0)));
        }
      }

      return result;
    }

    template<class SweeperTrait, typename Enabled>
    size_t
    AdvecDiff<SweeperTrait, Enabled>::get_num_dofs() const
    {
      return this->get_encap_factory()->size();
    }


    template<class SweeperTrait, typename Enabled>
    shared_ptr<typename SweeperTrait::encap_type>
    AdvecDiff<SweeperTrait, Enabled>::compute_error(const shared_ptr<typename SweeperTrait::encap_type> q,
                                                    const typename SweeperTrait::time_type& t)
    {
      return pfasst::encap::axpy(-1.0, this->exact(t), q);
    }

    template<class SweeperTrait, typename Enabled>
    shared_ptr<typename SweeperTrait::encap_type>
    AdvecDiff<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                        const shared_ptr<typename SweeperTrait::encap_type> u)
    {
      UNUSED(t);

      time_type c = (-1.0 * this->_v) / time_type(this->get_num_dofs());

      auto* z = this->_fft.forward(u);
      for (size_t i = 0; i < this->get_num_dofs(); ++i) {
        z[i] *= c * this->_ddx[i];
      }

      auto result = this->get_encap_factory()->create();
      this->_fft.backward(result);

      this->num_expl_f_evals++;

      return result;
    }

    template<class SweeperTrait, typename Enabled>
    shared_ptr<typename SweeperTrait::encap_type>
    AdvecDiff<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                        const shared_ptr<typename SweeperTrait::encap_type> u)
    {
      UNUSED(t);

      time_type c = this->_nu / time_type(this->get_num_dofs());

      auto* z = this->_fft.forward(u);
      for (size_t i = 0; i < this->get_num_dofs(); ++i) {
        z[i] *= c * this->_lap[i];
      }

      auto result = this->get_encap_factory()->create();
      this->_fft.backward(result);

      this->num_impl_f_evals++;

      return result;
    }

    template<class SweeperTrait, typename Enabled>
    void
    AdvecDiff<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                                     shared_ptr<typename SweeperTrait::encap_type> u,
                                                     const typename SweeperTrait::time_type& t,
                                                     const typename SweeperTrait::time_type& dt,
                                                     const shared_ptr<typename SweeperTrait::encap_type> rhs)
    {
      UNUSED(t);

      time_type c = this->_nu * dt;

      auto* z = this->_fft.forward(rhs);
      for (size_t i = 0; i < this->get_num_dofs(); ++i) {
        z[i] /= (1.0 - c * this->_lap[i]) * time_type(this->get_num_dofs());
      }
      this->_fft.backward(u);

      for (size_t i = 0; i < this->get_num_dofs(); ++i) {
        f->data()[i] = (u->data()[i] - rhs->data()[i]) / dt;
      }
    }
  }
}
