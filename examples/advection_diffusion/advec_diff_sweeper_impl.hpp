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
#include <pfasst/logging.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace advec_diff
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
      void
      AdvecDiff<SweeperTrait, Enabled>::post_sweep()
      {
        IMEX<SweeperTrait, Enabled>::post_sweep();

        assert(this->get_quadrature() != nullptr);
        const time_type t = this->get_status()->get_time();

        this->compute_residuals();
        auto error = this->compute_error(this->get_end_state(), t);

        CLOG(INFO, "SWEEPER") << "at t_end:";
        CLOG(INFO, "SWEEPER") << "  norm_0(residual): " << encap::norm0(this->get_residuals().back());

        CLOG(INFO, "SWEEPER") << "  norm_0(error):    " << encap::norm0(error);
      }

      template<class SweeperTrait, typename Enabled>
      void
      AdvecDiff<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        CLOG(INFO, "SWEEPER") << "number function evaluations:";
        CLOG(INFO, "SWEEPER") << "  expl: " << this->num_expl_f_evals;
        CLOG(INFO, "SWEEPER") << "  impl: " << this->num_impl_f_evals;

        this->num_expl_f_evals = 0;
        this->num_impl_f_evals = 0;
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
        CVLOG(2, "SWEEPER") << "evaluating EXPLICIT part at t=" << t;
        CVLOG(5, "SWEEPER") << "\tu:   " << to_string(u);

        time_type c = (-1.0 * this->_v) / time_type(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_ddx[i];
        }

        auto result = this->get_encap_factory()->create();
        this->_fft.backward(result);

        this->num_expl_f_evals++;

        CVLOG(5, "SWEEPER") << "\t  -> " << to_string(result);
        return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      AdvecDiff<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                          const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        CVLOG(2, "SWEEPER") << "evaluating IMPLICIT part at t=" << t;
        CVLOG(5, "SWEEPER") << "\tu:   " << to_string(u);

        time_type c = this->_nu / time_type(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_lap[i];
        }

        auto result = this->get_encap_factory()->create();
        this->_fft.backward(result);

        this->num_impl_f_evals++;

        CVLOG(5, "SWEEPER") << "\t  -> " << to_string(result);
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
        CVLOG(2, "SWEEPER") << "implicit spacial solve at t=" << t << " with dt=" << dt;
        CVLOG(5, "SWEEPER") << "\tf:   " << to_string(f);
        CVLOG(5, "SWEEPER") << "\tu:   " << to_string(u);
        CVLOG(5, "SWEEPER") << "\trhs: " << to_string(rhs);

        time_type c = this->_nu * dt;

        auto* z = this->_fft.forward(rhs);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] /= (1.0 - c * this->_lap[i]) * time_type(this->get_num_dofs());
        }
        this->_fft.backward(u);

        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          f->data()[i] = (u->data()[i] - rhs->data()[i]) / dt;
        }

        CVLOG(5, "SWEEPER") << "\t->";
        CVLOG(5, "SWEEPER") << "\t  f: " << to_string(f);
        CVLOG(5, "SWEEPER") << "\t  u: " << to_string(u);
      }
    }  // ::pfasst::examples::advec_diff
  }  // ::pfasst::examples
}  // ::pfasst
