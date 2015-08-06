#include "heat1d_sweeper.hpp"

#include <cmath>
#include <complex>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <boost/math/constants/constants.hpp>
#include <leathers/pop>
using boost::math::constants::pi;
using boost::math::constants::two_pi;
using boost::math::constants::pi_sqr;

#include <pfasst/globals.hpp>
#include <pfasst/util.hpp>
#include <pfasst/logging.hpp>


namespace pfasst
{
  namespace examples
  {
    namespace heat1d
    {
      template<class SweeperTrait, typename Enabled>
      Heat1D<SweeperTrait, Enabled>::Heat1D(const size_t& ndofs)
        :   IMEX<SweeperTrait, Enabled>()
          , _t0(0.0)
          , _nu(0.02)
          , _lap(ndofs)
      {
        this->encap_factory()->set_size(ndofs);

        for (size_t i = 0; i < ndofs; ++i) {
          spacial_type kx = two_pi<spacial_type>()
                            * ((i <= ndofs / 2) ? spacial_type(i)
                                                : spacial_type(i) - spacial_type(ndofs));
          this->_lap[i] = pfasst::almost_zero(kx * kx) ? 0.0 : -kx * kx;
        }
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat1D<SweeperTrait, Enabled>::exact(const typename SweeperTrait::time_type& t)
      {
        auto result = this->get_encap_factory()->create();

        // taken from pySDC
        //   np.sin(2*np.pi*xvalues)*np.exp(-t*(2*np.pi)**2*self.nu)
        const spacial_type dx = 1.0 / (spacial_type(this->get_num_dofs()) + 1);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          result->data()[i] = sin(pi<spacial_type>() * (i + 1) * dx) * exp(-t * pow(two_pi<spacial_type>(), 2) * this->_nu);
        }

        return result;
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat1D<SweeperTrait, Enabled>::post_predict()
      {
        IMEX<SweeperTrait, Enabled>::post_predict();

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
      Heat1D<SweeperTrait, Enabled>::post_sweep()
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
      Heat1D<SweeperTrait, Enabled>::post_step()
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
      Heat1D<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory()->size();
      }


      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat1D<SweeperTrait, Enabled>::compute_error(const shared_ptr<typename SweeperTrait::encap_type> q,
                                                   const typename SweeperTrait::time_type& t)
      {
        return pfasst::encap::axpy(-1.0, this->exact(t), q);
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat1D<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        CVLOG(2, "SWEEPER") << "evaluating EXPLICIT part at t=" << t;
        CVLOG(5, "SWEEPER") << "\tu:   " << to_string(u);

        auto result = this->get_encap_factory()->create();

        // taken form pySDC
        const spacial_type PI = pi<spacial_type>();
        const spacial_type PIsqr = pi_sqr<spacial_type>();
        const spacial_type dx = 1.0 / (spacial_type(this->get_num_dofs()) + 1);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          result->data()[i] = -1.0 * sin(PI * (i+1) * dx) * (sin(t) - this->_nu * PIsqr * cos(t));
        }

        this->num_expl_f_evals++;

        CVLOG(5, "SWEEPER") << "\t  -> " << to_string(result);
        return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat1D<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        CVLOG(2, "SWEEPER") << "evaluating IMPLICIT part at t=" << t;
        CVLOG(5, "SWEEPER") << "\tu:   " << to_string(u);

        // taken from Matt's original Advection-Diffusion example from old PFASST++
        spacial_type c = this->_nu / spacial_type(this->get_num_dofs());

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
      Heat1D<SweeperTrait, Enabled>::implicit_solve(shared_ptr<typename SweeperTrait::encap_type> f,
                                                    shared_ptr<typename SweeperTrait::encap_type> u,
                                                    const typename SweeperTrait::time_type& t,
                                                    const typename SweeperTrait::time_type& dt,
                                                    const shared_ptr<typename SweeperTrait::encap_type> rhs)
      {
        CVLOG(2, "SWEEPER") << "implicit spacial solve at t=" << t << " with dt=" << dt;
        CVLOG(5, "SWEEPER") << "\tf:   " << to_string(f);
        CVLOG(5, "SWEEPER") << "\tu:   " << to_string(u);
        CVLOG(5, "SWEEPER") << "\trhs: " << to_string(rhs);

        spacial_type c = this->_nu * dt;

        auto* z = this->_fft.forward(rhs);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] /= (1.0 - c * this->_lap[i]) * spacial_type(this->get_num_dofs());
        }
        this->_fft.backward(u);

        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          f->data()[i] = (u->get_data()[i] - rhs->get_data()[i]) / dt;
        }

        CVLOG(5, "SWEEPER") << "\t->";
        CVLOG(5, "SWEEPER") << "\t  f: " << to_string(f);
        CVLOG(5, "SWEEPER") << "\t  u: " << to_string(u);
      }
    }  // ::pfasst::examples::advec_diff
  }  // ::pfasst::examples
}  // ::pfasst
