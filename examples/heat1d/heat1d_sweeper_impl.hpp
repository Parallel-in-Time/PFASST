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
        //   xvalues = np.array([(i) * self.dx for i in range(self.nvars)])
        //   me.values = np.sin(2 * np.pi * xvalues) * np.exp(-t * (2 * np.pi)**2 * self.nu)
        const spacial_type dx = 1.0 / spacial_type(this->get_num_dofs());
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          result->data()[i] = sin(two_pi<spacial_type>() * i * dx) * exp(-t * pow(two_pi<spacial_type>(), 2) * this->_nu);
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
        const time_type dt = this->get_status()->get_dt();

        this->compute_residuals();
        auto error = this->compute_error(this->get_end_state(), t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_type(t));

        for (size_t m = 0; m < num_nodes + 1; ++m) {
          const time_type ds = dt * (nodes[m+1] - nodes[m]);
          CLOG(INFO, "USER") << "t["<<m<<"]=" << LOG_FIXED << (ds * nodes[m]);
          CLOG(INFO, "USER") << "  |residual| = " << LOG_FLOAT << encap::norm0(this->get_residuals()[m]);
          CLOG(INFO, "USER") << "  |error|    = " << LOG_FLOAT << encap::norm0(error[m]);
        }
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat1D<SweeperTrait, Enabled>::post_sweep()
      {
        IMEX<SweeperTrait, Enabled>::post_sweep();

        assert(this->get_status() != nullptr);
        const time_type t = this->get_status()->get_time();
        const time_type dt = this->get_status()->get_dt();

        this->compute_residuals();
        auto error = this->compute_error(this->get_end_state(), t);

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_type(t));

        for (size_t m = 0; m < num_nodes + 1; ++m) {
          const time_type ds = dt * (nodes[m+1] - nodes[0]);
          CLOG(INFO, "USER") << "t["<<m<<"]=" << LOG_FIXED << (ds * nodes[m]);
          CLOG(INFO, "USER") << "  |residual| = " << LOG_FLOAT << encap::norm0(this->get_residuals()[m]);
          CLOG(INFO, "USER") << "  |error|    = " << LOG_FLOAT << encap::norm0(error[m]);
        }
      }

      template<class SweeperTrait, typename Enabled>
      void
      Heat1D<SweeperTrait, Enabled>::post_step()
      {
        IMEX<SweeperTrait, Enabled>::post_step();

        CLOG(INFO, "USER") << "number function evaluations:";
        CLOG(INFO, "USER") << "  expl:        " << this->_num_expl_f_evals;
        CLOG(INFO, "USER") << "  impl:        " << this->_num_impl_f_evals;
        CLOG(INFO, "USER") << "  impl solves: " << this->_num_impl_solves;

        this->_num_expl_f_evals = 0;
        this->_num_impl_f_evals = 0;
        this->_num_impl_solves = 0;
      }

      template<class SweeperTrait, typename Enabled>
      size_t
      Heat1D<SweeperTrait, Enabled>::get_num_dofs() const
      {
        return this->get_encap_factory()->size();
      }


      template<class SweeperTrait, typename Enabled>
      vector<shared_ptr<typename SweeperTrait::encap_type>>
      Heat1D<SweeperTrait, Enabled>::compute_error(const shared_ptr<typename SweeperTrait::encap_type> q,
                                                   const typename SweeperTrait::time_type& t)
      {
        assert(this->get_status() != nullptr);
        const time_type dt = this->get_status()->get_dt();

        assert(this->get_quadrature() != nullptr);
        auto nodes = this->get_quadrature()->get_nodes();
        const auto num_nodes = this->get_quadrature()->get_num_nodes();
        nodes.insert(nodes.begin(), time_type(t));

        vector<shared_ptr<encap_type>> error;
        error.resize(num_nodes + 1);
        generate(error.begin(), error.end(),
                 bind(&encap_type::factory_type::create, this->encap_factory()));

        for (size_t m = 1; m < num_nodes + 1; ++m) {
          const time_type ds = dt * (nodes[m] - nodes[0]);
          error[m] = pfasst::encap::axpy(-1.0, this->exact(t + ds), q);
        }

        return error;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat1D<SweeperTrait, Enabled>::evaluate_rhs_expl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        CVLOG(2, "USER") << "evaluating EXPLICIT part at t=" << t;
        CVLOG(5, "USER") << "\tu:   " << to_string(u);

        auto result = this->get_encap_factory()->create();

        // taken form pySDC
        //   # xvalues = np.array([(i+1)*self.dx for i in range(self.nvars)])
        //   fexpl.values = np.zeros(self.nvars)  # -np.sin(np.pi * xvalues) * (np.sin(t) - self.nu * np.pi**2 * np.cos(t))
//         const spacial_type PI = pi<spacial_type>();
//         const spacial_type PIsqr = pi_sqr<spacial_type>();
//         const spacial_type dx = 1.0 / (spacial_type(this->get_num_dofs()) + 1);
//         for (size_t i = 0; i < this->get_num_dofs(); ++i) {
//           result->data()[i] = -1.0 * sin(PI * (i + 1) * dx) * (sin(t) - this->_nu * PIsqr * cos(t));
//         }
        result->zero();

        this->_num_expl_f_evals++;

        CVLOG(5, "USER") << "\t  -> " << to_string(result);
        return result;
      }

      template<class SweeperTrait, typename Enabled>
      shared_ptr<typename SweeperTrait::encap_type>
      Heat1D<SweeperTrait, Enabled>::evaluate_rhs_impl(const typename SweeperTrait::time_type& t,
                                                       const shared_ptr<typename SweeperTrait::encap_type> u)
      {
        CVLOG(2, "USER") << "evaluating IMPLICIT part at t=" << t;
        CVLOG(5, "USER") << "\tu:   " << to_string(u);

        // taken from Matt's original Advection-Diffusion example from old PFASST++
        spacial_type c = this->_nu / spacial_type(this->get_num_dofs());

        auto* z = this->_fft.forward(u);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] *= c * this->_lap[i];
        }

        auto result = this->get_encap_factory()->create();
        this->_fft.backward(result);

        this->_num_impl_f_evals++;

        CVLOG(5, "USER") << "\t  -> " << to_string(result);
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
        CVLOG(2, "USER") << "IMPLICIT spacial SOLVE at t=" << t << " with dt=" << dt;
        CVLOG(5, "USER") << "\tf:   " << to_string(f);
        CVLOG(5, "USER") << "\tu:   " << to_string(u);
        CVLOG(5, "USER") << "\trhs: " << to_string(rhs);

        spacial_type c = this->_nu * dt;

        auto* z = this->_fft.forward(rhs);
        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          z[i] /= (1.0 - c * this->_lap[i]) * spacial_type(this->get_num_dofs());
        }
        this->_fft.backward(u);

        for (size_t i = 0; i < this->get_num_dofs(); ++i) {
          f->data()[i] = (u->get_data()[i] - rhs->get_data()[i]) / dt;
        }

        this->_num_impl_solves++;

        CVLOG(5, "USER") << "\t->";
        CVLOG(5, "USER") << "\t  f: " << to_string(f);
        CVLOG(5, "USER") << "\t  u: " << to_string(u);
      }
    }  // ::pfasst::examples::advec_diff
  }  // ::pfasst::examples
}  // ::pfasst
