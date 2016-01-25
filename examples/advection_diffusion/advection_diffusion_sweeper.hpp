/**
 * @defgroup AdvectionDiffusionFiles Files
 * @ingroup AdvectionDiffusion
 *
 * @file examples/advection_diffusion/advection_diffusion_sweeper.hpp
 * @since v0.1.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__ADVECTION_DIFFUSION_SWEEPER_HPP_
#define _EXAMPLES__ADVEC_DIFF__ADVECTION_DIFFUSION_SWEEPER_HPP_

#include <cassert>
#include <complex>
#include <cstdlib>
#include <map>
#include <ostream>
#include <vector>
#include <utility>
using namespace std;

#include <pfasst/globals.hpp>
#include <pfasst/config.hpp>
#include <pfasst/logging.hpp>
#include <pfasst/encap/imex_sweeper.hpp>
using pfasst::encap::Encapsulation;
using pfasst::encap::as_vector;

#include "fft.hpp"

#ifndef PI
#define PI 3.1415926535897932385
#endif


namespace pfasst
{
  namespace examples
  {
    /**
     * Advection-Diffusion example
     *
     * @defgroup AdvectionDiffusion Advection Diffusion
     * @ingroup Examples
     *
     * This directory contains several implementations of an advection/diffusion solver using the
     * PFASST framework.
     *
     * All of the solvers use the SDC sweeper defined in `advection_diffusion_sweeper.hpp`, and the FFT
     * routines in `fft.hpp`.
     *
     * The implementations are, in order of complexity:
     *
     *  - `vanilla_sdc.cpp` - basic example that uses an encapsulated IMEX sweeper.
     *
     *  - `serial_mlsdc.cpp` - basic multi-level version that uses polynomial interpolation in time and
     *    spectral interpolation in space, as defined in `specrtal_transfer_1d.hpp`.
     *
     *  - `serial_mlsdc_autobuild.cpp` - same as above, but uses the "auto build" feature to shorten
     *    `main`.
     */
    namespace advection_diffusion
    {
      /**
       * @name Containers for errors/residuals etc.
       * @{
       */
      typedef map<tuple<size_t, size_t>, double> error_map; // step, iteration -> error
      typedef map<size_t, error_map> residual_map; // level, (step, iteration) -> residual
      //! @}

      typedef error_map::value_type vtype;
      typedef error_map::key_type ktype;

      /**
       * advection-diffusion sweeper with semi-implicit time-integration.
       * @ingroup AdvectionDiffusion
       */
      template<typename time = pfasst::time_precision>
      class AdvectionDiffusionSweeper
        : public encap::IMEXSweeper<time>
      {
        public:
          static void init_opts()
          {
            pfasst::config::options::add_option<size_t>("Adv/Diff Sweeper", "spatial_dofs", "Number of spatial degrees of freedom");
          }

          static void init_logs()
          {
            pfasst::log::add_custom_logger("Advec");
          }

        private:
          //! @{
          FFT fft;
          vector<complex<double>> ddx, lap;
          //! @}

          //! @{
          error_map errors;
          error_map residuals;
          //! @}

          //! @{
          double v  = 1.0;
          time   t0 = 1.0;
          double nu = 0.02;
          size_t nf1evals = 0;
          //! @}

        public:
          //! @{
          explicit AdvectionDiffusionSweeper(size_t nvars)
          {
            this->ddx.resize(nvars);
            this->lap.resize(nvars);
            for (size_t i = 0; i < nvars; i++) {
              double kx = 2 * PI * ((i <= nvars / 2) ? int(i) : int(i) - int(nvars));
              this->ddx[i] = complex<double>(0.0, 1.0) * kx;
              this->lap[i] = (kx * kx < 1e-13) ? 0.0 : -kx * kx;
            }
          }

          AdvectionDiffusionSweeper() = default;

          virtual ~AdvectionDiffusionSweeper()
          {
            ML_CLOG(INFO, "Advec", "number of f1 evals: " << this->nf1evals);
          }
          //! @}

          //! @{
          void exact(shared_ptr<Encapsulation<time>> q, time t)
          {
            this->exact(as_vector<double, time>(q), t);
          }

          void exact(DVectorT& q, time t)
          {
            size_t n = q.size();
            double a = 1.0 / sqrt(4 * PI * nu * (t + t0));

            for (size_t i = 0; i < n; i++) {
              q[i] = 0.0;
            }

            for (int ii = -2; ii < 3; ii++) {
              for (size_t i = 0; i < n; i++) {
                double x = double(i) / n - 0.5 + ii - t * v;
                q[i] += a * exp(-x * x / (4 * nu * (t + t0)));
              }
            }
          }

          void echo_error(time t)
          {
            auto& qend = as_vector<double, time>(this->get_end_state());
            DVectorT qex(qend.size());

            this->exact(qex, t);

            double max = 0.0;
            for (size_t i = 0; i < qend.size(); i++) {
              double d = abs(qend[i] - qex[i]);
              if (d > max) { max = d; }
            }

            auto n = this->get_controller()->get_step();
            auto k = this->get_controller()->get_iteration();

            this->errors.insert(vtype(ktype(n, k), max));
          }

          void echo_residual()
          {
            vector<shared_ptr<Encapsulation<time>>> residuals;

            for (size_t m = 0; m < this->get_nodes().size(); m++) {
                residuals.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
            this->residual(this->get_controller()->get_step_size(), residuals);

            vector<time> rnorms;
            for (auto r: residuals) {
              rnorms.push_back(r->norm0());
            }
            auto rmax = *std::max_element(rnorms.begin(), rnorms.end());

            auto n = this->get_controller()->get_step();
            auto k = this->get_controller()->get_iteration();

            auto err = this->errors[ktype(n, k)];

            ML_CLOG(INFO, "Advec", (boost::format(this->FORMAT_STR) % (n+1) % k % this->get_nodes().size() % as_vector<double, time>(this->state[0]).size() % rmax % err));

            this->residuals[ktype(n, k)] = rmax;
          }

          error_map get_errors()
          {
            return this->errors;
          }

          error_map get_residuals()
          {
            return this->residuals;
          }
          //! @}

          //! @{
          /**
           * @copybrief pfasst::encap::IMEXSweeper::predict()
           */
          void post_predict() override
          {
            time t  = this->get_controller()->get_time();
            time dt = this->get_controller()->get_step_size();
            this->echo_error(t + dt);
            this->echo_residual();
          }

          /**
           * @copybrief pfasst::encap::IMEXSweeper::sweep()
           */
          void post_sweep() override
          {
            time t  = this->get_controller()->get_time();
            time dt = this->get_controller()->get_step_size();
            this->echo_error(t + dt);
            this->echo_residual();
          }
          //! @}

          //! @{
          /**
           * @copybrief pfasst::encap::IMEXSweeper::f_expl_eval()
           */
          void f_expl_eval(shared_ptr<Encapsulation<time>> f_expl_encap,
                           shared_ptr<Encapsulation<time>> u_encap,
                           time t) override
          {
            UNUSED(t);
            auto& u = as_vector<double, time>(u_encap);
            auto& f_expl = as_vector<double, time>(f_expl_encap);

            double c = -v / double(u.size());

            auto* z = this->fft.forward(u);
            for (size_t i = 0; i < u.size(); i++) {
              z[i] *= c * this->ddx[i];
            }
            this->fft.backward(f_expl);

            this->nf1evals++;
          }

          /**
           * @copybrief pfasst::encap::IMEXSweeper::f_impl_eval()
           */
          void f_impl_eval(shared_ptr<Encapsulation<time>> f_impl_encap,
                           shared_ptr<Encapsulation<time>> u_encap,
                           time t) override
          {
            UNUSED(t);
            auto& u = as_vector<double, time>(u_encap);
            auto& f_impl = as_vector<double, time>(f_impl_encap);

            double c = nu / double(u.size());

            auto* z = this->fft.forward(u);
            for (size_t i = 0; i < u.size(); i++) {
              z[i] *= c * this->lap[i];
            }
            this->fft.backward(f_impl);
          }

          /**
           * @copybrief pfasst::encap::IMEXSweeper::impl_solve()
           */
          void impl_solve(shared_ptr<Encapsulation<time>> f_impl_encap,
                          shared_ptr<Encapsulation<time>> u_encap,
                          time t, time dt,
                          shared_ptr<Encapsulation<time>> rhs_encap) override
          {
            UNUSED(t);
            auto& u = as_vector<double, time>(u_encap);
            auto& f_impl = as_vector<double, time>(f_impl_encap);
            auto& rhs = as_vector<double, time>(rhs_encap);

            double c = nu * double(dt);

            auto* z = this->fft.forward(rhs);
            for (size_t i = 0; i < u.size(); i++) {
              z[i] /= (1.0 - c * this->lap[i]) * double(u.size());
            }
            this->fft.backward(u);

            for (size_t i = 0; i < u.size(); i++) {
              f_impl[i] = (u[i] - rhs[i]) / double(dt);
            }
          }
          //! @}
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__ADVEC_DIFF__ADVECTION_DIFFUSION_SWEEPER_HPP_
