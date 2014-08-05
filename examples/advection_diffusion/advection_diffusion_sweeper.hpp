/*
 * Advection/diffusion sweeper.
 */

#ifndef _ADVECTION_DIFFUSION_SWEEPER_HPP_
#define _ADVECTION_DIFFUSION_SWEEPER_HPP_

#include <cstdlib>
#include <cassert>
#include <complex>
#include <map>
#include <ostream>
#include <vector>
#include <utility>

#include <pfasst/globals.hpp>
#include <pfasst/encap/imex_sweeper.hpp>

#include "fft.hpp"

#ifndef PI
#define PI 3.1415926535897932385
#endif

using namespace std;
using pfasst::encap::Encapsulation;
using pfasst::encap::as_vector;


/**
 * errors at different iterations and time nodes
 * 
 * Mapping a pair of step/iteration indices onto the error of the solution.
 */
typedef map<pair<size_t, size_t>, double> error_map;

template<typename time = pfasst::time_precision>
class AdvectionDiffusionSweeper
  : public pfasst::encap::IMEXSweeper<time>
{
    //! @{
    FFT fft;
    vector<complex<double>> ddx, lap;
    //! @}

    //! @{
    error_map errors;
    //! @}

    //! @{
    double v  = 1.0;
    time   t0 = 1.0;
    double nu = 0.02;
    size_t nf1evals = 0;
    //! @}

  public:
    //! @{
    AdvectionDiffusionSweeper(size_t nvars)
    {
      this->ddx.resize(nvars);
      this->lap.resize(nvars);
      for (size_t i = 0; i < nvars; i++) {
        double kx = 2 * PI * ((i <= nvars / 2) ? int(i) : int(i) - int(nvars));
        this->ddx[i] = complex<double>(0.0, 1.0) * kx;
        this->lap[i] = (kx * kx < 1e-13) ? 0.0 : -kx * kx;
      }
    }

    virtual ~AdvectionDiffusionSweeper()
    {
      cout << "number of f1 evals: " << this->nf1evals << endl;
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

    void echo_error(time t, bool predict = false)
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
      cout << "err: " << n << " " << k << " " << scientific << max
           << " (" << qend.size() << ", " << predict << ")"
           << endl;

      this->errors.insert(pair<pair<size_t, size_t>, double>(pair<size_t, size_t>(n, k), max));
    }

    /**
     * retrieve errors at iterations and time nodes
     */
    error_map get_errors()
    {
      return this->errors;
    }
    //! @}

    //! @{
    /**
     * @copybrief pfasst::encap::IMEXSweeper::predict()
     */
    void predict(bool initial) override
    {
      pfasst::encap::IMEXSweeper<time>::predict(initial);
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_time_step();
      this->echo_error(t + dt, true);
    }

    /**
     * @copybrief pfasst::encap::IMEXSweeper::sweep()
     */
    void sweep() override
    {
      pfasst::encap::IMEXSweeper<time>::sweep();
      time t  = this->get_controller()->get_time();
      time dt = this->get_controller()->get_time_step();
      this->echo_error(t + dt);
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

#endif
