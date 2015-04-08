/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fft.hpp
 * @since v0.1.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__FFT_HPP_
#define _EXAMPLES__ADVEC_DIFF__FFT_HPP_

#include <map>
#include <memory>
#include <cstddef>

#include <pfasst/encap/vector.hpp>
typedef pfasst::encap::VectorEncapsulation<double> DVectorT;

#include <fftw3.h>


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * FFT helper class.
       *
       * @warning Side affects galore! This is not my best work ...
       *
       * @ingroup AdvectionDiffusion
       */
      class FFT
      {
          struct workspace {
            fftw_plan        ffft;
            fftw_plan        ifft;
            fftw_complex*    wk;
            complex<double>* z;
          };

          map<size_t, shared_ptr<workspace>> workspaces;

        public:
          ~FFT()
          {
            for (auto& x : workspaces) {
              shared_ptr<workspace> wk = std::get<1>(x);
              fftw_free(wk->wk);
              fftw_destroy_plan(wk->ffft);
              fftw_destroy_plan(wk->ifft);
            }
            workspaces.clear();
          }

          shared_ptr<workspace> get_workspace(size_t ndofs)
          {
            if (workspaces.find(ndofs) == workspaces.end()) {
              shared_ptr<workspace> wk = make_shared<workspace>();
              wk->wk = fftw_alloc_complex(ndofs);
              wk->ffft = fftw_plan_dft_1d(ndofs, wk->wk, wk->wk, FFTW_FORWARD, FFTW_ESTIMATE);
              wk->ifft = fftw_plan_dft_1d(ndofs, wk->wk, wk->wk, FFTW_BACKWARD, FFTW_ESTIMATE);
              wk->z = reinterpret_cast<complex<double>*>(wk->wk);
              workspaces.insert(pair<size_t, shared_ptr<workspace>>(ndofs, wk));
            }

            return workspaces[ndofs];
          }

          complex<double>* forward(const DVectorT& x)
          {
            shared_ptr<workspace> wk = get_workspace(x.size());
            for (size_t i = 0; i < x.size(); i++) {
              wk->z[i] = x[i];
            }
            fftw_execute_dft(wk->ffft, wk->wk, wk->wk);
            return wk->z;
          }

          void backward(DVectorT& x)
          {
            shared_ptr<workspace> wk = get_workspace(x.size());
            fftw_execute_dft(wk->ifft, wk->wk, wk->wk);
            for (size_t i = 0; i < x.size(); i++) {
              x[i] = real(wk->z[i]);
            }
          }
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__ADVEC_DIFF__FFT_HPP_
