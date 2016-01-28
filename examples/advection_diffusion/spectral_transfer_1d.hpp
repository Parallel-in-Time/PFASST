/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/spectral_transfer_1d.hpp
 * @since v0.1.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__SPECTRAL_TRANSFER_1D_HPP_
#define _EXAMPLES__ADVEC_DIFF__SPECTRAL_TRANSFER_1D_HPP_

#include <cassert>
#include <cstdlib>
#include <memory>
using namespace std;

#include <pfasst/encap/vector.hpp>
#include <pfasst/encap/poly_interp.hpp>

#include "fft_manager.hpp"
#include "fftw_workspace_dft1d.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Spectral (FFT) transfer routines.
       *
       * @ingroup AdvectionDiffusion
       */
      template<typename time = pfasst::time_precision>
      class SpectralTransfer1D
        : public encap::PolyInterpMixin<time>
      {
          using Encapsulation = encap::Encapsulation<double>;

          FFTManager<FFTWWorkspaceDFT1D<encap::VectorEncapsulation<double>>> _fft;

        public:
          void interpolate(shared_ptr<Encapsulation> dst,
                           shared_ptr<const Encapsulation> src) override
          {
            auto& fine = encap::as_vector<double, time>(dst);
            auto& crse = encap::as_vector<double, time>(src);

            auto* crse_z = this->_fft.get_workspace(crse.size())->forward(crse);
            auto* fine_z = this->_fft.get_workspace(fine.size())->z_ptr();

            for (size_t i = 0; i < fine.size(); i++) {
              fine_z[i] = 0.0;
            }

            double c = 1.0 / crse.size();

            for (size_t i = 0; i < crse.size() / 2; i++) {
              fine_z[i] = c * crse_z[i];
            }

            for (size_t i = 1; i < crse.size() / 2; i++) {
              fine_z[fine.size() - crse.size() / 2 + i] = c * crse_z[crse.size() / 2 + i];
            }

            this->_fft.get_workspace(fine.size())->backward(fine);
          }

          void restrict(shared_ptr<Encapsulation> dst,
                        shared_ptr<const Encapsulation> src) override
          {
            auto& fine = encap::as_vector<double, time>(src);
            auto& crse = encap::as_vector<double, time>(dst);

            size_t xrat = fine.size() / crse.size();

            for (size_t i = 0; i < crse.size(); i++) {
              crse[i] = fine[xrat*i];
            }
          }
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#endif  // _EXAMPLES__ADVEC_DIFF__SPECTRAL_TRANSFER_1D_HPP_
