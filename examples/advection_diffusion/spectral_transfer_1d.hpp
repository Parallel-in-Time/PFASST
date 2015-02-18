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

#include "fft.hpp"


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
          typedef encap::Encapsulation<double> Encapsulation;

          FFT fft;

        public:
          void interpolate(shared_ptr<Encapsulation> dst, shared_ptr<const Encapsulation> src) override
          {
            auto& fine = encap::as_vector<double, time>(dst);
            auto& crse = encap::as_vector<double, time>(src);

            auto* crse_z = this->fft.forward(crse);
            auto* fine_z = this->fft.get_workspace(fine.size())->z;

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

            this->fft.backward(fine);
          }

          void restrict(shared_ptr<Encapsulation> dst, shared_ptr<const Encapsulation> src) override
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
