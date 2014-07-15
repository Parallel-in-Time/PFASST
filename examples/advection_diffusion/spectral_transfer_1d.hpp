/*
 * Spectral (FFT) transfer routines.
 */

#ifndef _SPECTRAL_TRANSFER_1D_HPP_
#define _SPECTRAL_TRANSFER_1D_HPP_

#include <cassert>
#include <memory>

#include <pfasst/encap/vector.hpp>
#include <pfasst/encap/poly_interp.hpp>

#include "fft.hpp"

template<typename time = pfasst::time_precision>
class SpectralTransfer1D
  : public pfasst::encap::PolyInterpMixin<time>
{
    typedef pfasst::encap::Encapsulation<double> Encapsulation;

    FFT fft;

  public:

    void interpolate(shared_ptr<Encapsulation> dst, shared_ptr<const Encapsulation> src)
    {
      auto& fine = as_vector<double,time>(dst);
      auto& crse = as_vector<double,time>(src);

      auto* crse_z = fft.forward(crse);
      auto* fine_z = fft.get_workspace(fine.size())->z;

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

      fft.backward(fine);
    }

    void restrict(shared_ptr<Encapsulation> dst, shared_ptr<const Encapsulation> src)
    {
      auto& fine = as_vector<double,time>(src);
      auto& crse = as_vector<double,time>(dst);

      size_t xrat = fine.size() / crse.size();

      for (size_t i = 0; i < crse.size(); i++) {
        crse[i] = fine[xrat*i];
      }
    }

};

#endif
