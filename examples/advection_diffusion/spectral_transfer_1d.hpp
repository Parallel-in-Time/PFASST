/*
 * Spectral (FFT) transfer routines.
 */

#ifndef _SPECTRAL_TRANSFER_1D_HPP_
#define _SPECTRAL_TRANSFER_1D_HPP_

#include <pfasst/encap/vector.hpp>
#include "fft.hpp"

template<typename scalar, typename time>
class SpectralTransfer1D : public pfasst::encap::PolyInterpMixin<scalar,time> {

  using dvector = pfasst::encap::VectorEncapsulation<scalar,time>;
  FFT<scalar,time> fft;

public:

  void interpolate(Encapsulation<scalar,time> *dst, const Encapsulation<scalar,time> *src) {
    auto& crse = *dynamic_cast<const dvector*>(src);
    auto& fine = *dynamic_cast<dvector*>(dst);

    auto* crse_z = fft.forward(crse);
    auto* fine_z = fft.get_workspace(fine.size())->z;

    for (int i=0; i<fine.size(); i++)
      fine_z[i] = 0.0;

    double c = 1.0 / crse.size();

    for (int i=0; i<crse.size()/2; i++)
      fine_z[i] = c * crse_z[i];

    for (int i=1; i<crse.size()/2; i++)
      fine_z[fine.size()-crse.size()/2+i] = c * crse_z[crse.size()/2+i];

    fft.backward(fine);
  }

  void restrict(Encapsulation<scalar,time> *dst, const Encapsulation<scalar,time> *src) {
    auto& crse = *dynamic_cast<dvector*>(dst);
    auto& fine = *dynamic_cast<const dvector*>(src);

    int xrat = fine.size() / crse.size();

    for (int i=0; i<crse.size(); i++)
      crse[i] = fine[xrat*i];
  }

};

#endif
