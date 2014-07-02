/*
 * Spectral (FFT) transfer routines.
 */

#ifndef _SPECTRAL_TRANSFER_1D_HPP_
#define _SPECTRAL_TRANSFER_1D_HPP_

#include <cassert>

#include <pfasst/encap/vector.hpp>
#include <pfasst/encap/poly_interp.hpp>

#include "fft.hpp"

template<typename time = pfasst::time_precision>
class SpectralTransfer1D : public pfasst::encap::PolyInterpMixin<time>
{
    typedef pfasst::encap::Encapsulation<double> Encapsulation;
    typedef pfasst::encap::VectorEncapsulation<double> DVectorT;

    FFT fft;

  public:
    void interpolate(Encapsulation* dst, const Encapsulation* src)
    {
      DVectorT* fine = dynamic_cast<DVectorT*>(dst);
      assert(fine != nullptr);
      const DVectorT* crse = dynamic_cast<const DVectorT*>(src);
      assert(crse != nullptr);

      this->interpolate(fine, crse);
    }

    void interpolate(DVectorT* fine, const DVectorT* crse)
    {
      auto* crse_z = fft.forward(*crse);
      auto* fine_z = fft.get_workspace(fine->size())->z;

      for (size_t i = 0; i < fine->size(); i++) {
        fine_z[i] = 0.0;
      }

      double c = 1.0 / crse->size();

      for (size_t i = 0; i < crse->size() / 2; i++) {
        fine_z[i] = c * crse_z[i];
      }

      for (size_t i = 1; i < crse->size() / 2; i++) {
        fine_z[fine->size() - crse->size() / 2 + i] = c * crse_z[crse->size() / 2 + i];
      }

      fft.backward(*fine);
    }

    void restrict(Encapsulation* dst, const Encapsulation* src)
    {
      DVectorT* crse = dynamic_cast<DVectorT*>(dst);
      assert(crse != nullptr);
      const DVectorT* fine = dynamic_cast<const DVectorT*>(src);
      assert(fine != nullptr);

      this->restrict(crse, fine);
    }

    void restrict(DVectorT* crse, const DVectorT* fine)
    {
      size_t xrat = fine->size() / crse->size();

      for (size_t i = 0; i < crse->size(); i++) {
        crse->at(i) = fine->at(xrat * i);
      }
    }

};

#endif
