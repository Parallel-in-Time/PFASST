#include "pfasst/transfer/spectral_1d.hpp"

#include <cassert>
#include <memory>
#include <vector>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"
#include "pfasst/quadrature.hpp"
#include "pfasst/encap/vector.hpp"


namespace pfasst
{
  template<class TransferTraits, typename Enabled>
  void
  Spectral1DTransfer<TransferTraits, Enabled>::interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                                                shared_ptr<typename TransferTraits::fine_encap_type> fine)
  {
    CVLOG(1, "TRANS") << "interpolate data";

    const size_t coarse_ndofs = coarse->get_data().size();
    const size_t fine_ndofs = fine->get_data().size();
    assert(coarse_ndofs > 0);
    assert(fine_ndofs >= coarse_ndofs);
    
    if (fine_ndofs == coarse_ndofs) {
      // do a shortcut
      CLOG(DEBUG, "TRANS") << "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT";
      fine->data() = coarse->get_data();

    } else {
      complex<fine_spacial_type> *coarse_z = this->fft.forward(coarse);
      complex<fine_spacial_type> *fine_z = this->fft.get_workspace(fine_ndofs)->z;

      for (size_t i = 0; i < fine_ndofs; i++) {
        fine_z[i] = 0.0;
      }

      // FFTW is not normalized
      coarse_spacial_type c = 1.0 / coarse_ndofs;

      // positive frequencies
      for (size_t i = 0; i < coarse_ndofs / 2; i++) {
        fine_z[i] = c * coarse_z[i];
      }

      // negative frequencies (in backward order)
      for (size_t i = 1; i < coarse_ndofs / 2; i++) {
        fine_z[fine_ndofs - coarse_ndofs / 2 + i] = c * coarse_z[coarse_ndofs / 2 + i];
      }

      this->fft.backward(fine);
    }
  }

  template<class TransferTraits, typename Enabled>
  void
  Spectral1DTransfer<TransferTraits, Enabled>::restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                                             shared_ptr<typename TransferTraits::coarse_encap_type> coarse)
  {
    CVLOG(1, "TRANS") << "restrict data";

    const size_t coarse_ndofs = coarse->get_data().size();
    const size_t fine_ndofs = fine->get_data().size();
    assert(coarse_ndofs > 0);
    assert(fine_ndofs >= coarse_ndofs);
    
    if (fine_ndofs == coarse_ndofs) {
      // do a shortcut
      CLOG(DEBUG, "TRANS") << "number dofs of fine and coarse are the same; doing a trivial copy and NO FFT";
      coarse->data() = fine->get_data();

    } else {
      const size_t factor = fine_ndofs / coarse_ndofs;

      for (size_t i = 0; i < coarse_ndofs; i++) {
        coarse->data()[i] = fine->get_data()[factor * i];
      }
    }
  }
}  // ::pfasst
