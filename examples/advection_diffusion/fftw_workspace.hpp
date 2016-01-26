/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fftw_workspace.hpp
 * @since v0.6.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__FFTW_WORKSPACE_HPP_
#define _EXAMPLES__ADVEC_DIFF__FFTW_WORKSPACE_HPP_

#include <complex>
using std::complex;

#include <fftw3.h>

#include <pfasst/encap/vector.hpp>
typedef pfasst::encap::VectorEncapsulation<double> DVectorT;


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      class FFTWWorkspace
      {
        protected:
          size_t           _size;
          fftw_plan        _ffft;
          fftw_plan        _ifft;
          fftw_complex*    _wk_ptr;
          complex<double>* _z_ptr;

        public:
          explicit FFTWWorkspace(const size_t ndofs);
          FFTWWorkspace(const FFTWWorkspace& other) = delete;
          FFTWWorkspace(FFTWWorkspace&& other) = delete;
          virtual ~FFTWWorkspace();

          FFTWWorkspace& operator=(const FFTWWorkspace& other) = delete;
          FFTWWorkspace& operator=(FFTWWorkspace&& other) = delete;

          size_t size() const;
          complex<double>* z_ptr();

          complex<double>* forward(const DVectorT& x);
          void backward(DVectorT& x);
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#include "fftw_workspace_impl.hpp"

#endif  // _EXAMPLES__ADVEC_DIFF__FFTW_WORKSPACE_HPP_
