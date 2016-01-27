/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fftw_workspace_dft1d.hpp
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
      /**
       * Workspace for DFT in 1D
       *
       * @implements FFTWWorkspace
       */
      class FFTWWorkspaceDFT1D
      {
        protected:
          //! @{
          size_t           _size;
          fftw_plan        _ffft;
          fftw_plan        _ifft;
          fftw_complex*    _wk_ptr;
          complex<double>* _z_ptr;
          //! @}

        public:
          //! @{
          /**
           * @param[in] ndofs number of DOFs of this FFTWWorkspace
           */
          explicit FFTWWorkspaceDFT1D(const size_t ndofs);
          FFTWWorkspaceDFT1D(const FFTWWorkspaceDFT1D& other) = delete;
          FFTWWorkspaceDFT1D(FFTWWorkspaceDFT1D&& other) = delete;
          virtual ~FFTWWorkspaceDFT1D();
          //! @}

          //! @{
          FFTWWorkspaceDFT1D& operator=(const FFTWWorkspaceDFT1D& other) = delete;
          FFTWWorkspaceDFT1D& operator=(FFTWWorkspaceDFT1D&& other) = delete;
          //! @}

          //! @{
          /**
           * Get number of DOFs
           *
           * @return number of degrees of freedom
           */
          size_t size() const;

          /**
           * Access values in Fourier space
           *
           * @return pointer to values in Fourier space
           */
          complex<double>* z_ptr();
          //! @}

          //! @{
          /**
           * Transforms problem data into Fourier space
           *
           * @param[in] x encapsulation holding data in problem space
           * @return pointer to values in Fourier space
           */
          complex<double>* forward(const DVectorT& x);

          /**
           * Back-transforms Fourier space data (z_ptr()) into problem space
           *
           * @param[in,out] x encapsulation to hold back-transformed data; existing data will get
           *                  overwritten
           */
          void backward(DVectorT& x);
          //! @}
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#include "fftw_workspace_dft1d_impl.hpp"

#endif  // _EXAMPLES__ADVEC_DIFF__FFTW_WORKSPACE_HPP_
