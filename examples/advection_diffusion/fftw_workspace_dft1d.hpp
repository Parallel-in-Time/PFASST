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


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Workspace for DFT in 1D
       *
       * @tparam DataT type of the encapsulation in the problem space; must provide public member
       *               `value_type` providing the type of single data points, public member
       *               function `size()` returning number of elements and public member subscript
       *               operator.
       * @implements FFTWWorkspace
       */
      template<class DataT>
      class FFTWWorkspaceDFT1D
      {
        public:
          using data_type = DataT;

        protected:
          //! @{
          size_t                                   _size;
          fftw_plan                                _ffft;
          fftw_plan                                _ifft;
          fftw_complex*                            _wk_ptr;
          complex<typename data_type::value_type>* _z_ptr;
          //! @}

        public:
          //! @{
          /**
           * @param[in] ndofs number of DOFs of this FFTWWorkspace
           */
          explicit FFTWWorkspaceDFT1D(const size_t ndofs);
          FFTWWorkspaceDFT1D(const FFTWWorkspaceDFT1D<DataT>& other) = delete;
          FFTWWorkspaceDFT1D(FFTWWorkspaceDFT1D<DataT>&& other) = delete;
          virtual ~FFTWWorkspaceDFT1D();
          //! @}

          //! @{
          FFTWWorkspaceDFT1D& operator=(const FFTWWorkspaceDFT1D<DataT>& other) = delete;
          FFTWWorkspaceDFT1D& operator=(FFTWWorkspaceDFT1D<DataT>&& other) = delete;
          //! @}

          //! @{
          /**
           * Calls final cleanup routine of FFTW.
           */
          static void finalize_cleanup();
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
          complex<typename DataT::value_type>* z_ptr();
          //! @}

          //! @{
          /**
           * Transforms problem data into Fourier space
           *
           * @param[in] x encapsulation holding data in problem space
           * @return pointer to values in Fourier space
           */
          complex<typename DataT::value_type>* forward(const DataT& x);

          /**
           * Back-transforms Fourier space data (z_ptr()) into problem space
           *
           * @param[in,out] x encapsulation to hold back-transformed data; existing data will get
           *                  overwritten
           */
          void backward(DataT& x);
          //! @}
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#include "fftw_workspace_dft1d_impl.hpp"

#endif  // _EXAMPLES__ADVEC_DIFF__FFTW_WORKSPACE_HPP_
