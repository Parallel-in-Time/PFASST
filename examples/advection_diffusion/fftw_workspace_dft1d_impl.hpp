/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fftw_workspace_dft1d_impl.hpp
 * @since v0.6.0
 */
#include "fftw_workspace_dft1d.hpp"

#include <cassert>
#include <complex>
using std::real;


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      template<class DataT>
      FFTWWorkspaceDFT1D<DataT>::FFTWWorkspaceDFT1D(const size_t ndofs)
        :   _size(ndofs)
          , _wk_ptr(fftw_alloc_complex(ndofs))
          , _z_ptr(reinterpret_cast<complex<typename DataT::value_type>*>(_wk_ptr))
      {
        this->_ffft = fftw_plan_dft_1d(ndofs, this->_wk_ptr, this->_wk_ptr,
                                       FFTW_FORWARD, FFTW_ESTIMATE);
        this->_ifft = fftw_plan_dft_1d(ndofs, this->_wk_ptr, this->_wk_ptr,
                                       FFTW_BACKWARD, FFTW_ESTIMATE);
      }

      template<class DataT>
      FFTWWorkspaceDFT1D<DataT>::~FFTWWorkspaceDFT1D()
      {
        fftw_free(this->_wk_ptr);
        fftw_destroy_plan(this->_ffft);
        fftw_destroy_plan(this->_ifft);
        this->_z_ptr = nullptr;
      }

      template<class DataT>
      void FFTWWorkspaceDFT1D<DataT>::finalize_cleanup()
      {
        fftw_cleanup();
      }

      template<class DataT>
      size_t FFTWWorkspaceDFT1D<DataT>::size() const
      {
        return this->_size;
      }

      template<class DataT>
      complex<typename DataT::value_type>* FFTWWorkspaceDFT1D<DataT>::z_ptr()
      {
        return this->_z_ptr;
      }

      template<class DataT>
      complex<typename DataT::value_type>* FFTWWorkspaceDFT1D<DataT>::forward(const DataT& x)
      {
        assert(this->size() == x.size());

        for (size_t i = 0; i < this->size(); ++i) {
          this->_z_ptr[i] = x[i];
        }

        fftw_execute_dft(this->_ffft, this->_wk_ptr, this->_wk_ptr);

        return this->_z_ptr;
      }

      template<class DataT>
      void FFTWWorkspaceDFT1D<DataT>::backward(DataT& x)
      {
        assert(this->size() == x.size());

        fftw_execute_dft(this->_ifft, this->_wk_ptr, this->_wk_ptr);

        for (size_t i = 0; i < this->size(); ++i) {
          x[i] = real(this->_z_ptr[i]);
        }
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst
