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
      FFTWWorkspaceDFT1D::FFTWWorkspaceDFT1D(const size_t ndofs)
        :   _size(ndofs)
          , _wk_ptr(fftw_alloc_complex(ndofs))
          , _z_ptr(reinterpret_cast<complex<double>*>(_wk_ptr))
      {
        this->_ffft = fftw_plan_dft_1d(ndofs, this->_wk_ptr, this->_wk_ptr, FFTW_FORWARD, FFTW_ESTIMATE);
        this->_ifft = fftw_plan_dft_1d(ndofs, this->_wk_ptr, this->_wk_ptr, FFTW_BACKWARD, FFTW_ESTIMATE);
      }

      FFTWWorkspaceDFT1D::~FFTWWorkspaceDFT1D()
      {
        fftw_free(this->_wk_ptr);
        fftw_destroy_plan(this->_ffft);
        fftw_destroy_plan(this->_ifft);
        this->_z_ptr = nullptr;
      }

      size_t FFTWWorkspaceDFT1D::size() const
      {
        return this->_size;
      }

      complex<double>* FFTWWorkspaceDFT1D::z_ptr()
      {
        return this->_z_ptr;
      }

      complex<double>* FFTWWorkspaceDFT1D::forward(const DVectorT& x)
      {
        assert(this->size() == x.size());

        for (size_t i = 0; i < this->size(); ++i) {
          this->_z_ptr[i] = x[i];
        }

        fftw_execute_dft(this->_ffft, this->_wk_ptr, this->_wk_ptr);

        return this->_z_ptr;
      }

      void FFTWWorkspaceDFT1D::backward(DVectorT & x)
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
