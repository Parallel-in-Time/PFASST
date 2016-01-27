/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fft_manager.hpp
 * @since v0.6.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__FFT_MANAGER_HPP_
#define _EXAMPLES__ADVEC_DIFF__FFT_MANAGER_HPP_

#include <map>
#include <memory>
using std::map;
using std::shared_ptr;


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Singleton to hold and query for FFTW workspaces.
       *
       * The @ref FFTManager holds all instances of @ref FFTWorkspace once queried through
       * FFTManager::get_workspace().
       *
       * On destruction a static cleanup routine on FFTWorkspace is called.
       *
       * @tparam WorkspaceT type of the @ref FFTWorkspace the manager instance is for
       */
      template<class WorkspaceT>
      class FFTManager
      {
        public:
          //! @{
          using workspace_t = WorkspaceT;
          //! @}

        private:
          //! @{
          FFTManager();
          FFTManager(const FFTManager&) = delete;
          FFTManager(FFTManager&&) = delete;
          FFTManager& operator=(const FFTManager&) = delete;
          FFTManager& operator=(FFTManager&&) = delete;
          //! @}

        protected:
          //! @{
          //! Storage of workspaces mapped to their number of DOFs
          map<size_t, shared_ptr<workspace_t>> _workspaces;
          //! @}

        public:
          //! @{
          virtual ~FFTManager();
          //! @}

          //! @{
          /**
           * Get the single manager instance
           *
           * @return the one and only manager instance
           */
          static FFTManager& get_instance();
          //! @}

          //! @{
          /**
           * Get the one @ref FFTWorkspace for given number of DOFs
           *
           * @param[in] ndofs number of degrees of freedom of FFTW workspace
           * @return shared pointer to the only workspace with given number of DOFs
           */
          shared_ptr<WorkspaceT> get_workspace(const size_t ndofs);
          //! @}
      };

      /**
       * @class pfasst::examples::advection_diffusion::FFTWorkspace
       * @brief Concept of a FFTW workspace
       *
       * A FFTW workspace, manageable by @ref FFTManager, provides a wrapper around calls to FFTW
       * and persists variables and memory between calls to FFTW.
       *
       * A @ref FFTWorkspace is expected to conform to the RAII principle.
       *
       * Usually, the initial setup of FFTW for a given number of degrees of freedom is done on
       * construction of a @ref FFTWorkspace, e.g. allocating memory for transformed values and
       * creating _plans_.
       *
       * Cleanup and freeing of this setup should be done on destruction.
       *
       * A FFTWorkspace has the following member functions implemented:
       *
       * * `complex<DataT::value_type>* FFTWorkspace::forward(const DataT&)`
       *
       *   Transforms given data encapsulation to Fourier space and returns pointer to transformed
       *   values.
       *   `DataT` is the encapsulation type, e.g. pfasst::encap::VectorEncapsulation<double>,
       *   and `DataT::value_type` the type of the data at a single point in space, e.g. `double`.
       *
       * * `void FFTWorkspace::backward(DataT&)`
       *
       *   Transformed values in Fourier space stored in `z_ptr()` back to problem space and stores
       *   them in given `DataT` object, overwriting existing data.
       *   `DataT` is the encapsulation type.
       *
       * * `size_t size() const`
       *
       *   Returns number of degrees of freedom of this workspace.
       *
       * * `complex<DataT::value_type>* z_ptr()`
       *
       *   Returns a pointer to the values in Fourier space for computations in Fourier space.
       *
       * * `static void finalize_cleanup()`
       *
       *   Routine to do static cleanup after freeing all other workspaces.
       *   Some FFT implementations require to call a cleanup function before program exit, e.g.
       *   `fftw_cleanup()`.
       *
       * @note This is only the specification of a concept and not available as code.
       * @ingroup Concepts
       */
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#include "fft_manager_impl.hpp"

#endif // _EXAMPLES__ADVEC_DIFF__FFT_MANAGER_HPP_
