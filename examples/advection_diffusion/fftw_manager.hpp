/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fftw_manager.hpp
 * @since v0.6.0
 */
#ifndef _EXAMPLES__ADVEC_DIFF__FFTW_MANAGER_HPP_
#define _EXAMPLES__ADVEC_DIFF__FFTW_MANAGER_HPP_

#include <map>
#include <memory>
using std::map;
using std::shared_ptr;

#include "fftw_workspace_dft1d.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      /**
       * Singleton to hold and query for FFTW workspaces.
       *
       * The @ref FFTWManager holds all instances of @ref FFTWWorkspace once queried through
       * FFTWManager::get_workspace().
       *
       * @tparam WorkspaceT type of the @ref FFTWWorkspace the manager instance is for
       */
      template<class WorkspaceT>
      class FFTWManager
      {
        public:
          //! @{
          using workspace_t = WorkspaceT;
          //! @}

        private:
          //! @{
          FFTWManager();
          FFTWManager(const FFTWManager&) = delete;
          FFTWManager(FFTWManager&&) = delete;
          FFTWManager& operator=(const FFTWManager&) = delete;
          FFTWManager& operator=(FFTWManager&&) = delete;
          //! @}

        protected:
          //! @{
          //! Storage of workspaces mapped to their number of DOFs
          map<size_t, shared_ptr<workspace_t>> _workspaces;
          //! @}

        public:
          //! @{
          virtual ~FFTWManager();
          //! @}

          //! @{
          /**
           * Get the single manager instance
           *
           * @return the one and only manager instance
           */
          static FFTWManager& get_instance();
          //! @}

          //! @{
          /**
           * Get the one @ref FFTWWorkspace for given number of DOFs
           *
           * @param[in] ndofs number of degrees of freedom of FFTW workspace
           * @return shared pointer to the only workspace with given number of DOFs
           */
          shared_ptr<WorkspaceT> get_workspace(const size_t ndofs);
          //! @}
      };

      /**
       * @class pfasst::examples::advection_diffusion::FFTWWorkspace
       * @brief Concept of a FFTW workspace
       *
       * A FFTW workspace, manageable by @ref FFTWManager, provides a wrapper around calls to FFTW
       * and persists variables and memory between calls to FFTW.
       *
       * A @ref FFTWWorkspace is expected to conform to the RAII principle.
       *
       * Usually, the initial setup of FFTW for a given number of degrees of freedom is done on
       * construction of a @ref FFTWWorkspace, e.g. allocating memory for transformed values and
       * creating _plans_.
       *
       * Cleanup and freeing of this setup should be done on destruction.
       *
       * A FFTWWorkspace has the following member functions implemented:
       *
       * * `complex<DataT::element_t>* FFTWWorkspace::forward(const DataT&)`
       *
       *   Transforms given data encapsulation to Fourier space and returns pointer to transformed
       *   values.
       *   `DataT` is the encapsulation type, e.g. pfasst::encap::VectorEncapsulation<double>,
       *   and `DataT::element_t` the type of the data at a single point in space, e.g. `double`.
       *
       * * `void FFTWWorkspace::backward(DataT&)`
       *
       *   Transformed values in Fourier space stored in `z_ptr()` back to problem space and stores
       *   them in given `DataT` object, overwriting existing data.
       *   `DataT` is the encapsulation type.
       *
       * * `size_t size() const`
       *
       *   Returns number of degrees of freedom of this workspace.
       *
       * * `complex<DataT::element_t>* z_ptr()`
       *
       *   Returns a pointer to the values in Fourier space for computations in Fourier space.
       *
       * @note This is only the specification of a concept and not available as code.
       * @ingroup Concepts
       */
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#include "fftw_manager_impl.hpp"

#endif // _EXAMPLES__ADVEC_DIFF__FFTW_MANAGER_HPP_
