/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fftw_manager_impl.hpp
 * @since v0.6.0
 */
#include "fftw_manager.hpp"

#include <utility>
using std::pair;

#include <fftw3.h>


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      template<class WorkspaceT>
      FFTWManager<WorkspaceT>::FFTWManager()
      {}

      template<class WorkspaceT>
      FFTWManager<WorkspaceT>::~FFTWManager()
      {
        fftw_cleanup();
      }

      template<class WorkspaceT>
      FFTWManager<WorkspaceT>& FFTWManager<WorkspaceT>::get_instance()
      {
        static FFTWManager<WorkspaceT> instance;
        return instance;
      }

      template<class WorkspaceT>
      shared_ptr<WorkspaceT> FFTWManager<WorkspaceT>::get_workspace(const size_t ndofs)
      {
        if (this->_workspaces.find(ndofs) == this->_workspaces.end()) {
          auto ws = make_shared<WorkspaceT>(ndofs);
          this->_workspaces.insert(pair<size_t, shared_ptr<WorkspaceT>>(ndofs, ws));
        }

        return this->_workspaces[ndofs];
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst
