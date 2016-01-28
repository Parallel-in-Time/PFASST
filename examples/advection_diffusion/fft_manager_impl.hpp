/**
 * @ingroup AdvectionDiffusionFiles
 * @file examples/advection_diffusion/fft_manager_impl.hpp
 * @since v0.6.0
 */
#include "fft_manager.hpp"

#include <utility>
using std::pair;


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      template<class WorkspaceT>
      shared_ptr<WorkspaceT> FFTManager<WorkspaceT>::get_workspace(const size_t ndofs)
      {
        if (this->_workspaces.find(ndofs) == this->_workspaces.end()) {
          auto ws = make_shared<WorkspaceT>(ndofs);
          this->_workspaces.insert(pair<size_t, shared_ptr<WorkspaceT>>(ndofs, ws));
        }

        return this->_workspaces[ndofs];
      }

      template<class WorkspaceT>
      void FFTManager<WorkspaceT>::finalize_cleanup()
      {
        WorkspaceT::finalize_cleanup();
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst
