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
      FFTWManager::FFTWManager()
      {}

      FFTWManager::~FFTWManager()
      {
        fftw_cleanup();
      }

      FFTWManager& FFTWManager::get_instance()
      {
        static FFTWManager instance;
        return instance;
      }

      shared_ptr<FFTWWorkspace> FFTWManager::get_workspace(const size_t ndofs)
      {
        if (this->_workspaces.find(ndofs) == this->_workspaces.end()) {
          auto ws = make_shared<FFTWWorkspace>(ndofs);
          this->_workspaces.insert(pair<size_t, shared_ptr<FFTWWorkspace>>(ndofs, ws));
        }

        return this->_workspaces[ndofs];
      }
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst
