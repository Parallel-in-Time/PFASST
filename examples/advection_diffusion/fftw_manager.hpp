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

#include "fftw_workspace.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace advection_diffusion
    {
      class FFTWManager
      {
        private:
          FFTWManager();
          FFTWManager(const FFTWManager&) = delete;
          FFTWManager(FFTWManager&&) = delete;
          FFTWManager& operator=(const FFTWManager&) = delete;
          FFTWManager& operator=(FFTWManager&&) = delete;

        protected:
          map<size_t, shared_ptr<FFTWWorkspace>> _workspaces;

        public:
          virtual ~FFTWManager();

          static FFTWManager& get_instance();

          shared_ptr<FFTWWorkspace> get_workspace(const size_t ndofs);
      };
    }  // ::pfasst::examples::advection_diffusion
  }  // ::pfasst::examples
}  // ::pfasst

#include "fftw_manager_impl.hpp"

#endif // _EXAMPLES__ADVEC_DIFF__FFTW_MANAGER_HPP_
