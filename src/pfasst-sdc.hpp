/*
 * Vanilla SDC controller.
 */

#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include "pfasst-controller.hpp"

namespace pfasst {

  template<typename timeT>
  class SDC : public Controller<timeT> {
  public:
    void run() {
      ISweeper& swp = *this->get_level(0);
      for (int nstep=0; nstep<this->nsteps; nstep++) {
	for (int niter=0; niter<this->niters; niter++) {
	  swp.sweep(0.0, this->dt);
	}
      }
    }
  };

}

#endif
