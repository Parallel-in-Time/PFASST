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
	timeT t = nstep * this->dt;
	swp.predict(t, this->dt);
	for (int niter=0; niter<this->niters; niter++) {
	  swp.sweep(t, this->dt);
	}
	swp.advance();
      }
    }
  };

}

#endif
