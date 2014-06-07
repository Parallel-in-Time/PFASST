/*
 * Vanilla SDC controller.
 */

#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include "controller.hpp"

namespace pfasst {

  template<typename time>
  class SDC : public Controller<time> {
  public:
    void run() {
      ISweeper& swp = *this->get_level(0);
      for (int nstep=0; nstep<this->nsteps; nstep++) {
	time t = nstep * this->dt;
	swp.predict(t, this->dt);
	for (int niter=0; niter<this->niters; niter++)
	  swp.sweep(t, this->dt);
	swp.advance();
      }
    }
  };

}

#endif
