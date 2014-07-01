/*
 * Vanilla SDC controller.
 */

#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include "controller.hpp"

namespace pfasst
{

  template<typename time>
  class SDC : public Controller<time>
  {
    public:

      void run()
      {
        auto* sweeper = this->get_level(0);

        for (int nstep = 0; nstep < this->nsteps; nstep++) {
          time t = nstep * this->dt;

          sweeper->predict(t, this->dt, nstep == 0);

          for (int niter = 1; niter < this->niters; niter++)
          { sweeper->sweep(t, this->dt); }

          sweeper->advance();
        }
      }

  };

}

#endif
