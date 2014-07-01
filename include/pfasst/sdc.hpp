/*
 * Vanilla SDC controller.
 */

#ifndef _PFASST_SDC_HPP_
#define _PFASST_SDC_HPP_

#include "controller.hpp"

namespace pfasst
{

  template<typename timeT>
  class SDC : public Controller<timeT>
  {
    public:

      void run()
      {
        auto* sweeper = this->get_level(0);

        for (int nstep = 0; nstep < this->nsteps; nstep++) {
          timeT t = nstep * this->dt;

          sweeper->predict(t, this->dt, nstep == 0);

          for (int niter = 1; niter < this->niters; niter++)
          { sweeper->sweep(t, this->dt); }

          sweeper->advance();
        }
      }

  };

}

#endif
