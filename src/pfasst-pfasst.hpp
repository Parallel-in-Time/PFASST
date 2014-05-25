/*
 * PFASST controller.
 */

#ifndef _PFASST_PFASST_HPP_
#define _PFASST_PFASST_HPP_

#include "pfasst-controller.hpp"

using namespace std;

namespace pfasst {

  template<typename timeT>
  class PFASST : public Controller<timeT> {

    void run(timeT dt, int nsteps) {

    }

    // leveliter cycle_down(leveliter levels, double t, double dt);
    // leveliter cycle_up(leveliter levels, double t, double dt);
    // leveliter cycle_bottom(leveliter levels, double t, double dt);
    // leveliter cycle_top(leveliter levels, double t, double dt);
    // leveliter cycle_v(leveliter levels, double t, double dt);
  };

}

#endif
