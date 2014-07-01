/*
 * Multi-level SDC controller.
 */

#ifndef _PFASST_MLSDC_HPP_
#define _PFASST_MLSDC_HPP_

#include <algorithm>
#include <iostream>
#include <vector>

#include "controller.hpp"

using namespace std;


namespace pfasst {

  template<typename time=double>
  class MLSDC : public Controller<time> {
    vector<size_t> nsweeps;
    bool predict, initial;

    using LevelIter = typename pfasst::Controller<time>::LevelIter;

    void perform_sweeps(LevelIter leviter, time t, time dt)
    {
      auto* sweeper = leviter.current();
      for(size_t s = 0; s < nsweeps[leviter.level]; s++)
	if(predict) {
	  sweeper->predict(t, dt, initial & predict);
	  predict = false;
	} else {
	  sweeper->sweep(t, dt);
	}
    }

  public:
    void setup()
    {
      nsweeps.resize(this->nlevels());
      fill(nsweeps.begin(), nsweeps.end(), 1);
      for(auto leviter = this->coarsest(); leviter <= this->finest(); ++leviter) {
	leviter.current()->setup(leviter!=this->finest());
      }
    }

    /**
     * Evolve ODE using MLSDC.
     *
     * This assumes that the user has set initial conditions on the
     * finest level.
     *
     * Currently uses a fixed number of iterations per step.
     */
    void run()
    {
      for(size_t nstep = 0; nstep < this->nsteps; nstep++) {
	time t = nstep * this->dt;

	predict = true;		// use predictor for first fine sweep of each step
	initial = nstep == 0;	// only evaluate node 0 functions on first step

	// iterate by performing v-cycles
	for(size_t niter = 0; niter < this->niters; niter++) {
	  cycle_v(this->finest(), t, this->dt);
        }

	// advance all levels
	for(auto leviter = this->coarsest(); leviter <= this->finest(); ++leviter) {
	  leviter.current()->advance();
        }
      }
    }

    /**
     * Cycle down: sweep on current (fine), restrict to coarse.
     */
    LevelIter cycle_down(LevelIter leviter, time t, time dt)
    {
      auto* fine = leviter.current();
      auto* crse = leviter.coarse();
      auto* trns = leviter.transfer();

      perform_sweeps(leviter, t, dt);

      trns->restrict(crse, fine, initial);
      trns->fas(dt, crse, fine);
      crse->save();

      return leviter - 1;
    }

    /**
     * Cycle up: interpolate coarse correction to fine, sweep on
     * current (fine).
     *
     * Note that if the fine level corresponds to the finest MLSDC
     * level, we don't perform a sweep.  In this case the only
     * operation that is performed here is interpolation.
     */
    LevelIter cycle_up(LevelIter leviter, time t, time dt)
    {
      auto* fine = leviter.current();
      auto* crse = leviter.coarse();
      auto* trns = leviter.transfer();

      trns->interpolate(fine, crse);

      if(leviter < this->finest()) {
	perform_sweeps(leviter, t, dt);
      }

      return leviter + 1;
    }

    /**
     * Cycle bottom: sweep on the current (coarsest) level.
     */
    LevelIter cycle_bottom(LevelIter leviter, time t, time dt)
    {
      perform_sweeps(leviter, t, dt);
      return leviter + 1;
    }

    /**
     * Perform an MLSDC V-cycle.
     */
    LevelIter cycle_v(LevelIter leviter, time t, time dt)
    {
      if(leviter.level == 0) {
      	leviter = cycle_bottom(leviter, t, dt);
      } else {
      	leviter = cycle_down(leviter, t, dt);
      	leviter = cycle_v(leviter, t, dt);
      	leviter = cycle_up(leviter, t, dt);
      }
      return leviter;
    }

  };

}

#endif
