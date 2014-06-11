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

  template<typename time>
  class MLSDC : public Controller<time> {
    vector<int> nsweeps;

  public:
    void setup() {
      nsweeps.resize(this->nlevels());
      fill(nsweeps.begin(), nsweeps.end(), 1);
      for (auto leviter=this->coarsest(); leviter<=this->finest(); ++leviter)
	leviter.current()->setup(leviter!=this->finest());
    }

    //
    // this assumes that the user has set appropriate initial
    // conditions on each level
    //
    void run() {
      for (int nstep=0; nstep<this->nsteps; nstep++) {
	time t = nstep * this->dt;

	// predictor: sweep and interpolate to finest level
	for (auto leviter=this->coarsest(); leviter<this->finest(); ++leviter) {
	  auto* sweeper = leviter.current();
	  auto* transfer = leviter.transfer();

	  if (leviter.level == 0) {
	    sweeper->predict(t, this->dt);
	    for (int s=1; s<nsweeps[leviter.level]; s++)
	      sweeper->sweep(t, this->dt);
	  } else {
	    for (int s=0; s<nsweeps[leviter.level]; s++)
	      sweeper->sweep(t, this->dt);
	  }

	  if (leviter < this->finest())
	    transfer->interpolate(leviter.fine(), leviter.current(), true);
	}

	// iterate by performing v-cycles
	for (int niter=0; niter<this->niters; niter++)
	  cycle_v(this->finest(), t, this->dt);

	// advance each level
	for (auto lev=this->coarsest(); lev<=this->finest(); ++lev)
	  lev.current()->advance();
      }
    }

    using LevelIter = typename pfasst::Controller<time>::LevelIter;

    /**
     * Cycle down: sweep on current (fine), restrict to coarse.
     */
    LevelIter cycle_down(LevelIter leviter, double t, double dt)
    {
      auto* fine = leviter.current();
      auto* crse = leviter.coarse();
      auto* trns = leviter.transfer();

      for (int s=0; s<nsweeps[leviter.level]; s++)
	fine->sweep(t, dt);

      trns->restrict(crse, fine);
      trns->fas(crse, fine);
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
    LevelIter cycle_up(LevelIter leviter, double t, double dt)
    {
      auto* fine = leviter.current();
      auto* crse = leviter.coarse();
      auto* trns = leviter.transfer();

      trns->interpolate(fine, crse, false);

      if (leviter < this->finest())
	for (int s=0; s<nsweeps[leviter.level]; s++)
	  fine->sweep(t, dt);

      return leviter + 1;
    }

    /**
     * Cycle bottom: sweep on the current (coarsest) level.
     */
    LevelIter cycle_bottom(LevelIter leviter, double t, double dt)
    {
      auto* crse = leviter.current();

      for (int s=0; s<nsweeps[leviter.level]; s++)
	crse->sweep(t, dt);

      return leviter + 1;
    }

    /**
     * Perform an MLSDC V-cycle.
     */
    LevelIter cycle_v(LevelIter leviter, double t, double dt)
    {
      if (leviter.level == 0) {
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
