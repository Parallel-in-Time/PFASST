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


namespace pfasst
{

  template<typename time = double>
  class MLSDC
    : public Controller<time>
  {
    protected:
      vector<size_t> nsweeps;
      bool predict, initial;

      typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      void perform_sweeps(LevelIter leviter)
      {
        auto sweeper = leviter.current();
        for (size_t s = 0; s < nsweeps[leviter.level]; s++) {
          if (predict) {
            sweeper->predict(initial & predict);
            predict = false;
          } else {
            sweeper->sweep();
          }
        }
      }

    public:
      void setup()
      {
        nsweeps.resize(this->nlevels());
        fill(nsweeps.begin(), nsweeps.end(), 1);
        for (auto leviter = this->coarsest(); leviter <= this->finest(); ++leviter) {
	  leviter.current()->set_controller(this);
          leviter.current()->setup(leviter != this->finest());
        }
      }

      /**
       * evolve ODE using MLSDC.
       *
       * This assumes that the user has set initial conditions on the finest level.
       * Currently uses a fixed number of iterations per step.
       */
      void run()
      {
        for (; this->get_time() < this->get_end_time(); this->advance_time()) {
          predict = true;   // use predictor for first fine sweep of each step
          initial = this->get_step() == 0; // only evaluate node 0 functions on first step

          // iterate by performing v-cycles
	  for (this->set_iteration(0); this->get_iteration() < this->get_max_iterations(); this->advance_iteration()) {
            cycle_v(this->finest());
          }

          // advance all levels
          for (auto leviter = this->coarsest(); leviter <= this->finest(); ++leviter) {
            leviter.current()->advance();
          }
        }
      }

      /**
       * cycle down: sweep on current (fine), restrict to coarse.
       */
      LevelIter cycle_down(LevelIter leviter)
      {
        auto fine = leviter.current();
        auto crse = leviter.coarse();
        auto trns = leviter.transfer();

        perform_sweeps(leviter);

        trns->restrict(crse, fine, initial);
        trns->fas(this->get_time_step(), crse, fine);
        crse->save();

        return leviter - 1;
      }

      /**
       * cycle up: interpolate coarse correction to fine, sweep on current (fine).
       *
       * Note that if the fine level corresponds to the finest MLSDC level, we don't perform a
       * sweep.
       * In this case the only operation that is performed here is interpolation.
       */
      LevelIter cycle_up(LevelIter leviter)
      {
        auto fine = leviter.current();
        auto crse = leviter.coarse();
        auto trns = leviter.transfer();

        trns->interpolate(fine, crse);

        if (leviter < this->finest()) {
          perform_sweeps(leviter);
        }

        return leviter + 1;
      }

      /**
       * cycle bottom: sweep on the current (coarsest) level.
       */
      LevelIter cycle_bottom(LevelIter leviter)
      {
        perform_sweeps(leviter);
        return leviter + 1;
      }

      /**
       * perform an MLSDC V-cycle.
       */
      LevelIter cycle_v(LevelIter leviter)
      {
        if (leviter.level == 0) {
          leviter = cycle_bottom(leviter);
        } else {
          leviter = cycle_down(leviter);
          leviter = cycle_v(leviter);
          leviter = cycle_up(leviter);
        }
        return leviter;
      }

  };

}

#endif
