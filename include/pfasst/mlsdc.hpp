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

  /**
   * Multilevel SDC controller.
   */
  template<typename time = pfasst::time_precision>
  class MLSDC
    : public Controller<time>
  {
    protected:
      vector<size_t> nsweeps;

      typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      bool predict; //<! whether to use a 'predict' sweep
      bool initial; //<! whether we're sweeping from a new initial condition

      void perform_sweeps(size_t level)
      {
        auto sweeper = this->get_level(level);
        for (size_t s = 0; s < this->nsweeps[level]; s++) {
          if (predict) {
            sweeper->predict(initial & predict);
            sweeper->post_predict();
            predict = false;
          } else {
            sweeper->sweep();
            sweeper->post_sweep();
          }
        }
      }

    public:
      /**
       * Solve ODE using MLSDC.
       *
       * This assumes that the user has set initial conditions on the finest level.
       * Currently uses a fixed number of iterations per step.
       */
      void run()
      {
        for (; this->get_time() < this->get_end_time(); this->advance_time()) {
          predict = true;
          initial = true;

          for (this->set_iteration(0);
               this->get_iteration() < this->get_max_iterations();
               this->advance_iteration()) {
            cycle_v(this->finest());
            initial = false;
          }

          perform_sweeps(this->finest().level);

          this->get_finest()->post_step();

          if (this->get_time() + this->get_time_step() < this->get_end_time()) {
            this->get_finest()->advance();
          }
        }
      }

      virtual void setup() override
      {
        nsweeps.resize(this->nlevels());
        fill(nsweeps.begin(), nsweeps.end(), 1);
        for (auto leviter = this->coarsest(); leviter <= this->finest(); ++leviter) {
          leviter.current()->set_controller(this);
          leviter.current()->setup(leviter != this->finest());
        }
      }

      void set_nsweeps(vector<size_t> nsweeps)
      {
        this->nsweeps = nsweeps;
      }

    private:
      /**
       * Cycle down: sweep on current (fine), restrict to coarse.
       */
      LevelIter cycle_down(LevelIter l)
      {
        auto fine = l.current();
        auto crse = l.coarse();
        auto trns = l.transfer();

        perform_sweeps(l.level);

        if (l == this->finest() && fine->converged()) {
          l.set_converged();
          return l;
        }

        trns->restrict(crse, fine, initial);
        trns->fas(this->get_time_step(), crse, fine);
        crse->save();

        return l - 1;
      }

      /**
       * Cycle up: interpolate coarse correction to fine, sweep on current (fine).
       *
       * Note that if the fine level corresponds to the finest MLSDC level, we don't perform a
       * sweep.
       * In this case the only operation that is performed here is interpolation.
       */
      LevelIter cycle_up(LevelIter l)
      {
        auto fine = l.current();
        auto crse = l.coarse();
        auto trns = l.transfer();

        trns->interpolate(fine, crse);

        if (l < this->finest()) {
          perform_sweeps(l.level);
        }

        return l + 1;
      }

      /**
       * Cycle bottom: sweep on the current (coarsest) level.
       */
      LevelIter cycle_bottom(LevelIter l)
      {
        perform_sweeps(l.level);
        return l + 1;
      }

      /**
       * Perform an MLSDC V-cycle.
       */
      LevelIter cycle_v(LevelIter l)
      {
        if (l.level == 0) {
          l = bind_if_not_converged(l, cycle_bottom);
        } else {
          l = bind_if_not_converged(l, cycle_down);
          l = bind_if_not_converged(l, cycle_v);
          l = bind_if_not_converged(l, cycle_up);
        }
        return l;
      }
  };

}  // ::pfasst

#endif
