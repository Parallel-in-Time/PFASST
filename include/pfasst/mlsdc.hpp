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
   * Converged exception.
   *
   * Thrown when the controller detects that the finest sweeper has converged.
   */
  class ConvergedException
    : public exception
  {
  };

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
            predict = false;
          } else {
            sweeper->sweep();
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

          try {

            for (this->set_iteration(0);
                 this->get_iteration() < this->get_max_iterations();
                 this->advance_iteration()) {
              cycle_v(this->finest());
              initial = false;
            }
            perform_sweeps(this->finest().level);

          } catch (ConvergedException& e) { }

          if (this->get_time() + this->get_time_step() < this->get_end_time()) {
            this->get_finest()->advance();
          }
        }
      }

      void setup() override
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
          throw ConvergedException();
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
          l = cycle_bottom(l);
        } else {
          l = cycle_down(l);
          l = cycle_v(l);
          l = cycle_up(l);
        }
        return l;
      }

  }; // MLSDC

}  // ::pfasst

#endif
