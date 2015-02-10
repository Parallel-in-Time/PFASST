#ifndef _PFASST_MLSDC_HPP_
#define _PFASST_MLSDC_HPP_

#include <vector>
using namespace std;

#include "controller.hpp"


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

      bool predict;   //<! whether to use a 'predict' sweep
      bool initial;   //<! whether we're sweeping from a new initial condition
      bool converged; //<! whether we've converged

      virtual void perform_sweeps(size_t level);

    public:
      /**
       * Solve ODE using MLSDC.
       *
       * This assumes that the user has set initial conditions on the finest level.
       * Currently uses a fixed number of iterations per step.
       */
      virtual void setup() override;
      virtual void set_nsweeps(vector<size_t> nsweeps);
      virtual void run();

    private:
      /**
       * Cycle down: sweep on current (fine), restrict to coarse.
       */
      virtual LevelIter cycle_down(LevelIter l);

      /**
       * Cycle up: interpolate coarse correction to fine, sweep on current (fine).
       *
       * Note that if the fine level corresponds to the finest MLSDC level, we don't perform a
       * sweep.
       * In this case the only operation that is performed here is interpolation.
       */
      virtual LevelIter cycle_up(LevelIter l);

      /**
       * Cycle bottom: sweep on the current (coarsest) level.
       */
      virtual LevelIter cycle_bottom(LevelIter l);

      /**
       * Perform an MLSDC V-cycle.
       */
      virtual LevelIter cycle_v(LevelIter l);
  };
}  // ::pfasst

#include "mlsdc_impl.hpp"

#endif
