#ifndef _PFASST_PFASST_HPP_
#define _PFASST_PFASST_HPP_

#include "pfasst/controller/mlsdc.hpp"


namespace pfasst
{
  /**
   * implementation of the PFASST algorithm as described in \cite emmett_pfasst_2012
   */
  template<typename time = pfasst::time_precision>
  class PFASST
    : public MLSDC<time>
  {
      ICommunicator* comm;

      typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      bool predict; //<! whether to use a 'predict' sweep

      virtual void perform_sweeps(size_t level);

    public:
      /**
       * Evolve ODE using PFASST.
       *
       * This assumes that the user has set initial conditions on the
       * finest level.
       *
       * Currently uses "block mode" PFASST with the standard predictor.
       */
      virtual void run();

      virtual void set_comm(ICommunicator* comm);

    private:
      /**
       * Cycle down: sweep on current (fine), restrict to coarse.
       */
      virtual LevelIter cycle_down(LevelIter l);

      /**
       * Cycle up: interpolate coarse correction to fine, sweep on
       * current (fine).
       *
       * Note that if the fine level corresponds to the finest MLSDC
       * level, we don't perform a sweep.  In this case the only
       * operation that is performed here is interpolation.
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

      /**
       * Predictor: restrict initial down, preform coarse sweeps, return to finest.
       */
      virtual void predictor();

      virtual void broadcast();

      virtual int tag(LevelIter l);

      virtual void post();
  };
}  // ::pfasst

#include "pfasst/controller/pfasst_impl.hpp"

#endif
