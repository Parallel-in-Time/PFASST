/**
 * @file controller.mlsdc.hpp
 * @since v0.1.0
 */
#ifndef _PFASST__CONTROLLER__MLSDC_HPP_
#define _PFASST__CONTROLLER__MLSDC_HPP_

#include <vector>
using namespace std;

#include "pfasst/controller/interface.hpp"


namespace pfasst
{
  /**
   * multilevel SDC controller.
   *
   * @tparam time time precision
   *
   * @ingroup Controllers
   */
  template<typename time = pfasst::time_precision>
  class MLSDC
    : public Controller<time>
  {
    protected:
      vector<size_t> nsweeps;  //!< how many sweeps should be done on the different levels

      typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      bool predict;   //!< whether to use a _predict_ sweep
      bool initial;   //!< whether we're sweeping from a new initial condition
      bool converged; //!< whether we've converged

      /**
       * perform pre-configured number of sweeps on level @p level.
       *
       * @param[in] level level index to perform pre-configured number of sweeps
       * @see set_nsweeps() for pre-configuring number of sweeps
       */
      virtual void perform_sweeps(size_t level);

    public:
      //! @copydoc Controller::setup()
      virtual void setup() override;

      /**
       * set desired number of sweeps for each level independently.
       *
       * @param[in] nsweeps vector with number of sweeps for each level.
       *
       * @note The same level indicing as with Controller::levels applies.
       */
      virtual void set_nsweeps(vector<size_t> nsweeps);

      /**
       * solve ODE using MLSDC.
       *
       * @pre It is assumed that the user has set initial conditions on the finest level.
       */
      virtual void run();

    private:
      /**
       * sweep on current (fine), restrict to coarse.
       *
       * @param[in] level_iter level iterator pointing to the finer level
       * @returns level iterator to current level if it has converged after the sweep; otherwise a
       *   level iterator to the next coarser level
       */
      virtual LevelIter cycle_down(LevelIter level_iter);

      /**
       * interpolate coarse correction to fine, sweep on current (fine).
       *
       * @param[in] level_iter level iterator pointing to the finer level
       * @returns level iterator pointing to the coarser level
       */
      virtual LevelIter cycle_up(LevelIter level_iter);

      /**
       * sweep on the current (coarsest) level.
       *
       * @param[in] level_iter level iterator pointing to the coarsest level
       * @returns level iterator pointing to the next finer level
       *
       * @pre It is assumed that @p level_iter currently points to the coarsest level.
       */
      virtual LevelIter cycle_bottom(LevelIter level_iter);

      /**
       * perform an MLSDC V-cycle.
       *
       * @param[in] level_iter level iterator pointing to a fine level
       * @returns level iterator pointing to either the coarsest level if @p level_iter is the 
       *   coarsest level or the same level as @p level_iter if a sweep on that level results in
       *   convergence or all sweeps on all coarser levels did not let to convergence
       */
      virtual LevelIter cycle_v(LevelIter level_iter);
  };
}  // ::pfasst

#include "pfasst/controller/mlsdc_impl.hpp"

#endif  // _PFASST__CONTROLLER__MLSDC_HPP_
