/**
 * @file controller/pfasst.hpp
 * @since v0.1.0
 */
#ifndef _PFASST__CONTROLLER__PFASST_HPP_
#define _PFASST__CONTROLLER__PFASST_HPP_

#include "pfasst/controller/mlsdc.hpp"


namespace pfasst
{
  /**
   * Implementation of the PFASST algorithm as described in \cite emmett_pfasst_2012 .
   *
   * @ingroup Controllers
   */
  template<typename time = pfasst::time_precision>
  class PFASST
    : public MLSDC<time>
  {
      typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      ICommunicator* comm;  //!< communicator to use
      bool predict;         //!< whether to use a _predict_ sweep

      void perform_sweeps(size_t level);

    public:
      /**
       * Solve ODE using PFASST.
       *
       * @pre It is assumed that the user has set initial conditions on the finest level.
       */
      virtual void run() override;

    private:
      //! @{
      /**
       * @copydoc MLSDC::cycle_down()
       */
      LevelIter cycle_down(LevelIter level_iter);

      /**
       * @copydoc MLSDC::cycle_up()
       */
      LevelIter cycle_up(LevelIter level_iter);

      /**
       * @copydoc MLSDC::cycle_bottom()
       */
      LevelIter cycle_bottom(LevelIter level_iter);

      /**
       * @copydoc MLSDC::cycle_v()
       */
      LevelIter cycle_v(LevelIter level_iter);

      /**
       * Predictor: restrict initial down, preform coarse sweeps, return to finest.
       */
      virtual void predictor();
      //! @}

      /**
       * @name Communication
       * @{
       */
      /**
       * Broadcast finest level to all processes of PFASST::comm.
       *
       * @see ISweeper::broadcast() and its implementations for details on how broadcasting levels is
       *   done.
       */
      virtual void broadcast();

      /**
       * Generate a unique tag for level iterator.
       *
       * @param[in] level_iter level iterator providing information to compute the communication tag
       */
      virtual int tag(LevelIter level_iter);
      virtual int stag(LevelIter level_iter);

      /**
       * Post current status and values to next processor.
       */
      virtual void post();

    public:
      /**
       * Set communicator.
       *
       * @param[in] comm ICommunicator to use
       */
      virtual void set_comm(ICommunicator* comm);
      //! @}
  };
}  // ::pfasst

#include "pfasst/controller/pfasst_impl.hpp"

#endif  // _PFASST__CONTROLLER__PFASST_HPP_
