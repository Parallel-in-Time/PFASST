/*
 * Interfaces for SDC/MLSDC/PFASST algorithms.
 */

#ifndef _PFASST_INTERFACES_HPP_
#define _PFASST_INTERFACES_HPP_

#include <cassert>
#include <deque>
#include <exception>
#include <iterator>
#include <memory>
#include <string>

#include "globals.hpp"

using namespace std;

namespace pfasst
{

  using time_precision = double;

  // forward declare for ISweeper
  template<typename time>
  class Controller;

  /**
   * not implemented yet exception.
   *
   * Used by PFASST to mark methods that are required for a particular algorithm (SDC/MLSDC/PFASST)
   * that may not be necessary for all others.
   */
  class NotImplementedYet
    : public exception
  {
      string msg;
    public:
      NotImplementedYet(string msg)
        : msg(msg)
      {}

      const char* what() const throw()
      {
        return (string("Not implemented/supported yet, required for: ") + this->msg).c_str();
      }
  };


  /**
   * value exception.
   *
   * Thrown when a PFASST routine is passed an invalid value.
   */
  class ValueError
    : public exception
  {
      string msg;
    public:
      ValueError(string msg)
        : msg(msg)
      {}

      const char* what() const throw()
      {
        return (string("ValueError: ") + this->msg).c_str();
      }
  };

  class ICommunicator
  {
    public:
      virtual ~ICommunicator() { }
      virtual int size() = 0;
      virtual int rank() = 0;
      virtual void set_converged(bool converged) = 0;
      virtual bool get_converged(int rank) = 0;
      virtual void clear_converged() = 0;
      virtual void fence_status() = 0;
  };

  /**
   * abstract SDC sweeper.
   * @tparam time time precision
   *     defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class ISweeper
  {
    protected:
      Controller<time>* controller;

    public:
      //! @{
      ISweeper()
        : controller(nullptr)
      {}

      virtual ~ISweeper()
      {}
      //! @}

      //! @{
      /**
       * set the sweepers controller.
       */
      void set_controller(Controller<time>* ctrl)
      {
        this->controller = ctrl;
      }

      Controller<time>* get_controller()
      {
        assert(this->controller);
        return this->controller;
      }
      //! @}

      //! @{
      /**
       * setup (allocate etc) the sweeper.
       * @param[in] coarse
       *     `true` if this sweeper exists on a coarsened MLSDC or PFASST level.
       *     This implies that space for an FAS correction and "saved" solutions are necessary.
       */
      virtual void setup(bool coarse = false)
      {
        UNUSED(coarse);
      }

      /**
       * perform a predictor sweep.
       *
       * Compute a provisional solution from the initial condition.
       * This is typically very similar to a regular SDC sweep, except that integral terms based on
       * previous iterations don't exist yet.
       * @param[in] initial
       *     `true` if function values at the first node need to be computed.
       *     `false` if functions values at the first node already exist (usually this is the case
       *     when advancing from one time step to the next).
       */
      virtual void predict(bool initial) = 0;

      /**
       * Perform one SDC sweep/iteration.
       *
       * Compute a correction and update solution values.  Note that this function can assume that
       * valid function values exist from a previous pfasst::ISweeper::sweep() or
       * pfasst::ISweeper::predict().
       */
      virtual void sweep() = 0;

      /**
       * Advance from one time step to the next.
       *
       * Essentially this means copying the solution and function values from the last node to the
       * first node.
       */
      virtual void advance() = 0;

      /**
       * Return convergence status.
       *
       * This is used by controllers to shortcircuit iterations.
       */
      virtual bool converged() { return false; }

      /**
       * Save states (and/or function values) at all nodes.
       *
       * This is typically done in MLSDC/PFASST immediately after a call to restrict.
       * The saved states are used to compute deltas during interpolation.
       *
       * @note This method must be implemented in derived sweepers.
       */
      virtual void save(bool initial_only=false)
      {
        UNUSED(initial_only);
        throw NotImplementedYet("mlsdc/pfasst");
      }

      virtual void spread()
      {
        throw NotImplementedYet("pfasst");
      }
      //! @}

      //! @{
      virtual void post_sweep() { }
      virtual void post_predict() { }
      virtual void post_step() { }
      //! @}

      //! @{
      virtual void post(ICommunicator* comm, int tag)
      {
        UNUSED(comm); UNUSED(tag);
      };

      virtual void send(ICommunicator* comm, int tag, bool blocking)
      {
        UNUSED(comm); UNUSED(tag); UNUSED(blocking);
        throw NotImplementedYet("pfasst");
      }

      virtual void recv(ICommunicator* comm, int tag, bool blocking)
      {
        UNUSED(comm); UNUSED(tag); UNUSED(blocking);
        throw NotImplementedYet("pfasst");
      }

      virtual void broadcast(ICommunicator* comm)
      {
        UNUSED(comm);
        throw NotImplementedYet("pfasst");
      }
      //! @}

  };

  /**
   * abstract time/space transfer (restrict/interpolate) class.
   * @tparam time time precision
   *     defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class ITransfer
  {
    public:
      //! @{
      virtual ~ITransfer()
      {}
      //! @}

      //! @{
      /**
       * Interpolate initial condition from the coarse sweeper to the fine sweeper.
       */
      virtual void interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                       shared_ptr<const ISweeper<time>> src)
      {
        UNUSED(dst); UNUSED(src);
        throw NotImplementedYet("pfasst");
      }

      /**
       * Interpolate, in time and space, from the coarse sweeper to the fine sweeper.
       * @param[in] interp_initial
       *     `true` if a delta for the initial condtion should also be computed (PFASST).
       */
      virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                               shared_ptr<const ISweeper<time>> src,
                               bool interp_initial = false) = 0;

      /**
       * Restrict initial condition from the fine sweeper to the coarse sweeper.
       * @param[in] restrict_initial
       *     `true` if the initial condition should also be restricted.
       */
      virtual void restrict_initial(shared_ptr<ISweeper<time>> dst,
                                    shared_ptr<const ISweeper<time>> src)
      {
        UNUSED(dst); UNUSED(src);
        throw NotImplementedYet("pfasst");
      }


      /**
       * Restrict, in time and space, from the fine sweeper to the coarse sweeper.
       * @param[in] restrict_initial
       *     `true` if the initial condition should also be restricted.
       */
      virtual void restrict(shared_ptr<ISweeper<time>> dst,
                            shared_ptr<const ISweeper<time>> src,
                            bool restrict_initial = false) = 0;


      /**
       * Compute FAS correction between the coarse and fine sweepers.
       */
      virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                       shared_ptr<const ISweeper<time>> src) = 0;
      //! @}
  };

}  // ::pfasst

#endif
