/*
 * Interfaces for SDC/MLSDC/PFASST algorithms.
 */

#ifndef _PFASST_INTERFACES_HPP_
#define _PFASST_INTERFACES_HPP_

#include <deque>
#include <exception>
#include <iterator>
#include <memory>
#include <string>

using namespace std;

namespace pfasst
{

  using time_precision = double;

  /**
   * not implemented yet exception.
   *
   * Used by PFASST to mark methods that are required for a particular algorithm (SDC/MLSDC/PFASST)
   * that may not be necessary for all others.
   */
  class NotImplementedYet : public exception
  {
      string msg;
    public:
      NotImplementedYet(string msg) : msg(msg) { }
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
  class ValueError : public exception
  {
      string msg;
    public:
      ValueError(string msg) : msg(msg) { }
      const char* what() const throw()
      {
        return (string("ValueError: ") + this->msg).c_str();
      }
  };

  class ICommunicator {
  public:
    virtual ~ICommunicator() { }
    virtual int size() = 0;
    virtual int rank() = 0;
  };

  /**
   * abstract SDC sweeper.
   * @tparam time time precision
   *     defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class ISweeper
  {
    public:
      virtual ~ISweeper() { }

      /**
       * setup (allocate etc) the sweeper.
       * @param[in] coarse
       *     `true` if this sweeper exists on a coarsened MLSDC or PFASST level.
       *     This implies that space for an FAS correction and "saved" solutions are necessary.
       */
      virtual void setup(bool coarse = false) { }

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
      virtual void predict(time t, time dt, bool initial) = 0;

      /**
       * perform one SDC sweep/iteration.
       * Compute a correction and update solution values.
       * Note that this function can assume that valid function values exist from a previous
       * pfasst::ISweeper::sweep() or pfasst::ISweeper::predict().
       */
      virtual void sweep(time t, time dt) = 0;

      /**
       * advance from one time step to the next.
       *
       * Essentially this means copying the solution and function values from the last node to the
       * first node.
       */
      virtual void advance() = 0;

      /**
       * save solutions (and/or function values) at all nodes.
       *
       * This is typically done in MLSDC/PFASST immediately after a call to restrict.
       * The saved states are used to compute deltas during interpolation.
       */
      virtual void save(bool initial_only=false) { NotImplementedYet("mlsdc/pfasst"); }

    virtual void post(ICommunicator* comm) { };
    virtual void send(ICommunicator* comm) { NotImplementedYet("pfasst"); }
    virtual void recv(ICommunicator* comm) { NotImplementedYet("pfasst"); }

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
      // XXX: pass level iterator to these routines as well
      virtual ~ITransfer() { }

      /**
       * interpolate, in time and space, from the coarse sweeper to the fine sweeper.
       * @param[in] interp_delta_from_initial
       *     `true` if the delta computed at each node should be relative to the initial condition.
       * @param[in] interp_initial
       *     `true` if a delta for the initial condtion should also be computed (PFASST).
       */
      virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                               shared_ptr<const ISweeper<time>> src,
                               bool interp_delta_from_initial = false,
                               bool interp_initial = false) = 0;

      /**
       * restrict, in time and space, from the fine sweeper to the coarse sweeper.
       * @param[in] restrict_initial
       *     `true` if the initial condition should also be restricted.
       */
      virtual void restrict(shared_ptr<ISweeper<time>> dst,
                            shared_ptr<const ISweeper<time>> src,
                            bool restrict_initial = false,
			    bool restrict_initial_only = false) = 0;

      /**
       * compute FAS correction between the coarse and fine sweepers.
       */
      virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                       shared_ptr<const ISweeper<time>> src) = 0;

  };

}

#endif
