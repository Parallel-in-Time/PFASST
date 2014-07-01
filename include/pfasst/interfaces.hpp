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

namespace pfasst {

  /**
   * Not implemented yet exception.
   *
   * Used by PFASST to mark methods that are required for a particular algorithm (SDC/MLSDC/PFASST)
   * that may not be necessary for all others.
   */
  class NotImplementedYet : public exception {
    string msg;
  public:
    NotImplementedYet(string msg) : msg(msg) { }
    const char *what() const throw() {
      return (string("Not implemented/supported yet, required for: ") + this->msg).c_str();
    }
  };

  /**
   * Abstract SDC sweeper.
   *
   * Note that, at this level, time is always represented as a `double`.  We do this to avoid a
   * complicated mess of templates, especially for MLSDC and PFASST.
   */
  class ISweeper {
  public:
    virtual ~ISweeper() { }

    /**
     * Setup (allocate etc) the sweeper.
     *
     * @param coarse True if this sweeper exists on a coarsened MLSDC or PFASST level.  This implies
     *     that space for an FAS correction and "saved" solutions are necessary.
     */
    virtual void setup(bool coarse=false) { }

    /**
     * Perform a predictor sweep.
     *
     * Compute a provisional solution from the initial condition.  This is typically very similar to
     * a regular SDC sweep, except that integral terms based on previous iterations don't exist yet.
     *
     * @param initial True if function values at the first node need to be computed.  False if
     *     functions values at the first node already exist (usually this is the case when advancing
     *     from one time step to the next).
     */
    virtual void predict(double t, double dt, bool initial) = 0;

    /**
     * Perform one SDC sweep/iteration.
     *
     * Compute a correction and update solution values.  Note that `sweep` can assume that valid
     * function values exist from a previous `sweep` or `predict`.
     */
    virtual void sweep(double t, double dt) = 0;

    /**
     * Advance from one time step to the next.
     *
     * Essentially this means copying the solution and function values from the last node to the
     * first node.
     */
    virtual void advance() = 0;

    /**
     * Save solutions (and/or function values) at all nodes.
     *
     * This is typically done in MLSDC/PFASST immediately after a call to restrict.  The saved
     * states are used to compute deltas during interpolation.
     */
    virtual void save() { NotImplementedYet("mlsdc/pfasst"); }

  };

  /**
   * Abstract time/space transfer (restrict/interpolate) class.
   */
  class ITransfer {
  public:
    // XXX: pass level iterator to these routines as well
    virtual ~ITransfer() { }

    /**
     * Interpolate, in time and space, from the coarse sweeper to the fine sweeper.
     *
     * @param interp_delta_from_initial True if the delta computed at each node should be relative
     *     to the initial condition.
     *
     * @param interp_initial True if a delta for the initial condtion should also be computed
     *     (PFASST).
     */
    virtual void interpolate(ISweeper *fine, const ISweeper *crse,
			     bool interp_delta_from_initial=false,
			     bool interp_initial=false) = 0;

    /**
     * Restrict, in time and space, from the fine sweeper to the coarse sweeper.
     *
     * @param restrict_initial True if the initial condition should also be restricted.
     */
    virtual void restrict(ISweeper *crse, const ISweeper *fine,
			  bool restrict_initial=false) = 0;

    /**
     * Compute FAS correction between the coarse and fine sweepers.
     */
    virtual void fas(double dt, ISweeper *crse, const ISweeper *fine) = 0;

  };

  /**
   * Abstract time communicator.
   */
  class ICommunicator {
  public:
    virtual void post() { }
    virtual void send() = 0;
    virtual void recv() = 0;
  };

}

#endif
