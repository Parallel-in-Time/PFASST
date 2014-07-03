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
   * @brief not implemented yet exception
   * @details Used by PFASST to mark methods that are required for a particular algorithm 
   *     (SDC/MLSDC/PFASST) that may not be necessary for all others.
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
   * @brief value exception
   * @details Thrown when a PFASST routine is passed an invalid value.
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

  /**
   * @brief abstract SDC sweeper
   * @tparam time time precision
   *     defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class ISweeper
  {
    public:
      virtual ~ISweeper() { }

      /**
       * @brief setup (allocate etc) the sweeper
       * @param[in] coarse
       *     `true` if this sweeper exists on a coarsened MLSDC or PFASST level.
       *     This implies that space for an FAS correction and "saved" solutions are necessary.
       */
      virtual void setup(bool coarse = false) { }

      /**
       * @brief perform a predictor sweep
       * @details Compute a provisional solution from the initial condition.
       *     This is typically very similar to a regular SDC sweep, except that integral terms based
       *     on previous iterations don't exist yet.
       * @param[in] initial
       *     `true` if function values at the first node need to be computed.
       *     `false` if functions values at the first node already exist (usually this is the case 
       *     when advancing from one time step to the next).
       */
      virtual void predict(time t, time dt, bool initial) = 0;

      /**
       * @brief perform one SDC sweep/iteration
       * @details Compute a correction and update solution values.
       *     Note that this function can assume that valid function values exist from a previous 
       *     pfasst::ISweeper::sweep() or pfasst::ISweeper::predict().
       */
      virtual void sweep(time t, time dt) = 0;

      /**
       * @brief advance from one time step to the next
       * @details Essentially this means copying the solution and function values from the last 
       *     node to the first node.
       */
      virtual void advance() = 0;

      /**
       * @brief save solutions (and/or function values) at all nodes
       * @details This is typically done in MLSDC/PFASST immediately after a call to restrict.
       *     The saved states are used to compute deltas during interpolation.
       */
      virtual void save() { NotImplementedYet("mlsdc/pfasst"); }

  };

  /**
   * @brief abstract time/space transfer (restrict/interpolate) class
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
       * @brief interpolate, in time and space, from the coarse sweeper to the fine sweeper
       * @param[in] interp_delta_from_initial
       *     `true` if the delta computed at each node should be relative to the initial condition.
       * @param[in] interp_initial
       *     `true` if a delta for the initial condtion should also be computed (PFASST).
       */
      virtual void interpolate(ISweeper<time>* dst, const ISweeper<time>* src,
                               bool interp_delta_from_initial = false,
                               bool interp_initial = false) = 0;

      /**
       * @brief restrict, in time and space, from the fine sweeper to the coarse sweeper
       * @param[in] restrict_initial
       *     `true` if the initial condition should also be restricted.
       */
      virtual void restrict(ISweeper<time>* dst, const ISweeper<time>* src,
                            bool restrict_initial = false) = 0;

      /**
       * @brief compute FAS correction between the coarse and fine sweepers
       */
      virtual void fas(time dt, ISweeper<time>* dst, const ISweeper<time>* src) = 0;

  };

  /**
   * @brief abstract time communicator
   */
  class ICommunicator
  {
    public:
      virtual void post() { }
      virtual void send() = 0;
      virtual void recv() = 0;
  };

}

#endif
