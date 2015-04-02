/**
 * interfaces for SDC/MLSDC/PFASST algorithms.
 *
 * @file pfasst/interfaces.hpp
 * @since v0.1.0
 */
#ifndef _PFASST_INTERFACES_HPP_
#define _PFASST_INTERFACES_HPP_

#include <memory>
#include <stdexcept>
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
   *
   * @since v0.1.0
   */
  class NotImplementedYet
    : public runtime_error
  {
    protected:
      string msg;

    public:
      /**
       * @param[in] msg component or algorithm the throwing function is required for
       */
      explicit NotImplementedYet(const string& msg);
      virtual const char* what() const throw();
  };


  /**
   * value exception.
   *
   * Thrown when a PFASST routine is passed an invalid value.
   *
   * @since v0.1.0
   *
   * @todo Consider deprecating this in favour of std::invalid_argument.
   */
  class ValueError
    : public invalid_argument
  {
    protected:
      string msg;

    public:
      explicit ValueError(const string& msg);
      virtual const char* what() const throw();
  };


  // forward declare for IStatus
  class IStatus;


  /**
   * abstract interface for communicators.
   *
   * The interface ensures a communicator provides the notion of the total number of processors
   * (i.e. `size()`) and and the ID of the _current_ processor (i.e. `rank()`) as well as the
   * current @ref ICommunicator::status "status" of the algorithm.
   */
  class ICommunicator
  {
    public:
      virtual ~ICommunicator();
      virtual int size() = 0;
      virtual int rank() = 0;

      shared_ptr<IStatus> status;
  };


  /**
   * abstract interface for the current status of the algorithm.
   *
   * The status requires a @ref ICommunicator "communicator" to enable sending and receiving stati
   * of other processors.
   */
  class IStatus
  {
    protected:
      ICommunicator* comm;

    public:
      virtual ~IStatus();

      //! @{
      /**
       * resetting status.
       *
       * @note Logic is implementation defined.
       */
      virtual void clear() = 0;

      /**
       * sets a new converged state.
       */
      virtual void set_converged(bool converged) = 0;

      /**
       * retreive converged state for specific processor.
       *
       * @param[in] rank ID of processor to check converged state for
       * @returns `true` if processor with ID @p rank has converged; `false` otherwise
       *
       * @note Inspection logic is implementation defined.
       */
      virtual bool get_converged(int rank) = 0;
      //! @}

      //! @{
      /**
       * set new communicator to use.
       */
      virtual void set_comm(ICommunicator* comm);

      /**
       * check whether previous processor is still iterating.
       *
       * @returns `true` if previous processor has converged; `false` if it is still iterating
       */
      virtual bool previous_is_iterating();

      /**
       * check whether this processor should keep iterating.
       *
       * @returns `true` if this processor should keep iterating; `false` if it should switch to
       *   `converged` state
       */
      virtual bool keep_iterating();
      //! @}

      //! @{
      virtual void post() = 0;
      virtual void send() = 0;
      virtual void recv() = 0;
      //! @}
  };


  // forward declare for ISweeper
  template<typename time>
  class Controller;


  /**
   * abstract SDC sweeper.
   *
   * @tparam time time precision; defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class ISweeper
  {
    protected:
      /**
       * backreference to the controller managing the sweeper instance.
       */
      Controller<time>* controller;

    public:
      //! @{
      ISweeper();
      virtual ~ISweeper();
      //! @}

      //! @{
      /**
       * set the sweepers controller.
       * 
       * @param[in] ctrl new controller to manage this sweeper
       */
      virtual void set_controller(Controller<time>* ctrl);

      /**
       * accessor to the controller managing this sweeper.
       *
       * @returns controller managing this sweeper
       */
      virtual Controller<time>* get_controller();
      //! @}

      //! @{
      /**
       * set options from command line etc.
       */
      virtual void set_options();

      /**
       * Setup (allocate etc) the sweeper.
       *
       * @param[in] coarse `true` if this sweeper exists on a coarsened MLSDC or PFASST level.
       *   This implies that space for an FAS correction and "saved" solutions are necessary.
       */
      virtual void setup(bool coarse = false);

      /**
       * perform a predictor sweep.
       *
       * Compute a provisional solution from the initial condition.
       * This is typically very similar to a regular SDC sweep, except that integral terms based on
       * previous iterations don't exist yet.
       *
       * @param[in] initial `true` if function values at the first node need to be computed.
       *   `false` if functions values at the first node already exist (usually this is the case
       *   when advancing from one time step to the next).
       */
      virtual void predict(bool initial) = 0;

      /**
       * perform one SDC sweep/iteration.
       *
       * Compute a correction and update solution values.
       * Note that this function can assume that valid function values exist from a previous 
       * pfasst::ISweeper::sweep() or pfasst::ISweeper::predict().
       */
      virtual void sweep() = 0;

      /**
       * advance from one time step to the next.
       *
       * Essentially this means copying the solution and function values from the last node to the
       * first node.
       */
      virtual void advance() = 0;

      /**
       * return convergence status.
       *
       * This is used by controllers to shortcircuit iterations.
       */
      virtual bool converged();

      /**
       * save states (and/or function values) at all nodes.
       *
       * This is typically done in MLSDC/PFASST immediately after a call to restrict.
       * The saved states are used to compute deltas during interpolation.
       *
       * @param[in] initial_only flag indicating whether only the initial state should be saved
       *
       * @note This method must be implemented in derived sweepers.
       */
      virtual void save(bool initial_only=false);

      /**
       * initialize solution values at all time nodes with meaningful values.
       */
      virtual void spread();
      //! @}

      //! @{
      /**
       * hook automatically run after each completed sweep.
       */
      virtual void post_sweep();

      /**
       * hook automatically run after each completed predict.
       */
      virtual void post_predict();

      /**
       * hook automatically run after each completed time step.
       */
      virtual void post_step();
      //! @}

      //! @{
      virtual void post(ICommunicator* comm, int tag);
      virtual void send(ICommunicator* comm, int tag, bool blocking);
      virtual void recv(ICommunicator* comm, int tag, bool blocking);
      virtual void broadcast(ICommunicator* comm);
      //! @}
  };


  /**
   * abstract time/space transfer (restrict/interpolate) class.
   *
   * @tparam time time precision; defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class ITransfer
  {
    public:
      //! @{
      virtual ~ITransfer();
      //! @}

      //! @{
      /**
       * interpolate initial condition from the coarse sweeper to the fine sweeper.
       *
       * @param[in,out] dst sweeper to interpolate onto (i.e. fine level)
       * @param[in]     src sweeper to interpolate from (i.e. coarse level)
       */
      virtual void interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                       shared_ptr<const ISweeper<time>> src);

      /**
       * interpolate, in time and space, from the coarse sweeper to the fine sweeper.
       *
       * @param[in,out] dst            sweeper to interpolate onto (i.e. fine level)
       * @param[in]     src            sweeper to interpolate from (i.e. coarse level)
       * @param[in]     interp_initial `true` if a delta for the initial condtion should also be 
       *   computed (PFASST)
       */
      virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                               shared_ptr<const ISweeper<time>> src,
                               bool interp_initial = false) = 0;
      //! @}

      //! @{
      /**
       * restrict initial condition from the fine sweeper to the coarse sweeper.
       *
       * @param[in,out] dst sweeper to restrict onto (i.e. coarse level)
       * @param[in]     src sweeper to restrict from (i.e. fine level)
       */
      virtual void restrict_initial(shared_ptr<ISweeper<time>> dst,
                                    shared_ptr<const ISweeper<time>> src);


      /**
       * restrict, in time and space, from the fine sweeper to the coarse sweeper.
       *
       * @param[in,out] dst              sweeper to restrict onto (i.e. coarse level)
       * @param[in]     src              sweeper to restrict from (i.e. fine level)
       * @param[in]     restrict_initial `true` if the initial condition should also be restricted
       */
      virtual void restrict(shared_ptr<ISweeper<time>> dst,
                            shared_ptr<const ISweeper<time>> src,
                            bool restrict_initial = false) = 0;
      //! @}

      //! @{
      /**
       * compute FAS correction between the coarse and fine sweepers.
       *
       * @param[in]     dt  width of the time step to compute FAS correction for
       * @param[in,out] dst sweeper to compute FAS correction for (i.e. coarse level)
       * @param[in]     src sweeper to compute FAS correction from (i.e. fine level)
       */
      virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                       shared_ptr<const ISweeper<time>> src) = 0;
      //! @}
  };
}  // ::pfasst

#include "pfasst/interfaces_impl.hpp"

#endif
