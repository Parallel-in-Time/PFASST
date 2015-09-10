/**
 * @file controller/interface.hpp
 * @since v0.1.0
 */
#ifndef _PFASST_CONTROLLER_HPP_
#define _PFASST_CONTROLLER_HPP_

#include <cassert>
#include <deque>
#include <iterator>
#include <memory>
using namespace std;

#include "pfasst/interfaces.hpp"


namespace pfasst
{
  /**
   * Base SDC/MLSDC/PFASST controller.
   *
   * Base controller (see also SDC, MLSDC, and PFASST controllers).
   *
   * @tparam time time precision; defaults to pfasst::time_precision
   *
   * @since v0.1.0
   *
   * @ingroup Controllers
   */
  template<typename time = time_precision>
  class Controller
  {
    protected:
      //! @{
      /**
       * Ordered list of all levels.
       *
       * A level is represented by a sweeper implementing the ISweeper interface.
       *
       * @note The levels are ordered from coarsest (index `0`) to finest.
       */
      deque<shared_ptr<ISweeper<time>>>  levels;

      /**
       * Ordered list of transfer operators for levels.
       *
       * A transfer operator for a level must implement the ITransfer interface.
       */
      deque<shared_ptr<ITransfer<time>>> transfer;
      //! @}

      //! @{
      /**
       * Current time step index.
       */
      size_t step;

      /**
       * Current iteration index on current time step.
       */
      size_t iteration;

      /**
       * Maximum iterations per time step.
       */
      size_t max_iterations;

      /**
       * \\( t_0 \\) of current time step.
       */
      time t;

      /**
       * Width of current time step (\\( \\Delta t \\)).
       */
      time dt;

      /**
       * \\( T_{end} \\) of last time step.
       */
      time tend;
      //! @}

    public:
      //! @{
      Controller();
      virtual ~Controller();
      //! @}

      //! @{
      /**
       * Set options from command line etc.
       *
       * @param[in] all_sweepers if given also calls ISweeper::set_options for all already added
       *   levels
       */
      virtual void set_options(bool all_sweepers = true);

      /**
       * Basic setup routine for this controller.
       *
       * @pre It is expected (i.e. the user has to make that sure) that all levels and transfer
       *   operators have been instantiated and added to this controller before calling
       *   Controller::setup().
       */
      virtual void setup();

      /**
       * Set basic time scope of the Controller.
       *
       * @param[in] t0     time start point of the first time step
       * @param[in] tend   time end point of the last time step
       * @param[in] dt     width of one time step
       * @param[in] niters maximum number of iterations per time step
       * The total number of time steps will get computed internally if required.
       */
      virtual void set_duration(time t0, time tend, time dt, size_t niters);

      /**
       * Adding a level to the controller.
       *
       * @param[in] sweeper  sweeper representing the level
       * @param[in] transfer corresponding transfer operator for the level
       * @param[in] coarse   whether to add this level as a coarser one to the list of existing
       */
      virtual void add_level(shared_ptr<ISweeper<time>> sweeper,
                             shared_ptr<ITransfer<time>> transfer = shared_ptr<ITransfer<time>>(nullptr),
                             bool coarse = true);
      //! @}

      //! @{
      /**
       * Total number of levels controlled by this Controller.
       */
      virtual size_t nlevels();

      /**
       * Get sweeper for level with index @p level.
       *
       * @tparam R type of the level sweeper
       * @returns sweeper of type @p R for requested level
       *
       * @internals
       * @note Asserts sweeper for level @p level can be interpreted as @p R if `NDEBUG` is not 
       *   defined.
       * @endinternals
       */
      template<typename R = ISweeper<time>>
      shared_ptr<R> get_level(size_t level)
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(levels[level]);
        assert(r);
        return r;
      }

      /**
       * Get coarsest level.
       *
       * @tparam R type of the level sweeper
       * @see Controller::get_level() with `level=(Controller::nlevels() - 1)`
       */
      template<typename R = ISweeper<time>>
      shared_ptr<R> get_finest()
      {
        return get_level<R>(nlevels() - 1);
      }

      /**
       * Get coarsest level.
       *
       * @tparam R type of the level sweeper
       * @see Controller::get_level() with `level=0`
       */
      template<typename R = ISweeper<time>>
      shared_ptr<R> get_coarsest()
      {
        return get_level<R>(0);
      }

      /**
       * Retreive transfer operator for level @p level.
       *
       * @tparam R type of the requested transfer operator
       * @param[in] level level index to retreive transfer operator for
       * @returns transfer operator of requested type @p R for desired level
       *
       * @internals
       * @note Asserts transfer operator for level @p level can be interpreted as @p R if `NDEBUG`
       *   is not defined.
       * @endinternals
       */
      template<typename R = ITransfer<time>>
      shared_ptr<R> get_transfer(size_t level)
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(transfer[level]);
        assert(r);
        return r;
      }
      //! @}

      //! @{
      /**
       * Get current time step index.
       *
       * The time step number is zero-based, i.e. the first time step has index `0`.
       *
       * @returns current time step index
       */
      virtual size_t get_step();

      /**
       * Set current time step index.
       *
       * @param[in] n index of new time step
       */
      virtual void   set_step(size_t n);

      /**
       * Get width of current time step.
       *
       * @returns \\( \\Delta t \\) of current time step
       */
      virtual time   get_step_size();

      /**
       * Get width of current time step (alias for get_step_size).
       *
       * @returns \\( \\Delta t \\) of current time step
       */
      time get_dt() { return this->get_step_size(); }

      /**
       * Get start time point of current time step.
       *
       * @returns \\( t_0 \\) of current time step
       */
      virtual time   get_time();

      /**
       * Get start time point of current time step (alias for get_time).
       *
       * @returns \\( t_0 \\) of current time step
       */
      time get_t() { return this->get_time(); }

      /**
       * Advance to a following time step
       *
       * @param[in] nsteps number of time steps to advance; `1` meaning the next step
       */
      virtual void   advance_time(size_t nsteps = 1);

      /**
       * Get end time point of last time step.
       *
       * @returns \\( T_{end} \\) of last time step.
       */
      virtual time   get_end_time();

      /**
       * Get current iteration index of current time step.
       *
       * @returns current iteration index of current time step.
       */
      virtual size_t get_iteration();

      /**
       * Set current iteration of current time step.
       *
       * @param[in] iter iteration index to set
       */
      virtual void   set_iteration(size_t iter);

      /**
       * Advance to the next iteration.
       *
       * This method may or may not trigger additional post-iteration procedures.
       */
      virtual void   advance_iteration();

      /**
       * Get maximum number of allowed iterations per time step.
       *
       * @returns maximum allowed iterations per time step.
       */
      virtual size_t get_max_iterations();
      //! @}

      /**
       * Level (MLSDC/PFASST) iterator.
       *
       * This iterator is used to walk through the MLSDC/PFASST hierarchy of sweepers.
       * It keeps track of the _current_ level, and has convenience routines to return the
       * LevelIter::current(), LevelIter::fine() (i.e. `current+1`), and LevelIter::coarse()
       * (`current-1`) sweepers.
       *
       * @since v0.1.0
       *
       * @internals
       * Under the hood it satisfies the requirements of
       * [std::random_access_iterator_tag](http://en.cppreference.com/w/cpp/iterator/iterator_tags),
       * thus implementing a
       * [RandomAccessIterator](http://en.cppreference.com/w/cpp/concept/RandomAccessIterator).
       * @endinternals
       */
      class LevelIter
        : std::iterator<random_access_iterator_tag, shared_ptr<ISweeper<time>>, int,
                        ISweeper<time>*, ISweeper<time>>
      {
        protected:
          /**
           * Controller this iterator is bound to.
           */
          Controller* ts;

        public:
          //! @{
          typedef int                        difference_type;
          typedef shared_ptr<ISweeper<time>> value_type;
          typedef ISweeper<time>*            pointer;
          typedef ISweeper<time>             reference;
          typedef random_access_iterator_tag iterator_category;
          //! @}

          //! @{
          int level;
          //! @}

          //! @{
          LevelIter(int level, Controller* ts);
          //! @}

          //! @{
          /**
           * Get level this iterator is currently pointing at.
           *
           * @tparam R type implementing ISweeper
           * @returns sweeper of type @p R this is currently pointing at
           */
          template<typename R = ISweeper<time>>
          shared_ptr<R> current()
          {
            return ts->template get_level<R>(level);
          }

          /**
           * Get the next finer level based on LevelIter::current()
           *
           * @tparam R type implementing ISweeper
           * @returns next finer sweeper with respect to current()
           */
          template<typename R = ISweeper<time>>
          shared_ptr<R> fine()
          {
            return ts->template get_level<R>(level + 1);
          }

          /**
           * Get the next coarser level based on LevelIter::current()
           *
           * @tparam R type implementing ISweeper
           * @returns next coarser sweeper with respect to current()
           */
          template<typename R = ISweeper<time>>
          shared_ptr<R> coarse()
          {
            return ts->template get_level<R>(level - 1);
          }

          /**
           * Get transfer operator for current level.
           *
           * @tparam R type implementing ITransfer
           * @returns transfer operator for current level
           */
          template<typename R = ITransfer<time>>
          shared_ptr<R> transfer()
          {
            return ts->template get_transfer<R>(level);
          }
          //! @}

          /**
           * @name Standard Operators for Iterators
           * @{
           */
          // required by std::iterator
          template<typename R = reference>
          shared_ptr<R> operator*()
          {
            return current<R>();
          }

          // required by std::input_iterator_tag
          template<typename R = reference>
          shared_ptr<R> operator->()
          {
            return current<R>();
          }

          virtual LevelIter     operator++();
          virtual bool          operator==(LevelIter i);
          virtual bool          operator!=(LevelIter i);
          // required by std::bidirectional_iterator_tag
          virtual LevelIter     operator--();
          // required by std::random_access_iterator_tag
          virtual LevelIter     operator- (difference_type i);
          virtual LevelIter     operator+ (difference_type i);
          virtual bool          operator<=(LevelIter i);
          virtual bool          operator>=(LevelIter i);
          virtual bool          operator< (LevelIter i);
          virtual bool          operator> (LevelIter i);
          //! @}
      };

      //! @{
      /**
       * Convenience accessor to the finest level.
       *
       * @returns iterator to the finest level.
       */
      virtual LevelIter finest();

      /**
       * Convenience accessor to the coarsest level.
       *
       * @returns iterator to the coarsest level.
       */
      virtual LevelIter coarsest();
      //! @}
  };
}  // ::pfasst

#include "pfasst/controller/interface_impl.hpp"

#endif
