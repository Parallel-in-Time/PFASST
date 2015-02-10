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
   * base SDC/MLSDC/PFASST controller.
   *
   * Base controller (see also SDC, MLSDC, and PFASST controllers).
   *
   * @tparam time time precision;
   *     defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class Controller
  {
    protected:
      //! @{
      deque<shared_ptr<ISweeper<time>>>  levels;
      deque<shared_ptr<ITransfer<time>>> transfer;
      //! @}

      //! @{
      size_t step, iteration, max_iterations;
      time t, dt, tend;
      //! @}

    public:
      Controller();
      virtual ~Controller();

      //! @{
      virtual void set_options(bool all_sweepers = true);
      virtual void setup();
      virtual void set_duration(time t0, time tend, time dt, size_t niters);
      virtual void add_level(shared_ptr<ISweeper<time>> swpr,
                             shared_ptr<ITransfer<time>> trnsfr = shared_ptr<ITransfer<time>>(nullptr),
                             bool coarse = true);
      //! @}

      //! @{
      virtual size_t nlevels();

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_level(size_t level)
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(levels[level]);
        assert(r);
        return r;
      }

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_finest()
      {
        return get_level<R>(nlevels() - 1);
      }

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_coarsest()
      {
        return get_level<R>(0);
      }

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
       * Get current time step number.
       */
      virtual size_t get_step();
      virtual void   set_step(size_t n);
      virtual time   get_time_step();
      virtual time   get_time();
      virtual void   advance_time(size_t nsteps = 1);
      virtual time   get_end_time();
      virtual size_t get_iteration();
      virtual void   set_iteration(size_t iter);
      virtual void   advance_iteration();
      virtual size_t get_max_iterations();
      //! @}

      /**
       * level (MLSDC/PFASST) iterator.
       *
       * This iterator is used to walk through the MLSDC/PFASST hierarchy of sweepers.
       * It keeps track of the _current_ level, and has convenience routines to return the
       * LevelIter::current(), LevelIter::fine() (i.e. `current+1`), and LevelIter::coarse()
       * (`current-1`) sweepers.
       *
       * Under the hood it satisfies the requirements of std::random_access_iterator_tag, thus
       * implementing a `RandomAccessIterator`.
       */
      class LevelIter
        : iterator<random_access_iterator_tag, shared_ptr<ISweeper<time>>, int,
                   ISweeper<time>*, ISweeper<time>>
      {
        protected:
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
          template<typename R = ISweeper<time>>
          shared_ptr<R> current()
          {
            return ts->template get_level<R>(level);
          }

          template<typename R = ISweeper<time>>
          shared_ptr<R> fine()
          {
            return ts->template get_level<R>(level + 1);
          }

          template<typename R = ISweeper<time>>
          shared_ptr<R> coarse()
          {
            return ts->template get_level<R>(level - 1);
          }

          template<typename R = ITransfer<time>>
          shared_ptr<R> transfer()
          {
            return ts->template get_transfer<R>(level);
          }
          //! @}

          //! @{
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
      virtual LevelIter finest();
      virtual LevelIter coarsest();
      //! @}
  };
}  // ::pfasst

#include "pfasst/controller/interface_impl.hpp"

#endif
