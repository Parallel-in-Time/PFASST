/*
 * Base controller (see also SDC, MLSDC, and PFASST controllers).
 */

#ifndef _PFASST_CONTROLLER_HPP_
#define _PFASST_CONTROLLER_HPP_

#include <deque>
#include <memory>
#include <cassert>
#include <iterator>

#include "interfaces.hpp"

namespace pfasst
{
  /**
   * base SDC/MLSDC/PFASST controller.
   * @tparam time time precision
   *     defaults to pfasst::time_precision
   */
  template<typename time = time_precision>
  class Controller
  {
    protected:
      deque<shared_ptr<ISweeper<time>>>  levels;
      deque<shared_ptr<ITransfer<time>>> transfer;

      size_t step, iteration, max_iterations;
      time t, dt, tend;

    public:
      //! @{
      virtual void setup()
      {
        for (auto l = coarsest(); l <= finest(); ++l) {
          l.current()->set_controller(this);
          l.current()->setup();
        }
      }

      void set_duration(time t0, time tend, time dt, size_t niters)
      {
        this->t = t0;
        this->tend = tend;
        this->dt = dt;
        this->step = 0;
        this->iteration = 0;
        this->max_iterations = niters;
      }

      void add_level(shared_ptr<ISweeper<time>> swpr,
                     shared_ptr<ITransfer<time>> trnsfr = shared_ptr<ITransfer<time>>(nullptr),
                     bool coarse = true)
      {
        if (coarse) {
          levels.push_front(swpr);
          transfer.push_front(trnsfr);
        } else {
          levels.push_back(swpr);
          transfer.push_back(trnsfr);
        }
      }
      //! @}

      //! @{
      template<typename R = ISweeper<time>>
      shared_ptr<R> get_level(size_t level)
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(levels[level]); assert(r);
        return r;
      }

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_finest()
      {
        return get_level<R>(nlevels()-1);
      }

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_coarsest()
      {
        return get_level<R>(0);
      }

      template<typename R = ITransfer<time>>
      shared_ptr<R> get_transfer(size_t level)
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(transfer[level]); assert(r);
        return r;
      }

      size_t nlevels()
      {
        return levels.size();
      }
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
        : iterator<random_access_iterator_tag, shared_ptr<ISweeper<time>>, size_t,
                   ISweeper<time>*, ISweeper<time>>
      {
          Controller* ts;

        public:
          typedef size_t                     difference_type;
          typedef shared_ptr<ISweeper<time>> value_type;
          typedef ISweeper<time>*            pointer;
          typedef ISweeper<time>             reference;
          typedef random_access_iterator_tag iterator_category;

          size_t level;

          //! @{
          LevelIter(size_t level, Controller* ts)
            : ts(ts), level(level)
          {}
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
          shared_ptr<R> operator*()                { return current<R>(); }
          LevelIter  operator++()                  { level++; return *this; }
          // required by std::input_iterator_tag
          template<typename R = reference>
          shared_ptr<R> operator->()               { return current<R>(); }
          bool       operator==(LevelIter i)       { return level == i.level; }
          bool       operator!=(LevelIter i)       { return level != i.level; }
          // required by std::bidirectional_iterator_tag
          LevelIter  operator--()                  { level--; return *this; }
          // required by std::random_access_iterator_tag
          LevelIter  operator- (difference_type i) { return LevelIter(level - i, ts); }
          LevelIter  operator+ (difference_type i) { return LevelIter(level + i, ts); }
          bool       operator<=(LevelIter i)       { return level <= i.level; }
          bool       operator>=(LevelIter i)       { return level >= i.level; }
          bool       operator< (LevelIter i)       { return level <  i.level; }
          bool       operator> (LevelIter i)       { return level >  i.level; }
          //! @}
      };

      //! @{
      LevelIter finest()   { return LevelIter(nlevels() - 1, this); }
      LevelIter coarsest() { return LevelIter(0, this); }
      //! @}


      /**
       * Get current time step number.
       */
      size_t get_step()
      {
        return step;
      }

      time get_time_step()
      {
        return dt;
      }

      time get_time()
      {
        return t;
      }

      void advance_time(size_t nsteps=1)
      {
        step += nsteps;
        t += nsteps*dt;
      }

      time get_end_time()
      {
        return tend;
      }

      size_t get_iteration()
      {
        return iteration;
      }

      void set_iteration(size_t iter)
      {
        this->iteration = iter;
      }

      void advance_iteration()
      {
        iteration++;
      }

      size_t get_max_iterations()
      {
        return max_iterations;
      }

  };
}

#endif
