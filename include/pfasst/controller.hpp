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

      size_t nsteps, niters;
      time   dt;

    public:
      //! @{
      void setup()
      {
        for (auto l = coarsest(); l <= finest(); ++l) {
          l.current()->setup();
        }
      }

      void set_duration(time dt, size_t nsteps, size_t niters)
      {
        this->dt = dt;
        this->nsteps = nsteps;
        this->niters = niters;
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
        shared_ptr<R> r = dynamic_pointer_cast<R>(levels[level]);
        assert(r);
        return r;
      }

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_finest()
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(levels.back());
        assert(r);
        return r;
      }

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_coarsest()
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(levels.front());
        assert(r);
        return r;
      }

      template<typename R = ITransfer<time>>
      shared_ptr<R> get_transfer(size_t level)
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(transfer[level]);
        assert(r);
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
  };
}

#endif
