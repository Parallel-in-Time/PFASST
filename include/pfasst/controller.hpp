/*
 * Base controller (see also SDC, MLSDC, and PFASST controllers).
 */

#ifndef _PFASST_CONTROLLER_HPP_
#define _PFASST_CONTROLLER_HPP_

#include <deque>
#include <memory>
#include <cassert>

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

      template<typename R = ISweeper<time>>
      shared_ptr<R> get_level(size_t level)
      {
        shared_ptr<R> r = dynamic_pointer_cast<R>(levels[level]);
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

      /**
       * level (MLSDC/PFASST) iterator.
       * 
       * This iterator is used to walk through the MLSDC/PFASST hierarchy of sweepers.
       * It keeps track of the _current_ level, and has convenience routines to return the 
       * LevelIter::current(), LevelIter::fine() (i.e. `current+1`), and LevelIter::coarse()
       * (`current-1`) sweepers.
       */
      class LevelIter
      {
          Controller* ts;

        public:
          size_t level;

          LevelIter(size_t level, Controller* ts) : ts(ts), level(level) {}

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

          shared_ptr<ISweeper<time>> operator*() { return current(); }
          bool operator==(LevelIter i) { return level == i.level; }
          bool operator!=(LevelIter i) { return level != i.level; }
          bool operator<=(LevelIter i) { return level <= i.level; }
          bool operator>=(LevelIter i) { return level >= i.level; }
          bool operator< (LevelIter i) { return level <  i.level; }
          bool operator> (LevelIter i) { return level >  i.level; }
          LevelIter operator- (size_t i) { return LevelIter(level - i, ts); }
          LevelIter operator+ (size_t i) { return LevelIter(level + i, ts); }
          void operator++() { level++; }
          void operator--() { level--; }
      };

      LevelIter finest()   { return LevelIter(nlevels() - 1, this); }
      LevelIter coarsest() { return LevelIter(0, this); }

  };

}

#endif
