/*
 * Base controller (see also SDC, MLSDC, and PFASST controllers).
 */

#ifndef _PFASST_CONTROLLER_HPP_
#define _PFASST_CONTROLLER_HPP_

#include "interfaces.hpp"

namespace pfasst {

  // class LevelIter;

  template<typename time>
  class Controller {
  protected:
    deque<shared_ptr<ISweeper>>  levels;
    deque<shared_ptr<ITransfer>> transfer;

    int    nsteps, niters;
    time  dt;

  public:

    void setup() {
      for (auto l=coarsest(); l<=finest(); ++l) {
	l.current()->setup();
      }
    }

    void set_duration(time dt, int nsteps, int niters) {
      this->dt = dt; this->nsteps = nsteps; this->niters = niters;
    }

    void add_level(ISweeper *swpr, ITransfer *trnsfr=NULL, bool coarse=true) {
      if (coarse) {
	levels.push_front(shared_ptr<ISweeper>(swpr));
	transfer.push_front(shared_ptr<ITransfer>(trnsfr));
      } else {
	levels.push_back(shared_ptr<ISweeper>(swpr));
	transfer.push_back(shared_ptr<ITransfer>(trnsfr));
      }
    }

    template<typename R=ISweeper> R* get_level(int level) {
      return dynamic_cast<R*>(levels[level].get());
    }

    template<typename R=ITransfer> R* get_transfer(int level) {
      return dynamic_cast<R*>(transfer[level].get());
    }

    int nlevels() {
      return levels.size();
    }


    class LevelIter {
      Controller *ts;

    public:
      int level;

      LevelIter(int level, Controller *ts) : ts(ts), level(level) {}

      template<typename R=ISweeper> R* current() {
	return ts->get_level<R>(level);
      }
      template<typename R=ISweeper> R* fine() {
	return ts->get_level<R>(level+1);
      }
      template<typename R=ISweeper> R* coarse() {
	return ts->get_level<R>(level-1);
      }
      template<typename R=ITransfer> R* transfer() {
	return ts->get_transfer<R>(level);
      }

      ISweeper *operator*() { return current(); }
      bool operator==(LevelIter i) { return level == i.level; }
      bool operator!=(LevelIter i) { return level != i.level; }
      bool operator<=(LevelIter i) { return level <= i.level; }
      bool operator>=(LevelIter i) { return level >= i.level; }
      bool operator< (LevelIter i) { return level <  i.level; }
      bool operator> (LevelIter i) { return level >  i.level; }
      LevelIter operator- (int i) { return LevelIter(level-1, ts); }
      LevelIter operator+ (int i) { return LevelIter(level+1, ts); }
      void operator++() { level++; }
      void operator--() { level--; }
    };

    LevelIter finest()   { return LevelIter(nlevels()-1, this); }
    LevelIter coarsest() { return LevelIter(0, this); }

  };

}

#endif
