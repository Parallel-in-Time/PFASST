/*
 * Base controller (see also SDC, MLSDC, and PFASST controllers).
 */

#ifndef _PFASST_CONTROLLER_HPP_
#define _PFASST_CONTROLLER_HPP_

#include "pfasst-interfaces.hpp"

namespace pfasst {

  template<typename timeT>
  class Controller {
  protected:
    deque<shared_ptr<ISweeper>> levels;

    int    nsteps, niters;
    timeT  dt;

  public:

    void setup() {
      for (auto l=coarsest(); l<=finest(); ++l) {
	l.current()->setup();
      }
    }

    void set_duration(timeT dt, int nsteps, int niters) {
      this->dt = dt; this->nsteps = nsteps; this->niters = niters;
    }

    void add_level(ISweeper *sweeper, bool coarse=true) {
      if (coarse)
	levels.push_front(shared_ptr<ISweeper>(sweeper));
      else
	levels.push_back(shared_ptr<ISweeper>(sweeper));
    }

    template<typename R=ISweeper> R* get_level(int level) {
      return dynamic_cast<R*>(levels[level].get());
    }

    int nlevels() {
      return levels.size();
    }

    struct leveliter {
      int level;
      Controller *ts;

      leveliter(int level, Controller *ts) : level(level), ts(ts) {}

      template<typename R=ISweeper> R* current() {
	return ts->get_level<R>(level);
      }
      template<typename R=ISweeper> R* fine() {
	return ts->get_level<R>(level+1);
      }
      template<typename R=ISweeper> R* coarse() {
	return ts->get_level<R>(level-1);
      }

      ISweeper *operator*() { return current(); }
      bool operator==(leveliter i) { return level == i.level; }
      bool operator!=(leveliter i) { return level != i.level; }
      bool operator<=(leveliter i) { return level <= i.level; }
      bool operator>=(leveliter i) { return level >= i.level; }
      bool operator< (leveliter i) { return level <  i.level; }
      bool operator> (leveliter i) { return level >  i.level; }
      leveliter operator- (int i) { return leveliter(level-1, ts); }
      leveliter operator+ (int i) { return leveliter(level+1, ts); }
      void operator++() { level++; }
      void operator--() { level--; }
    };

    leveliter finest()   { return leveliter(nlevels()-1, this); }
    leveliter coarsest() { return leveliter(0, this); }

  };

}

#endif
