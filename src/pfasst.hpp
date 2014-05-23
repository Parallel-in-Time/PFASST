/*
 * Interfaces for SDC/MLSDC/PFASST algorithms.
 */

#ifndef _PFASST_HPP_
#define _PFASST_HPP_

#include <deque>
#include <exception>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>

using namespace std;

namespace pfasst {

  //
  // sdc sweeper interface
  //

  struct isweeper {
    virtual void setup() { }
    virtual ~isweeper() { }

    // required for all sdc schemes
    virtual void sweep(double t, double dt) = 0;
    virtual void predict(double t, double dt) = 0;

    // required for multi-level sdc and pfasst schemes
    virtual void interpolate(const isweeper*) { }
    virtual void restrict(const isweeper*) { }

    // required for pfasst schemes
    virtual void post() { }
    virtual void send() { }
    virtual void recv() { }
  };

  //
  // pfasst controller
  //

  class pfasst {
    deque<shared_ptr<isweeper>> levels;

  public:
    int nstep, niter;
    double dt, t;

    void add_level(isweeper *sweeper, bool coarse) {
      if (coarse)
	levels.push_front(shared_ptr<isweeper>(sweeper));
      else
	levels.push_back(shared_ptr<isweeper>(sweeper));
    }

    template<typename R=isweeper> R* get_level(int level) {
      return dynamic_cast<R*>(levels[level].get());
    }

    int nlevels() { return levels.size(); }
    void setup();
    void run();

    void predictor();
    void iteration();

    struct leveliter {
      int level;
      pfasst *pf;

      leveliter(int level, pfasst *pf) : level(level), pf(pf) {}

      template<typename R=isweeper> R* current() {
	return pf->get_level<R>(level);
      }
      template<typename R=isweeper> R* fine() {
	return pf->get_level<R>(level+1);
      }
      template<typename R=isweeper> R* coarse() {
	return pf->get_level<R>(level-1);
      }

      isweeper *operator*() { return current(); }
      bool operator==(leveliter i) { return level == i.level; }
      bool operator!=(leveliter i) { return level != i.level; }
      bool operator<=(leveliter i) { return level <= i.level; }
      bool operator>=(leveliter i) { return level >= i.level; }
      bool operator< (leveliter i) { return level <  i.level; }
      bool operator> (leveliter i) { return level >  i.level; }
      leveliter operator- (int i) { return leveliter(level-1, pf); }
      leveliter operator+ (int i) { return leveliter(level+1, pf); }
      void operator++() { level++; }
      void operator--() { level--; }
    };

    leveliter finest()   { return leveliter(nlevels()-1, this); }
    leveliter coarsest() { return leveliter(0, this); }

    leveliter cycle_down(leveliter levels, double t, double dt);
    leveliter cycle_up(leveliter levels, double t, double dt);
    leveliter cycle_bottom(leveliter levels, double t, double dt);
    leveliter cycle_top(leveliter levels, double t, double dt);
    leveliter cycle_v(leveliter levels, double t, double dt);
  };

}

#endif
