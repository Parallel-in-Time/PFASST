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

namespace pfasst {

  class NotImplementedYet : public exception {
    string msg;
  public:
    NotImplementedYet(string msg) : msg(msg) { }
    const char *what() const throw() {
      return (string("Not implemented/supported yet, required for: ")
	      + this->msg).c_str();
    }
  };

  class ISweeper {
  public:
    virtual void setup() { }
    virtual ~ISweeper() { }

    virtual void sweep(double t, double dt) = 0;
    virtual void predict(double t, double dt) = 0;
    virtual void integrate(double t, double dt) = 0;
    virtual void residual(double t, double dt) = 0;
  };

  class ITransfer {
  public:
    // XXX: pass level iterator to these routines as well
    virtual void interpolate(ISweeper *dst, const ISweeper *src) = 0;
    virtual void restrict(ISweeper *dst, const ISweeper *src) = 0;
  };

  class ICommunicator {
  public:
    virtual void post() { }
    virtual void send() = 0;
    virtual void recv() = 0;
  };

}

#endif
