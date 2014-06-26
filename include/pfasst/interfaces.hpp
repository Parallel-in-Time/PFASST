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

  using time_precision = double;

  class NotImplementedYet : public exception {
    string msg;
  public:
    NotImplementedYet(string msg) : msg(msg) { }
    const char *what() const throw() {
      return (string("Not implemented/supported yet, required for: ")
	      + this->msg).c_str();
    }
  };

  template<typename time=time_precision>
  class ISweeper {
  public:
    virtual ~ISweeper() { }
    virtual void setup(bool coarse=false) { }

    virtual void sweep(time t, time dt) = 0;
    virtual void predict(time t, time dt, bool initial) = 0;
    virtual void advance() = 0;
    virtual void save() { NotImplementedYet("mlsdc/pfasst"); }
  };

  template<typename time=time_precision>
  class ITransfer {
  public:
    // XXX: pass level iterator to these routines as well
    // XXX: these needs to be templated
    virtual ~ITransfer() { }
    virtual void interpolate(ISweeper<time> *dst, const ISweeper<time> *src,
			     bool interp_delta_from_initial=false,
			     bool interp_initial=false) = 0;
    virtual void restrict(ISweeper<time> *dst, const ISweeper<time> *src,
			  bool restrict_initial=false) = 0;
    virtual void fas(time dt, ISweeper<time> *dst, const ISweeper<time> *src) = 0;
  };

  class ICommunicator {
  public:
    virtual void post() { }
    virtual void send() = 0;
    virtual void recv() = 0;
  };

}

#endif
