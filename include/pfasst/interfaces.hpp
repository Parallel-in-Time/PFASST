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

  class ICommunicator {
  public:
    virtual ~ICommunicator() { }
    virtual int size() = 0;
    virtual int rank() = 0;
  };

  class ISweeper {
  public:
    virtual ~ISweeper() { }
    virtual void setup(bool coarse=false) { }

    virtual void sweep(double t, double dt) = 0; // XXX: this needs to be a templated
    virtual void predict(double t, double dt, bool initial) = 0; // XXX: this needs to be templated
    virtual void advance() = 0;
    virtual void save(bool initial_only=false) { NotImplementedYet("mlsdc/pfasst"); }
    virtual void post(ICommunicator* comm) { };
    virtual void send(ICommunicator* comm) { NotImplementedYet("pfasst"); }
    virtual void recv(ICommunicator* comm) { NotImplementedYet("pfasst"); }
  };

  class ITransfer {
  public:
    // XXX: pass level iterator to these routines as well
    // XXX: these needs to be templated
    virtual ~ITransfer() { }
    virtual void interpolate(ISweeper *dst, const ISweeper *src,
			     bool interp_delta_from_initial=false,
			     bool interp_initial=false) = 0;
    virtual void restrict(ISweeper *dst, const ISweeper *src,
			  bool restrict_initial=false,
			  bool restrict_initial_only=false) = 0;
    virtual void fas(double dt, ISweeper *dst, const ISweeper *src) = 0;
  };

}

#endif
