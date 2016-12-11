#include "pfasst/interfaces.hpp"

#include <cassert>
#include <memory>
#include <string>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  NotImplementedYet::NotImplementedYet(const string& m)
    : runtime_error(m)
  {
    this->msg = string("Not implemented/supported yet, required for: ") + string(runtime_error::what());
  }

  const char* NotImplementedYet::what() const throw()
  {
    return this->msg.c_str();
  }


  ValueError::ValueError(const string& msg)
    : invalid_argument(msg)
  {
    this->msg = string("ValueError: ") + string(invalid_argument::what());
  }

  const char* ValueError::what() const throw()
  {
    return this->msg.c_str();
  }


  ICommunicator::~ICommunicator()
  {}


  IStatus::~IStatus()
  {}

  void IStatus::set_comm(ICommunicator* comm)
  {
    this->comm = comm;
  }

  bool IStatus::previous_is_iterating()
  {
    if (this->comm->rank() == 0) {
      return false;
    }
    return !this->get_converged(this->comm->rank() - 1);
  }

  bool IStatus::keep_iterating()
  {
    if (this->comm->rank() == 0) {
      return !this->get_converged(0);
    }
    bool keep_iterating = !this->get_converged(this->comm->rank() - 1) || !this->get_converged(this->comm->rank());
    ML_CLOG(DEBUG, "Controller", "previous converged: " << boolalpha << this->get_converged(this->comm->rank() - 1)
                                 << "; this converged: " << boolalpha << this->get_converged(this->comm->rank())
                                 << " --> keep iterating: " << boolalpha << keep_iterating);
    return keep_iterating;
  }


  template<typename time>
  ISweeper<time>::ISweeper()
    : controller(nullptr)
  {}

  template<typename time>
  ISweeper<time>::~ISweeper()
  {}

  template<typename time>
  void ISweeper<time>::set_controller(Controller<time>* ctrl)
  {
    this->controller = ctrl;
  }

  template<typename time>
  Controller<time>* ISweeper<time>::get_controller()
  {
    assert(this->controller);
    return this->controller;
  }

  template<typename time>
  void ISweeper<time>::setup(bool coarse)
  {
    UNUSED(coarse);
  }

  template<typename time>
  bool ISweeper<time>::converged()
  {
    return false;
  }

  template<typename time>
  void ISweeper<time>::save(bool initial_only)
  {
    UNUSED(initial_only);
    throw NotImplementedYet("mlsdc/pfasst");
  }

  template<typename time>
  void ISweeper<time>::spread()
  {
    throw NotImplementedYet("pfasst");
  }

  template<typename time>
  void ISweeper<time>::post_sweep()
  {}

  template<typename time>
  void ISweeper<time>::post_predict()
  {}

  template<typename time>
  void ISweeper<time>::post_step()
  {}

  template<typename time>
  void ISweeper<time>::post(ICommunicator* comm, int tag)
  {
    UNUSED(comm); UNUSED(tag);
  }

  template<typename time>
  void ISweeper<time>::send(ICommunicator* comm, int tag, bool blocking)
  {
    UNUSED(comm); UNUSED(tag); UNUSED(blocking);
    throw NotImplementedYet("pfasst");
  }

  template<typename time>
  void ISweeper<time>::recv(ICommunicator* comm, int tag, bool blocking)
  {
    UNUSED(comm); UNUSED(tag); UNUSED(blocking);
    throw NotImplementedYet("pfasst");
  }

  template<typename time>
  void ISweeper<time>::broadcast(ICommunicator* comm)
  {
    UNUSED(comm);
    throw NotImplementedYet("pfasst");
  }


  template<typename time>
  ITransfer<time>::~ITransfer()
  {}

  template<typename time>
  void ITransfer<time>::interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                            shared_ptr<const ISweeper<time>> src)
  {
    UNUSED(dst); UNUSED(src);
    throw NotImplementedYet("pfasst");
  }

  template<typename time>
  void ITransfer<time>::restrict_initial(shared_ptr<ISweeper<time>> dst,
                                         shared_ptr<const ISweeper<time>> src)
  {
    UNUSED(dst); UNUSED(src);
    throw NotImplementedYet("pfasst");
  }

}  // ::pfasst
