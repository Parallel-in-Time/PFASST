#include "pfasst/interfaces.hpp"

#include <cassert>
#include <memory>
#include <string>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  NotImplementedYet::NotImplementedYet(const string& msg)
    : runtime_error(msg)
  {}

  /**
   * @internals
   * The message string is prepended by the string `Not implemented/supported yet, required for: `
   * @endinternals
   */
  const char* NotImplementedYet::what() const throw()
  {
    return (string("Not implemented/supported yet, required for: ") + string(runtime_error::what())).c_str();
  }


  ValueError::ValueError(const string& msg)
    : invalid_argument(msg)
  {}

  /**
   * @internals
   * The message string is prepended by the string `ValueError: `
   * @endinternals
   */
  const char* ValueError::what() const throw()
  {
    return (string("ValueError: ") + string(invalid_argument::what())).c_str();
  }


  ICommunicator::~ICommunicator()
  {}


  IStatus::~IStatus()
  {}

  //! @todo Consider asserting validity of the given pointer to the communicator.
  void IStatus::set_comm(ICommunicator* comm)
  {
    this->comm = comm;
  }

  //! @todo Consider asserting presence of valid communicator.
  bool IStatus::previous_is_iterating()
  {
    if (this->comm->rank() == 0) {
      return false;
    }
    return !this->get_converged(this->comm->rank() - 1);
  }

  /**
   * @internals
   * Returning logic depends on current process' rank.
   * In case it is not the master process, both the converged state of this and the previous rank
   * are checked.
   *
   * @see `IStatus::get_converged()`
   * @endinternals
   *
   * @todo Consider asserting presence of valid communicator.
   */
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

  /**
   * @todo Consider asserting presence of the given pointer to the controller.
   */
  template<typename time>
  void ISweeper<time>::set_controller(Controller<time>* ctrl)
  {
    this->controller = ctrl;
  }

  /**
   * @internals
   * @note Asserts presense of a controller if `NDEBUG` is not defined.
   * @endinternals
   */
  template<typename time>
  Controller<time>* ISweeper<time>::get_controller()
  {
    assert(this->controller);
    return this->controller;
  }

  /**
   * @internals
   * @note Unless overwritten by implementations, this is a no-op.
   * @endinternals
   */
  template<typename time>
  void ISweeper<time>::set_options()
  {}

  template<typename time>
  void ISweeper<time>::setup(bool coarse)
  {
    UNUSED(coarse);
  }

  /**
   * @returns Unless overwritten by implementations, this will always return `false`.
   */
  template<typename time>
  bool ISweeper<time>::converged()
  {
    return false;
  }

  /**
   * @throws NotImplementedYet This function is required by MLSDC and PFASST
   */
  template<typename time>
  void ISweeper<time>::save(bool initial_only)
  {
    UNUSED(initial_only);
    throw NotImplementedYet("mlsdc/pfasst");
  }

  /**
   * @throws NotImplementedYet This function is required by PFASST
   */
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

  /**
   * @throws NotImplementedYet This function is required by PFASST
   */
  template<typename time>
  void ISweeper<time>::send(ICommunicator* comm, int tag, bool blocking)
  {
    UNUSED(comm); UNUSED(tag); UNUSED(blocking);
    throw NotImplementedYet("pfasst");
  }

  /**
   * @throws NotImplementedYet This function is required by PFASST
   */
  template<typename time>
  void ISweeper<time>::recv(ICommunicator* comm, int tag, bool blocking)
  {
    UNUSED(comm); UNUSED(tag); UNUSED(blocking);
    throw NotImplementedYet("pfasst");
  }

  /**
   * @throws NotImplementedYet This function is required by PFASST
   */
  template<typename time>
  void ISweeper<time>::broadcast(ICommunicator* comm)
  {
    UNUSED(comm);
    throw NotImplementedYet("pfasst");
  }


  template<typename time>
  ITransfer<time>::~ITransfer()
  {}

  /**
   * @throws NotImplementedYet This function is required by PFASST
   */
  template<typename time>
  void ITransfer<time>::interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                            shared_ptr<const ISweeper<time>> src)
  {
    UNUSED(dst); UNUSED(src);
    throw NotImplementedYet("pfasst");
  }

  /**
   * @throws NotImplementedYet This function is required by PFASST
   */
  template<typename time>
  void ITransfer<time>::restrict_initial(shared_ptr<ISweeper<time>> dst,
                                         shared_ptr<const ISweeper<time>> src)
  {
    UNUSED(dst); UNUSED(src);
    throw NotImplementedYet("pfasst");
  }
}  // ::pfasst
