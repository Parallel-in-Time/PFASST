#include "pfasst/exceptions.hpp"


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
  const char* NotImplementedYet::what() const noexcept
  {
    string message = "Not implemented/supported yet, required for: ";
    message += string(runtime_error::what());
    return message.c_str();
  }


  ValueError::ValueError(const string& msg)
    : invalid_argument(msg)
  {}

  /**
   * @internals
   * The message string is prepended by the string `Value Error: `
   * @endinternals
   */
  const char* ValueError::what() const noexcept
  {
    string message = "Value Error: ";
    message += string(invalid_argument::what());
    return message.c_str();
  }
}  // ::pfasst
