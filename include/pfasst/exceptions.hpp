/**
 * @file pfasst/exceptions.hpp
 * @since v0.6.0
 */
#ifndef _PFASST__EXCEPTIONS_HPP_
#define _PFASST__EXCEPTIONS_HPP_

#include <stdexcept>
#include <string>
using namespace std;


namespace pfasst
{
  /**
   * Not implemented yet exception.
   *
   * Used by PFASST to mark methods that are required for a particular algorithm (SDC/MLSDC/PFASST)
   * that may not be necessary for all others.
   *
   * @since v0.1.0
   */
  class NotImplementedYet
    : public runtime_error
  {
    public:
      /**
       * @param[in] msg component or algorithm the throwing function is required for
       */
      explicit NotImplementedYet(const string& msg);
      virtual const char* what() const noexcept;
  };


  /**
   * Value exception.
   *
   * Thrown when a PFASST routine is passed an invalid value.
   *
   * @since v0.1.0
   *
   * @todo Consider deprecating this in favour of std::invalid_argument.
   */
  class ValueError
    : public invalid_argument
  {
    public:
      explicit ValueError(const string& msg);
      virtual const char* what() const noexcept;
  };
}  // ::pfasst

#include "pfasst/exceptions_impl.hpp"

#endif  // _PFASST__EXCEPTIONS_HPP_
