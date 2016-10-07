/**
 * @file pfasst/globals.hpp
 * @since v0.2.0
 */
#ifndef _GLOBALS__HPP_
#define _GLOBALS__HPP_

/**
 * Denote unused function parameters for omitting compiler warnings.
 *
 * To denote an unused function parameter use this macro:
 * @code{.cpp}
 * void foo(T bar) {
 *   UNUSED(bar);
 *   // some logic not using parameter `bar`
 * }
 * @endcode
 *
 * @param[in] expr parameter to be denoted as unused
 * @since v0.2.0
 * @ingroup Utilities
 */
#define UNUSED(expr) \
  (void)(expr)

#endif  // _GLOBALS__HPP_
