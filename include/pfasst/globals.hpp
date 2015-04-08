/**
 * @file pfasst/globals.hpp
 * @since v0.2.0
 */
#ifndef _GLOBALS__HPP_
#define _GLOBALS__HPP_

/**
 * Denoting unused function parameters for omitting compiler warnings.
 *
 * To denote an unused function parameter just use this makro on it in the function body:
 * @code{.cpp}
 * void foo(T bar) {
 *   UNUSED(bar);
 *   // some logic not using parameter `bar`
 * }
 * @endcode
 * which renders to
 * @code{.cpp}
 * void foo(T bar) {
 *   (void)(bar);
 *   // some logic not using parameter `bar`
 * }
 * @endcode
 * which is the standard and compiler independet way of omitting warnings on unused parameters while
 * still being able to fully document the parameter with Doxygen.
 *
 * @param[in] expr parameter to be denoted as unused
 * @since v0.2.0
 * @ingroup Utilities
 */
#define UNUSED(expr) \
  (void)(expr)


#endif  // _GLOBALS__HPP_
