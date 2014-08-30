#ifndef _GLOBALS__HPP_
#define _GLOBALS__HPP_

/**
 * denoting unused function parameters for omitting compiler warnings
 *
 * To denote an unused function parameter just use this makro on it in the function body:
 *
 * ~~~{.cpp}
 * void func(auto p) {
 *   UNUSED(p);
 * }
 * ~~~
 *
 * which renders to
 *
 * ~~~{.cpp}
 * void func(auto p) {
 *   (void)(p);
 * }
 * ~~~
 *
 * which is the standard and compiler independet way of omitting warnings on unused parameters while
 * still being able to fully document the parameter with Doxygen.
 *
 * @param[in] expr parameter to be denoted as unused
 */
#define UNUSED(expr) (void)(expr)

#endif  // _GLOBALS__HPP_
