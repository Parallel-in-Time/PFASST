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

/*
 * the following is taken from SO:
 * http://stackoverflow.com/a/8990275/588243
 */
#if defined(__GNUC__)
#define DEPRECATE(foo, msg) foo __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
#define DEPRECATE(foo, msg) __declspec(deprecated(msg)) foo
#else
#error This compiler is not supported
#endif

#define PP_CAT(x,y) PP_CAT1(x,y)
#define PP_CAT1(x,y) x##y

namespace detail
{
    struct true_type {};
    struct false_type {};
    template <int test> struct converter : public true_type {};
    template <> struct converter<0> : public false_type {};
}

#define STATIC_WARNING(cond, msg) \
struct PP_CAT(static_warning,__LINE__) { \
  DEPRECATE(void _(::detail::false_type const& ),msg) {}; \
  void _(::detail::true_type const& ) {}; \
  PP_CAT(static_warning,__LINE__)() {_(::detail::converter<(cond)>());} \
}

// Note: using STATIC_WARNING_TEMPLATE changes the meaning of a program in a small way.
// It introduces a member/variable declaration.  This means at least one byte of space
// in each structure/class instantiation.  STATIC_WARNING should be preferred in any 
// non-template situation.
//  'token' must be a program-wide unique identifier.
#define STATIC_WARNING_TEMPLATE(token, cond, msg) \
    STATIC_WARNING(cond, msg) PP_CAT(PP_CAT(_localvar_, token),__LINE__)


namespace pfasst
{
  typedef double time_precision;
}

#endif  // _GLOBALS__HPP_
