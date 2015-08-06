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
 *
 * @since v0.2.0
 *
 * @ingroup Utilities
 */
#define UNUSED(expr) \
  (void)(expr)

//! @{
//! @ingroup Internals

#if defined(__GNUC__)
  //! @see [Stackoverflow Q&A](http://stackoverflow.com/a/8990275/588243) for reference
  #define DEPRECATE(foo, msg) foo __attribute__((deprecated(msg)))
#elif defined(_MSC_VER)
  //! @see [Stackoverflow Q&A](http://stackoverflow.com/a/8990275/588243) for reference
  #define DEPRECATE(foo, msg) __declspec(deprecated(msg)) foo
#else
  #error This compiler is not supported
#endif

//! @see [Stackoverflow Q&A](http://stackoverflow.com/a/8990275/588243) for reference
#define PP_CAT(x,y) PP_CAT1(x,y)

//! @see [Stackoverflow Q&A](http://stackoverflow.com/a/8990275/588243) for reference
 *   for reference
#define PP_CAT1(x,y) x##y

namespace detail
{
    struct true_type {};
    struct false_type {};
    template <int test> struct converter : public true_type {};
    template <> struct converter<0> : public false_type {};
}
//! @}

/**
 * Produces compiler warnings when static condition @p cond is not met at compile time.
 *
 * This is similar to `static_assert` with the variation that evaluation of @p cond to `std::false_type` will not
 * yield compile errors but deprecation warnings.
 *
 * @param[in] cond
 * @param[in] msg
 *
 * @since v0.6.0
 * @see [Stackoverflow Q&A](http://stackoverflow.com/a/8990275/588243)
 *   for reference
 *
 * @ingroup Utilities
 */
#define STATIC_WARNING(cond, msg) \
struct PP_CAT(static_warning,__LINE__) { \
  DEPRECATE(void _(::detail::false_type const& ),msg) {}; \
  void _(::detail::true_type const& ) {}; \
  PP_CAT(static_warning,__LINE__)() {_(::detail::converter<(cond)>());} \
}

/**
 * @param[in] token must be a program-wide unique identifier.
 *
 * @note Using `STATIC_WARNING_TEMPLATE` changes the meaning of a program in a small way.
 *  It introduces a member/variable declaration.
 *  This means at least one byte of space in each structure/class instantiation.
 *  @ref STATIC_WARNING should be preferred in any non-template situation.
 * @since v0.6.0
 * @see [Stackoverflow Q&A](http://stackoverflow.com/a/8990275/588243)
 *   for reference
 *
 * @ingroup Internals
 */
#define STATIC_WARNING_TEMPLATE(token, cond, msg) \
    STATIC_WARNING(cond, msg) PP_CAT(PP_CAT(_localvar_, token),__LINE__)


namespace pfasst
{
  typedef double time_precision;
}

#endif  // _GLOBALS__HPP_
