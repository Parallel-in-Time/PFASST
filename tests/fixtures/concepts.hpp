#ifndef _TESTS__FIXTURES__CONCEPTS_HPP_
#define _TESTS__FIXTURES__CONCEPTS_HPP_

#include "fixtures/test_helpers.hpp"

#include <type_traits>
using namespace std;


template<class T>
class Concepts
  : public ::testing::Test
{};

TYPED_TEST_CASE_P(Concepts);

TYPED_TEST_P(Concepts, default_constructible)
{
  EXPECT_TRUE(is_default_constructible<TypeParam>::value);
  EXPECT_TRUE(is_destructible<TypeParam>::value);
}

TYPED_TEST_P(Concepts, move_and_copy_construtible)
{
  EXPECT_TRUE(is_copy_constructible<TypeParam>::value);
  EXPECT_TRUE(is_move_constructible<TypeParam>::value);
}

TYPED_TEST_P(Concepts, assignable)
{
  const bool assignable = is_assignable<TypeParam, TypeParam>::value;
  EXPECT_TRUE(assignable);
}

TYPED_TEST_P(Concepts, move_and_copy_assignable)
{
  EXPECT_TRUE(is_copy_assignable<TypeParam>::value);
  EXPECT_TRUE(is_move_assignable<TypeParam>::value);
}

REGISTER_TYPED_TEST_CASE_P(Concepts,
                           default_constructible, move_and_copy_construtible,
                           assignable, move_and_copy_assignable);

#endif // _TESTS__FIXTURES__CONCEPTS_HPP_
