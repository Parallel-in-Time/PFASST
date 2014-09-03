#ifndef _TESTS__FIXTURES__CONCEPTS_HPP_
#define _TESTS__FIXTURES__CONCEPTS_HPP_

#include <type_traits>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

template<typename T>
class ConceptsTest
  : public ::testing::Test
{};

TYPED_TEST_CASE_P(ConceptsTest);

TYPED_TEST_P(ConceptsTest, Constructors)
{
  EXPECT_TRUE(is_default_constructible<TypeParam>::value);
  EXPECT_TRUE(is_destructible<TypeParam>::value);
  EXPECT_TRUE(is_copy_constructible<TypeParam>::value);
  EXPECT_TRUE(is_move_constructible<TypeParam>::value);
}

TYPED_TEST_P(ConceptsTest, Assignments)
{
  EXPECT_TRUE(is_copy_assignable<TypeParam>::value);
  EXPECT_TRUE(is_move_assignable<TypeParam>::value);
}

REGISTER_TYPED_TEST_CASE_P(ConceptsTest,
                           Constructors, Assignments);

#endif // _TESTS__FIXTURES__CONCEPTS_HPP_
