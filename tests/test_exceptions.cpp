#include "fixtures/test_helpers.hpp"

#include <pfasst/exceptions.hpp>


TEST(ExceptionsTest, NotImplementedYet) {
  auto except = pfasst::NotImplementedYet("method");

  EXPECT_THAT(string(except.what()),
              StrEq("Not implemented/supported yet, required for: method"));
}


TEST(ExceptionsTest, ValueError) {
  auto except = pfasst::ValueError("wrong value");

  EXPECT_THAT(string(except.what()),
              StrEq("Value Error: wrong value"));
}


TEST_MAIN()
