/*
 * Tests for the scalar example solving the test equation
 */ 
 
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#define PFASST_UNIT_TESTING
#include "../examples/scalar/scalar_sdc.cpp"
#undef PFASST_UNIT_TESTING

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}