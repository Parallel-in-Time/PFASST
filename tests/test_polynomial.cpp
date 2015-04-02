/*
 * Tests for Polynomials
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>

using namespace ::testing;

#include <pfasst/quadrature/polynomial.hpp>


TEST(PolynomialTest, Legendre)
{
  auto l0 = pfasst::quadrature::Polynomial<double>::legendre(0);
  EXPECT_EQ(l0.order(), 0);
  EXPECT_EQ(l0[0], 1.0);

  auto l1 = pfasst::quadrature::Polynomial<double>::legendre(1);
  EXPECT_EQ(l1.order(), 1);
  EXPECT_EQ(l1[0], 0.0);
  EXPECT_EQ(l1[1], 1.0);

  auto l2 = pfasst::quadrature::Polynomial<double>::legendre(2);
  EXPECT_EQ(l2.order(), 2);
  EXPECT_EQ(l2[0], -0.5);
  EXPECT_EQ(l2[1],  0.0);
  EXPECT_EQ(l2[2],  1.5);

  auto l2d = l2.differentiate();
  EXPECT_EQ(l2d.order(), 1);
  EXPECT_EQ(l2d[0], 0.0);
  EXPECT_EQ(l2d[1], 3.0);

  auto l2i = l2.integrate();
  EXPECT_EQ(l2i.order(), 3);
  EXPECT_EQ(l2i[0], 0.0);
  EXPECT_EQ(l2i[1], -0.5);
  EXPECT_EQ(l2i[2], 0.0);
  EXPECT_EQ(l2i[3], 0.5);

  double a1 = l2.evaluate(1.0);
  EXPECT_EQ(a1, 1.0);
}


int main(int argc, char** argv)
{
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
