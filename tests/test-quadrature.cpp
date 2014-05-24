/*
 * Tests for quadrature related routines: polynomials, nodes, and matrices.
 */

#include <iostream>
#include <tuple>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <pfasst-quadrature.hpp>

using namespace std;

TEST(PolyTest, LegendrePolys) {
  auto l0 = pfasst::polynomial<double>::legendre(0);
  EXPECT_EQ(l0.order(), 0);
  EXPECT_EQ(l0[0], 1.0);

  auto l1 = pfasst::polynomial<double>::legendre(1);
  EXPECT_EQ(l1.order(), 1);
  EXPECT_EQ(l1[0], 0.0);
  EXPECT_EQ(l1[1], 1.0);

  auto l2 = pfasst::polynomial<double>::legendre(2);
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

MATCHER(DoubleNear, "") {
  return abs(get<0>(arg) - get<1>(arg)) < 1e-15;
}

TEST(NodesTest, GaussLegendreNodes) {
  const long double l3e[3] = { 0.11270166537925831,
			       0.5,
			       0.8872983346207417 };

  const long double l5e[5] = { 0.046910077030668004,
			       0.23076534494715845,
			       0.5,
			       0.7692346550528415,
			       0.953089922969332 };

  const long double l7e[7] = { 0.025446043828620736,
			       0.12923440720030277,
			       0.2970774243113014,
			       0.5,
			       0.7029225756886985,
			       0.8707655927996972,
			       0.9745539561713793 };

  auto l3 = pfasst::compute_nodes<long double>(3, "gauss-legendre");
  EXPECT_THAT(l3, testing::Pointwise(DoubleNear(), l3e));

  auto l5 = pfasst::compute_nodes<long double>(5, "gauss-legendre");
  EXPECT_THAT(l5, testing::Pointwise(DoubleNear(), l5e));

  auto l7 = pfasst::compute_nodes<long double>(7, "gauss-legendre");
  EXPECT_THAT(l7, testing::Pointwise(DoubleNear(), l7e));
}

TEST(NodesTest, GaussLobattoNodes) {
  const long double l2e[2] = { 0.0,
			       1.0 };

  const long double l3e[3] = { 0.0,
			       0.5,
			       1.0 };

  const long double l5e[5] = { 0.0,
			       0.17267316464601143,
			       0.5,
			       0.8273268353539885,
			       1.0 };

  const long double l7e[7] = { 0.0,
			       0.08488805186071653,
			       0.2655756032646429,
			       0.5,
			       0.7344243967353571,
			       0.9151119481392834,
			       1.0 };

  const long double l9e[9] = { 0.0,
			       0.05012100229426992,
			       0.16140686024463113,
			       0.3184412680869109,
			       0.5,
			       0.6815587319130891,
			       0.8385931397553689,
			       0.94987899770573,
			       1.0 };

  auto l2 = pfasst::compute_nodes<long double>(2, "gauss-lobatto");
  EXPECT_THAT(l2, testing::Pointwise(DoubleNear(), l2e));

  auto l3 = pfasst::compute_nodes<long double>(3, "gauss-lobatto");
  EXPECT_THAT(l3, testing::Pointwise(DoubleNear(), l3e));

  auto l5 = pfasst::compute_nodes<long double>(5, "gauss-lobatto");
  EXPECT_THAT(l5, testing::Pointwise(DoubleNear(), l5e));

  auto l7 = pfasst::compute_nodes<long double>(7, "gauss-lobatto");
  EXPECT_THAT(l7, testing::Pointwise(DoubleNear(), l7e));

  auto l9 = pfasst::compute_nodes<long double>(9, "gauss-lobatto");
  EXPECT_THAT(l9, testing::Pointwise(DoubleNear(), l9e));
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
