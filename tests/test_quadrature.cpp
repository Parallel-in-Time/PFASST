/*
 * Tests for quadrature related routines: polynomials, nodes, and matrices.
 */

#include <iostream>
#include <tuple>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <pfasst/quadrature.hpp>

using namespace std;

TEST(PolyTest, LegendrePolys)
{
  auto l0 = pfasst::Polynomial<double>::legendre(0);
  EXPECT_EQ(l0.order(), 0);
  EXPECT_EQ(l0[0], 1.0);

  auto l1 = pfasst::Polynomial<double>::legendre(1);
  EXPECT_EQ(l1.order(), 1);
  EXPECT_EQ(l1[0], 0.0);
  EXPECT_EQ(l1[1], 1.0);

  auto l2 = pfasst::Polynomial<double>::legendre(2);
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

MATCHER(DoubleNear, "")
{
  return abs(get<0>(arg) - get<1>(arg)) < 1e-15;
}

TEST(NodesTest, GaussLegendreNodes)
{
  const long double l3e[3] = { 0.11270166537925831,
                               0.5,
                               0.8872983346207417
                             };

  const long double l5e[5] = { 0.046910077030668004,
                               0.23076534494715845,
                               0.5,
                               0.7692346550528415,
                               0.953089922969332
                             };

  const long double l7e[7] = { 0.025446043828620736,
                               0.12923440720030277,
                               0.2970774243113014,
                               0.5,
                               0.7029225756886985,
                               0.8707655927996972,
                               0.9745539561713793
                             };

  auto l3 = pfasst::compute_nodes<long double>(3, "gauss-legendre");
  EXPECT_THAT(l3, testing::Pointwise(DoubleNear(), l3e));

  auto l5 = pfasst::compute_nodes<long double>(5, "gauss-legendre");
  EXPECT_THAT(l5, testing::Pointwise(DoubleNear(), l5e));

  auto l7 = pfasst::compute_nodes<long double>(7, "gauss-legendre");
  EXPECT_THAT(l7, testing::Pointwise(DoubleNear(), l7e));
}

TEST(NodesTest, GaussLobattoNodes)
{
  const long double l2e[2] = { 0.0,
                               1.0
                             };

  const long double l3e[3] = { 0.0,
                               0.5,
                               1.0
                             };

  const long double l5e[5] = { 0.0,
                               0.17267316464601143,
                               0.5,
                               0.8273268353539885,
                               1.0
                             };

  const long double l7e[7] = { 0.0,
                               0.08488805186071653,
                               0.2655756032646429,
                               0.5,
                               0.7344243967353571,
                               0.9151119481392834,
                               1.0
                             };

  const long double l9e[9] = { 0.0,
                               0.05012100229426992,
                               0.16140686024463113,
                               0.3184412680869109,
                               0.5,
                               0.6815587319130891,
                               0.8385931397553689,
                               0.94987899770573,
                               1.0
                             };

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

TEST(NodesTest, ClenshawCurtisNodes)
{
  const long double cc2e[2] = { 0.0,
                                1.0
                              };

  const long double cc3e[3] = { 0.0,
                                0.5,
                                1.0
                              };

  const long double cc5e[5] = { 0.0,
                                0.14644660940672623779957781894757548,
                                0.5,
                                0.85355339059327376220042218105242452,
                                1.0
                              };

  const long double cc7e[7] = { 0.0,
                                0.066987298107780676618138414623531908,
                                0.25,
                                0.5,
                                0.75,
                                0.93301270189221932338186158537646809,
                                1.0
                              };

  const long double cc9e[9] = { 0.0,
                                0.038060233744356621935908405301605857,
                                0.14644660940672623779957781894757548,
                                0.30865828381745511413577000798480057,
                                0.5,
                                0.69134171618254488586422999201519943,
                                0.85355339059327376220042218105242452,
                                0.96193976625564337806409159469839414,
                                1.0
                              };

  auto cc2 = pfasst::compute_nodes<long double>(2, "clenshaw-curtis");
  EXPECT_THAT(cc2, testing::Pointwise(DoubleNear(), cc2e));

  auto cc3 = pfasst::compute_nodes<long double>(3, "clenshaw-curtis");
  EXPECT_THAT(cc3, testing::Pointwise(DoubleNear(), cc3e));

  auto cc5 = pfasst::compute_nodes<long double>(5, "clenshaw-curtis");
  EXPECT_THAT(cc5, testing::Pointwise(DoubleNear(), cc5e));

  auto cc7 = pfasst::compute_nodes<long double>(7, "clenshaw-curtis");
  EXPECT_THAT(cc7, testing::Pointwise(DoubleNear(), cc7e));

  auto cc9 = pfasst::compute_nodes<long double>(9, "clenshaw-curtis");
  EXPECT_THAT(cc9, testing::Pointwise(DoubleNear(), cc9e));
}

TEST(NodesTest, UniformNodes)
{
  const long double u2e[2] = { 0.0,
                               1.0
                             };

  const long double u3e[3] = { 0.0,
                               0.5,
                               1.0
                             };

  const long double u5e[5] = { 0.0,
                               0.25,
                               0.5,
                               0.75,
                               1.0
                             };

  auto u2 = pfasst::compute_nodes<long double>(2, "uniform");
  EXPECT_THAT(u2, testing::Pointwise(DoubleNear(), u2e));

  auto u3 = pfasst::compute_nodes<long double>(3, "uniform");
  EXPECT_THAT(u3, testing::Pointwise(DoubleNear(), u3e));

  auto u5 = pfasst::compute_nodes<long double>(5, "uniform");
  EXPECT_THAT(u5, testing::Pointwise(DoubleNear(), u5e));
}

TEST(QuadratureTest, GaussLobattoNodes)
{
  auto l3 = pfasst::compute_nodes<long double>(3, "gauss-lobatto");
  auto a3 = pfasst::augment_nodes(l3);
  auto s3 = pfasst::compute_quadrature(get<0>(a3), get<0>(a3), get<1>(a3), pfasst::QuadratureMatrix::S);
  const long double s3e[6] = { 0.20833333333333333,
                               0.33333333333333333,
                               -0.04166666666666666,
                               -0.04166666666666666,
                               0.33333333333333333,
                               0.20833333333333333
                             };
    
  EXPECT_THAT(s3.data(), testing::Pointwise(DoubleNear(), s3e));

  auto l5 = pfasst::compute_nodes<long double>(5, "gauss-lobatto");
  auto a5 = pfasst::augment_nodes(l5);
  auto s5 = pfasst::compute_quadrature(get<0>(a5), get<0>(a5), get<1>(a5), pfasst::QuadratureMatrix::S);
  const long double s5e[] = { 0.067728432186156897969267419174073482,
                              0.11974476934341168251615379970493965,
                              -0.021735721866558113665511351745074292,
                              0.010635824225415491883105056997129926,
                              -0.0037001392424145306021611522544979462,
                              -0.027103432186156897969267419174073483,
                              0.1834394139796310955018131986775051,
                              0.19951349964433589144328912952285207,
                              -0.041597785326236047678849833157352459,
                              0.013075139242414530602161152254497946,
                              0.013075139242414530602161152254497944,
                              -0.041597785326236047678849833157352467,
                              0.19951349964433589144328912952285207,
                              0.1834394139796310955018131986775051,
                              -0.027103432186156897969267419174073483,
                              -0.0037001392424145306021611522544979483,
                              0.010635824225415491883105056997129916,
                              -0.021735721866558113665511351745074289,
                              0.11974476934341168251615379970493965,
                              0.067728432186156897969267419174073482
                            };
  EXPECT_THAT(s5.data(), testing::Pointwise(DoubleNear(), s5e));
}

TEST(QuadratureTest, ClenshawCurtisNodes)
{
  auto c4 = pfasst::compute_nodes<long double>(4, "clenshaw-curtis");
  auto a4 = pfasst::augment_nodes(c4);
  auto s4 = pfasst::compute_quadrature(get<0>(a4), get<0>(a4), get<1>(a4), pfasst::QuadratureMatrix::S);
  const long double s4e[] = { 0.10243055555555555555555555555555556,
                              0.16319444444444444444444444444444444,
                              -0.024305555555555555555555555555555556,
                              0.0086805555555555555555555555555555557,
                              -0.055555555555555555555555555555555556,
                              0.30555555555555555555555555555555556,
                              0.30555555555555555555555555555555556,
                              -0.055555555555555555555555555555555556,
                              0.0086805555555555555555555555555555545,
                              -0.024305555555555555555555555555555554,
                              0.16319444444444444444444444444444444,
                              0.10243055555555555555555555555555556
                            };
  EXPECT_THAT(s4.data(), testing::Pointwise(DoubleNear(), s4e));
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
