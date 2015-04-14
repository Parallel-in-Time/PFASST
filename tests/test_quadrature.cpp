/*
 * Tests for quadrature related routines: nodes and matrices.
 */
#include <iostream>
#include <algorithm>
#include <iostream>
#include <tuple>
#include <type_traits>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <Eigen/Dense>

using namespace ::testing;

#include <pfasst/quadrature.hpp>

#include "fixtures/concepts.hpp"

using namespace std;
using namespace pfasst::quadrature;

MATCHER(DoubleNear, "")
{
  return abs(get<0>(arg) - get<1>(arg)) < 1e-15;
}

TEST(NodesTest, GaussLegendreNodes)
{
  EXPECT_FALSE(GaussLegendre<>(3).left_is_node());
  EXPECT_FALSE(GaussLegendre<>(3).right_is_node());

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

  auto l3 = compute_nodes<long double>(3, QuadratureType::GaussLegendre);
  EXPECT_THAT(l3, Pointwise(DoubleNear(), l3e));

  auto l5 = compute_nodes<long double>(5, QuadratureType::GaussLegendre);
  EXPECT_THAT(l5, Pointwise(DoubleNear(), l5e));

  auto l7 = compute_nodes<long double>(7, QuadratureType::GaussLegendre);
  EXPECT_THAT(l7, Pointwise(DoubleNear(), l7e));
}

TEST(NodesTest, GaussLobattoNodes)
{
  EXPECT_TRUE(GaussLobatto<>(3).left_is_node());
  EXPECT_TRUE(GaussLobatto<>(3).right_is_node());

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

  auto l2 = compute_nodes<long double>(2, QuadratureType::GaussLobatto);
  EXPECT_THAT(l2, Pointwise(DoubleNear(), l2e));

  auto l3 = compute_nodes<long double>(3, QuadratureType::GaussLobatto);
  EXPECT_THAT(l3, Pointwise(DoubleNear(), l3e));

  auto l5 = compute_nodes<long double>(5, QuadratureType::GaussLobatto);
  EXPECT_THAT(l5, Pointwise(DoubleNear(), l5e));

  auto l7 = compute_nodes<long double>(7, QuadratureType::GaussLobatto);
  EXPECT_THAT(l7, Pointwise(DoubleNear(), l7e));

  auto l9 = compute_nodes<long double>(9, QuadratureType::GaussLobatto);
  EXPECT_THAT(l9, Pointwise(DoubleNear(), l9e));
}

TEST(NodesTest, ClenshawCurtisNodes)
{
  EXPECT_TRUE(ClenshawCurtis<>(3).left_is_node());
  EXPECT_TRUE(ClenshawCurtis<>(3).right_is_node());

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

  auto cc2 = compute_nodes<long double>(2, QuadratureType::ClenshawCurtis);
  EXPECT_THAT(cc2, Pointwise(DoubleNear(), cc2e));

  auto cc3 = compute_nodes<long double>(3, QuadratureType::ClenshawCurtis);
  EXPECT_THAT(cc3, Pointwise(DoubleNear(), cc3e));

  auto cc5 = compute_nodes<long double>(5, QuadratureType::ClenshawCurtis);
  EXPECT_THAT(cc5, Pointwise(DoubleNear(), cc5e));

  auto cc7 = compute_nodes<long double>(7, QuadratureType::ClenshawCurtis);
  EXPECT_THAT(cc7, Pointwise(DoubleNear(), cc7e));

  auto cc9 = compute_nodes<long double>(9, QuadratureType::ClenshawCurtis);
  EXPECT_THAT(cc9, Pointwise(DoubleNear(), cc9e));
}

TEST(NodesTest, UniformNodes)
{
  EXPECT_TRUE(Uniform<>(3).left_is_node());
  EXPECT_TRUE(Uniform<>(3).right_is_node());

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

  auto u2 = compute_nodes<long double>(2, QuadratureType::Uniform);
  EXPECT_THAT(u2, Pointwise(DoubleNear(), u2e));

  auto u3 = compute_nodes<long double>(3, QuadratureType::Uniform);
  EXPECT_THAT(u3, Pointwise(DoubleNear(), u3e));

  auto u5 = compute_nodes<long double>(5, QuadratureType::Uniform);
  EXPECT_THAT(u5, Pointwise(DoubleNear(), u5e));
}


TEST(QuadratureTest, GaussLobattoNodes)
{
  Eigen::Matrix<double, 3, 3> s_mat_3 = GaussLobatto<double>(3).get_s_mat();
  Eigen::Matrix<double, 3, 3> s_mat_3_expected;
  s_mat_3_expected <<  0.0,                 0.0,                  0.0,
                       0.20833333333333333, 0.33333333333333333, -0.04166666666666666,
                      -0.04166666666666666, 0.33333333333333333, 0.20833333333333333;
  for (size_t row = 0; row < 3; ++row) {
    for (size_t col = 0; col < 3; ++col) {
      EXPECT_THAT(s_mat_3(row, col), ::testing::DoubleNear((s_mat_3_expected(row, col)), (double)(1e-14)));
    }
  }

  Eigen::Matrix<double, 5, 5> s_mat_5 = GaussLobatto<double>(5).get_s_mat();
  Eigen::Matrix<double, 5, 5> s_mat_5_expected;
  s_mat_5_expected <<  0.0, 0.0, 0.0, 0.0, 0.0,
                       0.067728432186156897969267419174073482, 0.11974476934341168251615379970493965,
                          -0.021735721866558113665511351745074292, 0.010635824225415491883105056997129926,
                          -0.0037001392424145306021611522544979462,
                      -0.027103432186156897969267419174073483, 0.1834394139796310955018131986775051,
                           0.19951349964433589144328912952285207, -0.041597785326236047678849833157352459,
                           0.013075139242414530602161152254497946,
                       0.013075139242414530602161152254497944, -0.041597785326236047678849833157352467,
                           0.19951349964433589144328912952285207, 0.1834394139796310955018131986775051,
                          -0.027103432186156897969267419174073483,
                      -0.0037001392424145306021611522544979483, 0.010635824225415491883105056997129916,
                          -0.021735721866558113665511351745074289, 0.11974476934341168251615379970493965,
                           0.067728432186156897969267419174073482;

  for (size_t row = 0; row < 5; ++row) {
    for (size_t col = 0; col < 5; ++col) {
      EXPECT_THAT(s_mat_5(row, col), ::testing::DoubleNear((s_mat_5_expected(row, col)), (double)(1e-14)));
    }
  }

  Eigen::Matrix<double, 5, 5> q_mat_5 = GaussLobatto<double>(5).get_q_mat();
  vector<double> q_vec = GaussLobatto<double>(5).get_q_vec();

  for (size_t col = 0; col < 5; ++col) {
    EXPECT_THAT(q_vec[col], ::testing::DoubleEq(q_mat_5(4, col)));
  }
}

TEST(QuadratureTest, ClenshawCurtisNodes)
{
  Eigen::Matrix<double, 4, 4> s_mat_4 = ClenshawCurtis<double>(4).get_s_mat();
  Eigen::Matrix<double, 4, 4> s_mat_4_expected;
  s_mat_4_expected <<  0.0, 0.0, 0.0, 0.0,
                       0.10243055555555555555555555555555556,  0.16319444444444444444444444444444444,
                         -0.024305555555555555555555555555555556,  0.0086805555555555555555555555555555557,
                      -0.055555555555555555555555555555555556,  0.30555555555555555555555555555555556,
                          0.30555555555555555555555555555555556, -0.055555555555555555555555555555555556,
                       0.0086805555555555555555555555555555545, -0.024305555555555555555555555555555554,
                          0.16319444444444444444444444444444444,  0.10243055555555555555555555555555556;

  for (size_t row = 0; row < 4; ++row) {
    for (size_t col = 0; col < 4; ++col) {
      EXPECT_THAT(s_mat_4(row, col), ::testing::DoubleNear((s_mat_4_expected(row, col)), (double)(1e-14)));
    }
  }
}


class QmatTest
  : public ::TestWithParam<tuple<size_t, QuadratureType>>
{
  protected:
    size_t nnodes;
    QuadratureType qtype;
    shared_ptr<IQuadrature<long double>> quad;

  public:
    virtual void SetUp()
    {
      nnodes = get<0>(GetParam());
      qtype = get<1>(GetParam());

      this->quad = quadrature_factory<long double>(nnodes, qtype);

      cout << "Quadrature type no. " << int(qtype) << " -- Number of nodes " << nnodes << endl;
    }
};

TEST_P(QmatTest, AllNodes)
{
  for (int m = 0; m < this->quad->get_q_mat().rows(); ++m) {
    long double qsum = 0;
    for (int j = 0; j < this->quad->get_q_mat().cols(); ++j) {
      qsum += this->quad->get_q_mat()(m,j);
    }
    EXPECT_NEAR(qsum, this->quad->get_nodes()[m], (long double)(3E-12));
  }
}


INSTANTIATE_TEST_CASE_P(Quadrature, QmatTest,
                        ::Combine(::Range<size_t>(2, 14),
                                  Values<QuadratureType>(QuadratureType::GaussLegendre,
                                                         QuadratureType::GaussLobatto,
                                                         QuadratureType::GaussRadau,
                                                         QuadratureType::ClenshawCurtis,
                                                         QuadratureType::Uniform)));


typedef ::testing::Types<GaussLegendre<>,
                         GaussLobatto<>,
                         GaussRadau<>,
                         ClenshawCurtis<>,
                         Uniform<>> QuadratureTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Quadrature, ConceptsTest, QuadratureTypes);


int main(int argc, char** argv)
{
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
