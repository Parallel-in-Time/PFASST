#include <memory>
#include <random>
#include <type_traits>
using namespace std;

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include "../examples/boris/particle.hpp"
#include "../examples/boris/particle_cloud.hpp"
#include "../examples/boris/particle_util.hpp"
using namespace pfasst::examples::boris;

#define DIMS 3
#define PRECISION double

struct RandomGenerator
{
  static default_random_engine rd_gen;
  static uniform_real_distribution<PRECISION> dist;

  static PRECISION roll()
  {
    return dist(rd_gen);
  }
};

default_random_engine RandomGenerator::rd_gen = default_random_engine(42);
uniform_real_distribution<PRECISION> RandomGenerator::dist = uniform_real_distribution<PRECISION>(-10.0, 10.0);


void fill_single(ParticleComponent<PRECISION>& vec)
{
  assert(vec.size() == DIMS);
  for (size_t i = 0; i < DIMS; ++i) {
    vec[i] = RandomGenerator::roll();
  }
}

#define create_single(name)\
  ParticleComponent<PRECISION> name(DIMS)
#define create_and_fill_single(name)\
  create_single(name);\
  fill_single(name);

void fill_cloud(ParticleCloudComponent<PRECISION>& cloud)
{
  create_single(temp);
  for (size_t i = 0; i < cloud.size(); ++i) {
    fill_single(temp);
    cloud[i] = temp;
  }
}

#define create_cloud(num, name)\
  ParticleCloudComponent<PRECISION> name;\
  for(size_t i=0; i < num; ++i) name.push_back(ParticleComponent<PRECISION>(DIMS));
#define create_and_fill_cloud(num, name)\
  create_cloud(num, name);\
  fill_cloud(name);


TEST(OperatorTests, AddSingleOnSingle)
{
  create_and_fill_single(first_single);
  create_and_fill_single(second_single);
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] + second_single[i];
  }
  ParticleComponent<PRECISION> result_single = first_single + second_single;
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, AddCloudOnCloud)
{
  create_and_fill_cloud(5, first_cloud);
  create_and_fill_cloud(5, second_cloud);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = first_cloud[i][j] + second_cloud[i][j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = first_cloud + second_cloud;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, AddSingleOnCloud)
{
  create_and_fill_cloud(5, cloud);
  create_and_fill_single(single);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] + single[j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = cloud + single;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, InplaceAddSingleOnSingle)
{
  create_and_fill_single(first_single);
  create_and_fill_single(second_single);
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] + second_single[i];
  }
  first_single += second_single;
  EXPECT_THAT(first_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceAddCloudOnCloud)
{
  create_and_fill_cloud(5, first_cloud);
  create_and_fill_cloud(5, second_cloud);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = first_cloud[i][j] + second_cloud[i][j];
    }
  }
  first_cloud += second_cloud;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(first_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, InplaceAddSingleOnCloud)
{
  create_and_fill_cloud(5, cloud);
  create_and_fill_single(single);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] + single[j];
    }
  }
  cloud += single;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, MinusSingleOnSingle)
{
  create_and_fill_single(first_single);
  create_and_fill_single(second_single);
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] - second_single[i];
  }
  ParticleComponent<PRECISION> result_single = first_single - second_single;
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, MinusCloudOnCloud)
{
  create_and_fill_cloud(5, first_cloud);
  create_and_fill_cloud(5, second_cloud);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = first_cloud[i][j] - second_cloud[i][j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = first_cloud - second_cloud;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, MinusSingleOnCloud)
{
  create_and_fill_cloud(5, cloud);
  create_and_fill_single(single);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] - single[j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = cloud - single;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, InplaceMinusSingleOnSingle)
{
  create_and_fill_single(first_single);
  create_and_fill_single(second_single);
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] - second_single[i];
  }
  first_single -= second_single;
  EXPECT_THAT(first_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceMinusCloudOnCloud)
{
  create_and_fill_cloud(5, first_cloud);
  create_and_fill_cloud(5, second_cloud);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = first_cloud[i][j] - second_cloud[i][j];
    }
  }
  first_cloud -= second_cloud;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(first_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, InplaceMinusSingleOnCloud)
{
  create_and_fill_cloud(5, cloud);
  create_and_fill_single(single);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] - single[j];
    }
  }
  cloud -= single;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, MulWithSingle)
{
  create_and_fill_single(single);
  PRECISION value = 2.0;
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = single[i] * value;
  }
  ParticleComponent<PRECISION> result_single1 = single * value;
  EXPECT_THAT(result_single1, Pointwise(Eq(), expected_single));
  ParticleComponent<PRECISION> result_single2 = value * single;
  EXPECT_THAT(result_single2, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, MulWithCloud)
{
  create_and_fill_cloud(5, cloud);
  PRECISION value = 2.0;
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] * value;
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud1 = cloud * value;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud1[i], Pointwise(Eq(), expected_cloud[i]));
  }
  ParticleCloudComponent<PRECISION> result_cloud2 = value * cloud;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud2[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, InplaceMulWithSingle)
{
  create_and_fill_single(single);
  PRECISION value = 2.0;
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = single[i] * value;
  }
  single *= value;
  EXPECT_THAT(single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceMulWithCloud)
{
  create_and_fill_cloud(5, cloud);
  PRECISION value = 2.0;
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] * value;
    }
  }
  cloud *= value;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, DivWithSingle)
{
  create_and_fill_single(single);
  PRECISION value = 2.0;
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = single[i] / value;
  }
  ParticleComponent<PRECISION> result_single = single / value;
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, DivWithCloud)
{
  create_and_fill_cloud(5, cloud);
  PRECISION value = 2.0;
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] / value;
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = cloud / value;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, InplaceDivWithSingle)
{
  create_and_fill_single(single);
  PRECISION value = 2.0;
  create_single(expected_single);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = single[i] / value;
  }
  single /= value;
  EXPECT_THAT(single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceDivWithCloud)
{
  create_and_fill_cloud(5, cloud);
  PRECISION value = 2.0;
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i][j] = cloud[i][j] / value;
    }
  }
  cloud /= value;
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, CrossProdSingleOnSingle)
{
  create_and_fill_single(first_single);
  create_and_fill_single(second_single);
  create_single(expected_single);
  expected_single[0] = first_single[1] * second_single[2] - first_single[2] * second_single[1];
  expected_single[1] = first_single[2] * second_single[0] - first_single[0] * second_single[2];
  expected_single[2] = first_single[0] * second_single[1] - first_single[1] * second_single[0];

  ParticleComponent<PRECISION> result_single = cross_prod(first_single, second_single);
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, CrossProdCloudOnCloud)
{
  create_and_fill_cloud(5, first_cloud);
  create_and_fill_cloud(5, second_cloud);
  create_cloud(5, expected_cloud)
  for (size_t i = 0; i < 5; ++i) {
    expected_cloud[i][0] = first_cloud[i][1] * second_cloud[i][2] - first_cloud[i][2] * second_cloud[i][1];
    expected_cloud[i][1] = first_cloud[i][2] * second_cloud[i][0] - first_cloud[i][0] * second_cloud[i][2];
    expected_cloud[i][2] = first_cloud[i][0] * second_cloud[i][1] - first_cloud[i][1] * second_cloud[i][0];
  }
  ParticleCloudComponent<PRECISION> result_cloud = cross_prod(first_cloud, second_cloud);
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}

TEST(OperatorTests, CrossProdSingleOnCloud)
{
  create_and_fill_cloud(5, cloud);
  create_and_fill_single(single);
  create_cloud(5, expected_cloud);
  for (size_t i = 0; i < 5; ++i) {
    expected_cloud[i][0] = cloud[i][1] * single[2] - cloud[i][2] * single[1];
    expected_cloud[i][1] = cloud[i][2] * single[0] - cloud[i][0] * single[2];
    expected_cloud[i][2] = cloud[i][0] * single[1] - cloud[i][1] * single[0];
  }
  ParticleCloudComponent<PRECISION> result_cloud = cross_prod(cloud, single);
  for(size_t i = 0; i < 5; ++i) {
    EXPECT_THAT(result_cloud[i], Pointwise(Eq(), expected_cloud[i]));
  }
}


TEST(UtilitiesTests, DistanceBetweenTwoParticles)
{
  Particle<PRECISION> first_particle(DIMS);
  create_and_fill_single(first_single);
  first_particle.pos() = first_single;
  Particle<PRECISION> second_particle(DIMS);
  create_and_fill_single(second_single);
  second_particle.pos() = second_single;
  PRECISION expected_distance = 0.0;
  for (size_t i = 0; i < DIMS; ++i) {
    expected_distance += (first_single[i] - second_single[i]) * (first_single[i] - second_single[i]);
  }
  expected_distance = sqrt(expected_distance);
  PRECISION dist = distance(first_particle, second_particle);
  EXPECT_THAT(dist, DoubleEq(expected_distance));
}

TEST(UtilitiesTests, DistanceOfCloudToReference)
{
  Particle<PRECISION> particle(DIMS);
  create_and_fill_single(single);
  particle.pos() = single;
  ParticleCloud<PRECISION> cloud(5, DIMS);
  create_single(temp);
  for (size_t i = 0; i < cloud.size(); ++i) {
    fill_single(temp);
    cloud.positions()[i] = temp;
  }
  vector<PRECISION> expected_distance(5, 0.0);
  for (size_t j = 0; j < 5; ++j) {
    for (size_t i = 0; i < DIMS; ++i) {
      expected_distance[j] += (cloud.positions()[j][i] - particle.pos()[i]) * (cloud.positions()[j][i] - particle.pos()[i]);
    }
    expected_distance[j] = sqrt(expected_distance[j]);
  }
  vector<PRECISION> dist = distance_to_reference(cloud, particle);
  EXPECT_THAT(dist, Pointwise(Eq(), expected_distance));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
