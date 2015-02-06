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
  for (auto&& elem : vec) {
    elem = RandomGenerator::roll();
  }
}

#define create_single(name, num)\
  ParticleComponent<PRECISION> name(num)
#define create_and_fill_single(name, num)\
  create_single(name, num);\
  fill_single(name)
#define create_cloud(name, num)\
  create_single(name, num * DIMS)
#define create_and_fill_cloud(name, num)\
  create_and_fill_single(name, num * DIMS)


TEST(OperatorTests, AddSingleOnSingle)
{
  create_and_fill_single(first_single, DIMS);
  create_and_fill_single(second_single, DIMS);
  create_single(expected_single, DIMS);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] + second_single[i];
  }
  ParticleComponent<PRECISION> result_single = first_single + second_single;
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, AddCloudOnCloud)
{
  create_and_fill_cloud(first_cloud, 5);
  create_and_fill_cloud(second_cloud, 5);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = first_cloud[i * DIMS + j] + second_cloud[i * DIMS + j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = first_cloud + second_cloud;
  EXPECT_THAT(result_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, AddSingleOnCloud)
{
  create_and_fill_cloud(cloud, 5);
  create_and_fill_single(single, DIMS);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] + single[j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = cloud + single;
  EXPECT_THAT(result_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, InplaceAddSingleOnSingle)
{
  create_and_fill_single(first_single, DIMS);
  create_and_fill_single(second_single, DIMS);
  create_single(expected_single, DIMS);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] + second_single[i];
  }
  first_single += second_single;
  EXPECT_THAT(first_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceAddCloudOnCloud)
{
  create_and_fill_cloud(first_cloud, 5);
  create_and_fill_cloud(second_cloud, 5);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = first_cloud[i * DIMS + j] + second_cloud[i * DIMS + j];
    }
  }
  first_cloud += second_cloud;
  EXPECT_THAT(first_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, InplaceAddSingleOnCloud)
{
  create_and_fill_cloud(cloud, 5);
  create_and_fill_single(single, DIMS);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] + single[j];
    }
  }
  cloud += single;
  EXPECT_THAT(cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, MinusSingleOnSingle)
{
  create_and_fill_single(first_single, DIMS);
  create_and_fill_single(second_single, DIMS);
  create_single(expected_single, DIMS);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] - second_single[i];
  }
  ParticleComponent<PRECISION> result_single = first_single - second_single;
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, MinusCloudOnCloud)
{
  create_and_fill_cloud(first_cloud, 5);
  create_and_fill_cloud(second_cloud, 5);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = first_cloud[i * DIMS + j] - second_cloud[i * DIMS + j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = first_cloud - second_cloud;
  EXPECT_THAT(result_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, MinusSingleOnCloud)
{
  create_and_fill_cloud(cloud, 5);
  create_and_fill_single(single, DIMS);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] - single[j];
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = cloud - single;
  EXPECT_THAT(result_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, InplaceMinusSingleOnSingle)
{
  create_and_fill_single(first_single, DIMS);
  create_and_fill_single(second_single, DIMS);
  create_single(expected_single, DIMS);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = first_single[i] - second_single[i];
  }
  first_single -= second_single;
  EXPECT_THAT(first_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceMinusCloudOnCloud)
{
  create_and_fill_cloud(first_cloud, 5);
  create_and_fill_cloud(second_cloud, 5);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = first_cloud[i * DIMS + j] - second_cloud[i * DIMS + j];
    }
  }
  first_cloud -= second_cloud;
  EXPECT_THAT(first_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, InplaceMinusSingleOnCloud)
{
  create_and_fill_cloud(cloud, 5);
  create_and_fill_single(single, DIMS);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] - single[j];
    }
  }
  cloud -= single;
  EXPECT_THAT(cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, MulWithSingle)
{
  create_and_fill_single(single, DIMS);
  PRECISION value = 2.0;
  create_single(expected_single, DIMS);
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
  create_and_fill_cloud(cloud, 5);
  PRECISION value = 2.0;
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] * value;
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud1 = cloud * value;
  EXPECT_THAT(result_cloud1, Pointwise(Eq(), expected_cloud));
  ParticleCloudComponent<PRECISION> result_cloud2 = value * cloud;
  EXPECT_THAT(result_cloud2, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, InplaceMulWithSingle)
{
  create_and_fill_single(single, DIMS);
  PRECISION value = 2.0;
  create_single(expected_single, DIMS);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = single[i] * value;
  }
  single *= value;
  EXPECT_THAT(single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceMulWithCloud)
{
  create_and_fill_cloud(cloud, 5);
  PRECISION value = 2.0;
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] * value;
    }
  }
  cloud *= value;
  EXPECT_THAT(cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, DivWithSingle)
{
  create_and_fill_single(single, DIMS);
  PRECISION value = 2.0;
  create_single(expected_single, DIMS);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = single[i] / value;
  }
  ParticleComponent<PRECISION> result_single = single / value;
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, DivWithCloud)
{
  create_and_fill_cloud(cloud, 5);
  PRECISION value = 2.0;
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] / value;
    }
  }
  ParticleCloudComponent<PRECISION> result_cloud = cloud / value;
  EXPECT_THAT(result_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, InplaceDivWithSingle)
{
  create_and_fill_single(single, DIMS);
  PRECISION value = 2.0;
  create_single(expected_single, DIMS);
  for(size_t i = 0; i < DIMS; ++i) {
    expected_single[i] = single[i] / value;
  }
  single /= value;
  EXPECT_THAT(single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, InplaceDivWithCloud)
{
  create_and_fill_cloud(cloud, 5);
  PRECISION value = 2.0;
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < DIMS; ++j) {
      expected_cloud[i * DIMS + j] = cloud[i * DIMS + j] / value;
    }
  }
  cloud /= value;
  EXPECT_THAT(cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, CrossProdSingleOnSingle)
{
  create_and_fill_single(first_single, DIMS);
  create_and_fill_single(second_single, DIMS);
  create_single(expected_single, DIMS);
  expected_single[0] = first_single[1] * second_single[2] - first_single[2] * second_single[1];
  expected_single[1] = first_single[2] * second_single[0] - first_single[0] * second_single[2];
  expected_single[2] = first_single[0] * second_single[1] - first_single[1] * second_single[0];

  ParticleComponent<PRECISION> result_single = cross_prod(first_single, second_single);
  EXPECT_THAT(result_single, Pointwise(Eq(), expected_single));
}

TEST(OperatorTests, CrossProdCloudOnCloud)
{
  create_and_fill_cloud(first_cloud, 5);
  create_and_fill_cloud(second_cloud, 5);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    expected_cloud[i * DIMS]     = first_cloud[i * DIMS + 1] * second_cloud[i * DIMS + 2] - first_cloud[i * DIMS + 2] * second_cloud[i * DIMS + 1];
    expected_cloud[i * DIMS + 1] = first_cloud[i * DIMS + 2] * second_cloud[i * DIMS]     - first_cloud[i * DIMS]     * second_cloud[i * DIMS + 2];
    expected_cloud[i * DIMS + 2] = first_cloud[i * DIMS]     * second_cloud[i * DIMS + 1] - first_cloud[i * DIMS + 1] * second_cloud[i * DIMS];
  }
  ParticleCloudComponent<PRECISION> result_cloud = cross_prod(first_cloud, second_cloud);
    EXPECT_THAT(result_cloud, Pointwise(Eq(), expected_cloud));
}

TEST(OperatorTests, CrossProdSingleOnCloud)
{
  create_and_fill_cloud(cloud, 5);
  create_and_fill_single(single, DIMS);
  create_cloud(expected_cloud, 5);
  for (size_t i = 0; i < 5; ++i) {
    expected_cloud[i * DIMS]     = cloud[i * DIMS + 1] * single[2] - cloud[i * DIMS + 2] * single[1];
    expected_cloud[i * DIMS + 1] = cloud[i * DIMS + 2] * single[0] - cloud[i * DIMS]     * single[2];
    expected_cloud[i * DIMS + 2] = cloud[i * DIMS]     * single[1] - cloud[i * DIMS + 1] * single[0];
  }
  ParticleCloudComponent<PRECISION> result_cloud = cross_prod(cloud, single);
  EXPECT_THAT(result_cloud, Pointwise(Eq(), expected_cloud));
}


TEST(UtilitiesTests, DistanceBetweenTwoParticles)
{
  Particle<PRECISION> first_particle(DIMS);
  create_and_fill_single(first_single, DIMS);
  first_particle.pos() = first_single;
  Particle<PRECISION> second_particle(DIMS);
  create_and_fill_single(second_single, DIMS);
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
  create_and_fill_single(single, DIMS);
  particle.pos() = single;
  ParticleCloud<PRECISION> cloud(5, DIMS);
  create_and_fill_cloud(temp, 5);
  cloud.positions() = temp;
  vector<PRECISION> expected_distance(5, 0.0);
  for (size_t j = 0; j < 5; ++j) {
    for (size_t i = 0; i < DIMS; ++i) {
      expected_distance[j] += (cloud.positions()[j * DIMS + i] - particle.pos()[i]) * (cloud.positions()[j * DIMS + i] - particle.pos()[i]);
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
