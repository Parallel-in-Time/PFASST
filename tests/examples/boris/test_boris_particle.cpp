#include <memory>
#include <type_traits>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace ::testing;

#include "../examples/boris/particle.hpp"
#include "../examples/boris/particle_3d.hpp"
using namespace pfasst::examples::boris;

typedef ParticleEncapsulation<double, double> ParticleT;

typedef Position3DEncapsulation<double, double> Position3D;
typedef Velocity3DEncapsulation<double, double> Velocity3D;
typedef Acceleration3DEncapsulation<double, double> Acceleration3D;
typedef Particle3DEncapsulation<double, double> Particle3DT;


/*
 * ParticleComponent3DEncapsulation
 */

template<typename T>
class OperatorTests :
  public ::testing::Test
{
  public:
    const double c1, c2, c3;
    const T val1;
    T val2;
    T val1_copy;

  public:
    OperatorTests()
      :   c1(0.1), c2(1.2), c3(-42.0)
        , val1( c1,  c2,  c3)
        , val2(-c1, -c2, -c3)
    {}
    virtual ~OperatorTests()
    {}
};

TYPED_TEST_CASE_P(OperatorTests);

TYPED_TEST_P(OperatorTests, Initialization)
{
  EXPECT_TRUE(is_default_constructible<TypeParam>::value);
  EXPECT_TRUE(is_destructible<TypeParam>::value);
  EXPECT_TRUE(is_copy_constructible<TypeParam>::value);
  EXPECT_TRUE(is_move_constructible<TypeParam>::value);

  TypeParam default_val;
  EXPECT_THAT(default_val.DIM, 3);
  EXPECT_THAT(default_val[0], DoubleEq(0.0));
  EXPECT_THAT(default_val[1], DoubleEq(0.0));
  EXPECT_THAT(default_val[2], DoubleEq(0.0));

  EXPECT_THAT(this->val1[0], DoubleEq(this->c1));
  EXPECT_THAT(this->val1[1], DoubleEq(this->c2));
  EXPECT_THAT(this->val1[2], DoubleEq(this->c3));
}

TYPED_TEST_P(OperatorTests, Copyable)
{
  EXPECT_TRUE(is_copy_assignable<TypeParam>::value);
  EXPECT_TRUE(is_move_assignable<TypeParam>::value);

  TypeParam copy1(this->val1);
  EXPECT_THAT(copy1[0], DoubleEq(this->val1[0]));
  EXPECT_THAT(copy1[1], DoubleEq(this->val1[1]));
  EXPECT_THAT(copy1[2], DoubleEq(this->val1[2]));

  TypeParam copy2;
  copy2 = this->val1;
  EXPECT_THAT(copy2[0], DoubleEq(this->val1[0]));
  EXPECT_THAT(copy2[1], DoubleEq(this->val1[1]));
  EXPECT_THAT(copy2[2], DoubleEq(this->val1[2]));
}

TYPED_TEST_P(OperatorTests, axpy)
{
  this->val2.saxpy(1.0, this->val1);
  EXPECT_THAT(this->val2[0], DoubleEq(0.0));
  EXPECT_THAT(this->val2[1], DoubleEq(0.0));
  EXPECT_THAT(this->val2[2], DoubleEq(0.0));

  shared_ptr<TypeParam> spr1 = make_shared<TypeParam>(this->c1, this->c2, this->c3);
  spr1->saxpy(0.0, this->val1);
  EXPECT_THAT(spr1->operator[](0), DoubleEq(this->c1));
  EXPECT_THAT(spr1->operator[](1), DoubleEq(this->c2));
  EXPECT_THAT(spr1->operator[](2), DoubleEq(this->c3));

  shared_ptr<TypeParam> spr2 = make_shared<TypeParam>(this->c1, this->c2, this->c3);
  spr2->saxpy(-1.0, this->val1);
  EXPECT_THAT(spr2->operator[](0), DoubleEq(0.0));
  EXPECT_THAT(spr2->operator[](1), DoubleEq(0.0));
  EXPECT_THAT(spr2->operator[](2), DoubleEq(0.0));
}


TYPED_TEST_P(OperatorTests, Addition)
{
  TypeParam add1 = this->val1 + this->val2;
  EXPECT_THAT(add1[0], DoubleEq(this->c1 + (- this->c1)));
  EXPECT_THAT(add1[1], DoubleEq(this->c2 + (- this->c2)));
  EXPECT_THAT(add1[2], DoubleEq(this->c3 + (- this->c3)));
  EXPECT_THAT(this->val2[0], DoubleEq(- this->c1));
  EXPECT_THAT(this->val2[1], DoubleEq(- this->c2));
  EXPECT_THAT(this->val2[2], DoubleEq(- this->c3));

  TypeParam add2 = this->val1 + 1.0;
  EXPECT_THAT(add2[0], DoubleEq(this->c1 + 1.0));
  EXPECT_THAT(add2[1], DoubleEq(this->c2 + 1.0));
  EXPECT_THAT(add2[2], DoubleEq(this->c3 + 1.0));

  TypeParam add3 = 1.0 + this->val1;
  EXPECT_THAT(add3[0], DoubleEq(this->c1 + 1.0));
  EXPECT_THAT(add3[1], DoubleEq(this->c2 + 1.0));
  EXPECT_THAT(add3[2], DoubleEq(this->c3 + 1.0));

  this->val1_copy = this->val1;
  this->val1_copy += 1.0;
  EXPECT_THAT(this->val1_copy[0], DoubleEq(this->c1 + 1.0));
  EXPECT_THAT(this->val1_copy[1], DoubleEq(this->c2 + 1.0));
  EXPECT_THAT(this->val1_copy[2], DoubleEq(this->c3 + 1.0));

  this->val1_copy = this->val1;
  this->val1_copy += this->val2;
  EXPECT_THAT(this->val1_copy[0], DoubleEq(this->c1 + (- this->c1)));
  EXPECT_THAT(this->val1_copy[1], DoubleEq(this->c2 + (- this->c2)));
  EXPECT_THAT(this->val1_copy[2], DoubleEq(this->c3 + (- this->c3)));
}

TYPED_TEST_P(OperatorTests, Subtraction)
{
  TypeParam sub1 = this->val1 - this->val2;
  EXPECT_THAT(sub1[0], DoubleEq(this->c1 - (- this->c1)));
  EXPECT_THAT(sub1[1], DoubleEq(this->c2 - (- this->c2)));
  EXPECT_THAT(sub1[2], DoubleEq(this->c3 - (- this->c3)));
  EXPECT_THAT(this->val2[0], DoubleEq(- this->c1));
  EXPECT_THAT(this->val2[1], DoubleEq(- this->c2));
  EXPECT_THAT(this->val2[2], DoubleEq(- this->c3));

  TypeParam sub2 = this->val1 - 1.0;
  EXPECT_THAT(sub2[0], DoubleEq(this->c1 - 1.0));
  EXPECT_THAT(sub2[1], DoubleEq(this->c2 - 1.0));
  EXPECT_THAT(sub2[2], DoubleEq(this->c3 - 1.0));

  TypeParam sub3 = 1.0 - this->val1;
  EXPECT_THAT(sub3[0], DoubleEq(this->c1 - 1.0));
  EXPECT_THAT(sub3[1], DoubleEq(this->c2 - 1.0));
  EXPECT_THAT(sub3[2], DoubleEq(this->c3 - 1.0));

  this->val1_copy = this->val1;
  this->val1_copy -= 1.0;
  EXPECT_THAT(this->val1_copy[0], DoubleEq(this->c1 - 1.0));
  EXPECT_THAT(this->val1_copy[1], DoubleEq(this->c2 - 1.0));
  EXPECT_THAT(this->val1_copy[2], DoubleEq(this->c3 - 1.0));

  this->val1_copy = this->val1;
  this->val1_copy -= this->val2;
  EXPECT_THAT(this->val1_copy[0], DoubleEq(this->c1 - (- this->c1)));
  EXPECT_THAT(this->val1_copy[1], DoubleEq(this->c2 - (- this->c2)));
  EXPECT_THAT(this->val1_copy[2], DoubleEq(this->c3 - (- this->c3)));
}

TYPED_TEST_P(OperatorTests, Multiplication)
{
  TypeParam twice = this->val1 * 2.0;
  EXPECT_THAT(twice[0], DoubleEq(2.0 * this->c1));
  EXPECT_THAT(twice[1], DoubleEq(2.0 * this->c2));
  EXPECT_THAT(twice[2], DoubleEq(2.0 * this->c3));

  TypeParam twice2 = 2.0 * this->val1;
  EXPECT_THAT(twice2[0], DoubleEq(2.0 * this->c1));
  EXPECT_THAT(twice2[1], DoubleEq(2.0 * this->c2));
  EXPECT_THAT(twice2[2], DoubleEq(2.0 * this->c3));
}

TYPED_TEST_P(OperatorTests, Division)
{
  TypeParam twice = this->val1 / 2.0;
  EXPECT_THAT(twice[0], DoubleEq(this->c1 / 2.0));
  EXPECT_THAT(twice[1], DoubleEq(this->c2 / 2.0));
  EXPECT_THAT(twice[2], DoubleEq(this->c3 / 2.0));
}

REGISTER_TYPED_TEST_CASE_P(OperatorTests,
                           Initialization, Copyable, axpy,
                           Addition, Subtraction, Multiplication, Division);

typedef ::testing::Types<Position3D, Velocity3D, Acceleration3D> Dim3Types;
INSTANTIATE_TYPED_TEST_CASE_P(ParticleComponent3D, OperatorTests, Dim3Types);

/*
 * Velocity3DEncapsulation
 */
TEST(BorisParticleVelocity3D, convert)
{
  Velocity3D vel1(0.1, 0.2, -0.1);

  Position3D converted = vel1.convert(dt<double>(2.0));
  EXPECT_THAT(converted.x, DoubleEq(0.2));
  EXPECT_THAT(converted.y, DoubleEq(0.4));
  EXPECT_THAT(converted.z, DoubleEq(-0.2));
}

/*
 * Acceleration3DEncapsulation
 */
TEST(BorisParticleAcceleration3D, convert)
{
  Acceleration3D accel1(0.1, 0.2, -0.5);

  Position3D converted_pos = accel1.convert(dtdt<double>(2.0));
  EXPECT_THAT(converted_pos.x, DoubleEq(0.2));
  EXPECT_THAT(converted_pos.y, DoubleEq(0.4));

  Velocity3D converted_vel = accel1.convert(dt<double>(2.0));
  EXPECT_THAT(converted_vel.u, DoubleEq(0.2));
  EXPECT_THAT(converted_vel.v, DoubleEq(0.4));
}


/*
 * Particle3DEncapsulation
 */
TEST(BorisParticle3D, Instantiation)
{
  EXPECT_TRUE(is_default_constructible<Particle3DT>::value);
  EXPECT_TRUE(is_copy_constructible<Particle3DT>::value);
  EXPECT_TRUE(is_move_constructible<Particle3DT>::value);
  EXPECT_TRUE(is_destructible<Particle3DT>::value);

  Particle3DT default_ctor;
  EXPECT_THAT(default_ctor.DIM, 3);
  EXPECT_THAT(default_ctor.pos().x, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.pos().y, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel().u, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel().v, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel().a, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel().b, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.mass(), DoubleEq(1.0));
  EXPECT_THAT(default_ctor.charge(), DoubleEq(1.0));

  Particle3DT special_ctor = Particle3DT(0.5, 1.0);
  EXPECT_THAT(special_ctor.pos().x, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.pos().y, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel().u, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel().v, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel().a, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel().b, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.mass(), DoubleEq(0.5));
  EXPECT_THAT(special_ctor.charge(), DoubleEq(1.0));
}

TEST(BorisParticle3D, Copyable)
{
  Particle3DT original = Particle3DT(0.1, 0.2);
  original.pos().x = 0.3;
  EXPECT_THAT(original.pos().x, DoubleEq(0.3));
  original.pos().y = 0.4;
  EXPECT_THAT(original.pos().y, DoubleEq(0.4));

  original.vel().u = 0.5;
  EXPECT_THAT(original.vel().u, DoubleEq(0.5));
  original.vel().v = 0.6;
  EXPECT_THAT(original.vel().v, DoubleEq(0.6));

  original.accel().a = 0.7;
  EXPECT_THAT(original.accel().a, DoubleEq(0.7));
  original.accel().b = 0.8;
  EXPECT_THAT(original.accel().b, DoubleEq(0.8));

  Particle3DT expl_copy;
  expl_copy = original;
  EXPECT_THAT(expl_copy.pos().x, DoubleEq(0.3));
  EXPECT_THAT(expl_copy.pos().y, DoubleEq(0.4));
  EXPECT_THAT(expl_copy.vel().u, DoubleEq(0.5));
  EXPECT_THAT(expl_copy.vel().v, DoubleEq(0.6));
  EXPECT_THAT(expl_copy.accel().a, DoubleEq(0.7));
  EXPECT_THAT(expl_copy.accel().b, DoubleEq(0.8));
  EXPECT_THAT(expl_copy.mass(), DoubleEq(0.1));
  EXPECT_THAT(expl_copy.charge(), DoubleEq(0.2));

//   Particle3DT original = Particle3DT(0.1, 0.2);
  original.pos().x = 0.3;
  original.pos().y = 0.4;
  original.vel().u = 0.5;
  original.vel().v = 0.6;
  original.accel().a = 0.7;
  original.accel().b = 0.8;
  Particle3DT impl_copy = original;
  EXPECT_THAT(impl_copy.pos().x, DoubleEq(0.3));
  EXPECT_THAT(impl_copy.pos().y, DoubleEq(0.4));
  EXPECT_THAT(impl_copy.vel().u, DoubleEq(0.5));
  EXPECT_THAT(impl_copy.vel().v, DoubleEq(0.6));
  EXPECT_THAT(impl_copy.accel().a, DoubleEq(0.7));
  EXPECT_THAT(impl_copy.accel().b, DoubleEq(0.8));
  EXPECT_THAT(impl_copy.mass(), DoubleEq(0.1));
  EXPECT_THAT(impl_copy.charge(), DoubleEq(0.2));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
