#include <memory>
#include <type_traits>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../examples/boris/particle.hpp"
#include "../examples/boris/particle_3d.hpp"

using namespace ::testing;

typedef ParticleEncapsulation<double, double> ParticleT;

typedef Position3DEncapsulation<double, double> Position3D;
typedef Velocity3DEncapsulation<double, double> Velocity3D;
typedef Acceleration3DEncapsulation<double, double> Acceleration3D;
typedef Particle3DEncapsulation<double, double> Particle3DT;


TEST(BorisParticle, Instantiation)
{
  ParticleT default_ctor;
  EXPECT_THAT(default_ctor.DIM, 0);
  EXPECT_THAT(default_ctor.charge, DoubleEq(1.0));
  EXPECT_THAT(default_ctor.mass, DoubleEq(1.0));

  ParticleT special_ctor = ParticleT(0.5, 1.0);
  EXPECT_THAT(special_ctor.mass, DoubleEq(0.5));
  EXPECT_THAT(special_ctor.charge, DoubleEq(1.0));
}

/*
 * Position3DEncapsulation
 */
TEST(BorisParticlePosition3D, Initialization)
{
  EXPECT_TRUE(is_default_constructible<Position3D>::value);
  EXPECT_TRUE(is_destructible<Position3D>::value);
  EXPECT_TRUE(is_copy_constructible<Position3D>::value);
  EXPECT_TRUE(is_move_constructible<Position3D>::value);

  Position3D default_pos;
  EXPECT_THAT(default_pos.DIM, 3);
  EXPECT_THAT(default_pos.x, DoubleEq(0.0));
  EXPECT_THAT(default_pos.y, DoubleEq(0.0));
  EXPECT_THAT(default_pos.z, DoubleEq(0.0));

  shared_ptr<const Position3D> pos2 = make_shared<const Position3D>(0.5, 1.0, -0.5);
  EXPECT_THAT(pos2->x, DoubleEq(0.5));
  EXPECT_THAT(pos2->y, DoubleEq(1.0));
  EXPECT_THAT(pos2->z, DoubleEq(-0.5));
}

TEST(BorisParticlePosition3D, Copyable)
{
  EXPECT_TRUE(is_copy_assignable<Position3D>::value);
  EXPECT_TRUE(is_move_assignable<Position3D>::value);

  shared_ptr<const Position3D> pos_original = make_shared<const Position3D>(0.5, 1.0, -0.5);

  shared_ptr<Position3D> pos_copy1 = make_shared<Position3D>();
  pos_copy1->copy(pos_original);
  EXPECT_THAT(pos_copy1->x, DoubleEq(0.5));
  EXPECT_THAT(pos_copy1->y, DoubleEq(1.0));

  Position3D pos_copy2;
  pos_copy2 = pos_original;
  EXPECT_THAT(pos_copy2.x, DoubleEq(0.5));
  EXPECT_THAT(pos_copy2.y, DoubleEq(1.0));

  Position3D pos_copy3;
  pos_copy3 = *(pos_original.get());
  EXPECT_THAT(pos_copy3.x, DoubleEq(0.5));
  EXPECT_THAT(pos_copy3.y, DoubleEq(1.0));
}

TEST(BorisParticlePosition3D, axpy)
{
  shared_ptr<Position3D> pos_2 = make_shared<Position3D>(0.5, 1.0, -0.5);
  shared_ptr<const Position3D> other_pos = make_shared<const Position3D>(0.5, 1.0, -0.5);

  pos_2->saxpy(1.0, other_pos);
  EXPECT_THAT(pos_2->x, DoubleEq(1.0));
  EXPECT_THAT(pos_2->y, DoubleEq(2.0));

  shared_ptr<Position3D> pos_1 = make_shared<Position3D>(0.5, 1.0, -0.5);
  pos_1->saxpy(0.0, other_pos);
  EXPECT_THAT(pos_1->x, DoubleEq(0.5));
  EXPECT_THAT(pos_1->y, DoubleEq(1.0));

  shared_ptr<Position3D> pos_0 = make_shared<Position3D>(0.5, 1.0, -0.5);
  pos_0->saxpy(-1.0, other_pos);
  EXPECT_THAT(pos_0->x, DoubleEq(0.0));
  EXPECT_THAT(pos_0->y, DoubleEq(0.0));
}

TEST(BorisParticlePosition3D, operators)
{
  Position3D pos1(0.1, 0.2, -0.1);
  Position3D pos2(-0.1, -0.2, 0.1);

  long double one = 1.0;

  Position3D twice = pos1 * 2.0;
  EXPECT_THAT(twice.x, DoubleEq(0.2));
  EXPECT_THAT(twice.y, DoubleEq(0.4));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));

  Position3D twice2 = 2.0 * pos1;
  EXPECT_THAT(twice2.x, DoubleEq(0.2));
  EXPECT_THAT(twice2.y, DoubleEq(0.4));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));

  Position3D add1 = pos1 + pos2;
  EXPECT_THAT(add1.x, DoubleEq(0.0));
  EXPECT_THAT(add1.y, DoubleEq(0.0));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));
  EXPECT_THAT(pos2.x, DoubleEq(-0.1));
  EXPECT_THAT(pos2.y, DoubleEq(-0.2));

  Position3D add2 = pos1 + one;
  EXPECT_THAT(add2.x, DoubleEq(1.1));
  EXPECT_THAT(add2.y, DoubleEq(1.2));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));

  Position3D add3 = one + pos1;
  EXPECT_THAT(add3.x, DoubleEq(1.1));
  EXPECT_THAT(add3.y, DoubleEq(1.2));
  EXPECT_THAT(pos1.x, DoubleEq(0.1));
  EXPECT_THAT(pos1.y, DoubleEq(0.2));
}


/*
 * Velocity3DEncapsulation
 */
TEST(BorisParticleVelocity3D, Initialization)
{
  EXPECT_TRUE(is_default_constructible<Velocity3D>::value);
  EXPECT_TRUE(is_copy_constructible<Velocity3D>::value);
  EXPECT_TRUE(is_move_constructible<Velocity3D>::value);
  EXPECT_TRUE(is_destructible<Velocity3D>::value);

  Velocity3D default_vel;
  EXPECT_THAT(default_vel.DIM, 3);
  EXPECT_THAT(default_vel.u, DoubleEq(0.0));
  EXPECT_THAT(default_vel.v, DoubleEq(0.0));

  shared_ptr<const Velocity3D> vel2 = make_shared<const Velocity3D>(0.5, 1.0, -0.5);
  EXPECT_THAT(vel2->u, DoubleEq(0.5));
  EXPECT_THAT(vel2->v, DoubleEq(1.0));
}

TEST(BorisParticleVelocity3D, Copyable)
{
  EXPECT_TRUE(is_copy_assignable<Velocity3D>::value);
  EXPECT_TRUE(is_move_assignable<Velocity3D>::value);

  shared_ptr<Velocity3D> vel_copy = make_shared<Velocity3D>();
  shared_ptr<const Velocity3D> vel2 = make_shared<const Velocity3D>(0.5, 1.0, -0.5);
  vel_copy->copy(vel2);
  EXPECT_THAT(vel_copy->u, DoubleEq(0.5));
  EXPECT_THAT(vel_copy->v, DoubleEq(1.0));
}

TEST(BorisParticleVelocity3D, axpy)
{
  shared_ptr<Velocity3D> vel_2 = make_shared<Velocity3D>(0.5, 1.0, -0.5);
  shared_ptr<const Velocity3D> vel2 = make_shared<const Velocity3D>(0.5, 1.0, -0.5);
  vel_2->saxpy(1.0, vel2);
  EXPECT_THAT(vel_2->u, DoubleEq(1.0));
  EXPECT_THAT(vel_2->v, DoubleEq(2.0));

  shared_ptr<Velocity3D> vel_1 = make_shared<Velocity3D>(0.5, 1.0, -0.5);
  vel_1->saxpy(0.0, vel2);
  EXPECT_THAT(vel_1->u, DoubleEq(0.5));
  EXPECT_THAT(vel_1->v, DoubleEq(1.0));

  shared_ptr<Velocity3D> vel_0 = make_shared<Velocity3D>(0.5, 1.0, -0.5);
  vel_0->saxpy(-1.0, vel2);
  EXPECT_THAT(vel_0->u, DoubleEq(0.0));
  EXPECT_THAT(vel_0->v, DoubleEq(0.0));
}

TEST(BorisParticleVelocity3D, convert)
{
  Velocity3D vel1(0.1, 0.2, -0.1);

  Position3D converted = vel1.convert(dt<double>(2.0));
  EXPECT_THAT(converted.x, DoubleEq(0.2));
  EXPECT_THAT(converted.y, DoubleEq(0.4));
}

TEST(BorisParticleVelocity3D, operators)
{
  Velocity3D vel1(0.1, 0.2, -0.1);
  Velocity3D vel2(-0.1, -0.2, 0.1);

  long double one = 1.0;

  Velocity3D twice = vel1 * 2.0;
  EXPECT_THAT(twice.u, DoubleEq(0.2));
  EXPECT_THAT(twice.v, DoubleEq(0.4));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));

  Velocity3D twice2 = 2.0 * vel1;
  EXPECT_THAT(twice2.u, DoubleEq(0.2));
  EXPECT_THAT(twice2.v, DoubleEq(0.4));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));

  Velocity3D add1 = vel1 + vel2;
  EXPECT_THAT(add1.u, DoubleEq(0.0));
  EXPECT_THAT(add1.v, DoubleEq(0.0));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));
  EXPECT_THAT(vel2.u, DoubleEq(-0.1));
  EXPECT_THAT(vel2.v, DoubleEq(-0.2));

  Velocity3D add2 = vel1 + one;
  EXPECT_THAT(add2.u, DoubleEq(1.1));
  EXPECT_THAT(add2.v, DoubleEq(1.2));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));

  Velocity3D add3 = one + vel1;
  EXPECT_THAT(add3.u, DoubleEq(1.1));
  EXPECT_THAT(add3.v, DoubleEq(1.2));
  EXPECT_THAT(vel1.u, DoubleEq(0.1));
  EXPECT_THAT(vel1.v, DoubleEq(0.2));
}


/*
 * Acceleration3DEncapsulation
 */
TEST(BorisParticleAcceleration3D, Initialization)
{
  EXPECT_TRUE(is_default_constructible<Acceleration3D>::value);
  EXPECT_TRUE(is_copy_constructible<Acceleration3D>::value);
  EXPECT_TRUE(is_move_constructible<Acceleration3D>::value);
  EXPECT_TRUE(is_destructible<Acceleration3D>::value);

  Acceleration3D default_accel;
  EXPECT_THAT(default_accel.DIM, 3);
  EXPECT_THAT(default_accel.a, DoubleEq(0.0));
  EXPECT_THAT(default_accel.b, DoubleEq(0.0));

  shared_ptr<const Acceleration3D> accel2 = make_shared<const Acceleration3D>(0.5, 1.0, -0.5);
  EXPECT_THAT(accel2->a, DoubleEq(0.5));
  EXPECT_THAT(accel2->b, DoubleEq(1.0));
}

TEST(BorisParticleAcceleration3D, Copyable)
{
  shared_ptr<Acceleration3D> accel_copy = make_shared<Acceleration3D>();
  shared_ptr<const Acceleration3D> accel2 = make_shared<const Acceleration3D>(0.5, 1.0, -0.5);
  accel_copy->copy(accel2);
  EXPECT_THAT(accel_copy->a, DoubleEq(0.5));
  EXPECT_THAT(accel_copy->b, DoubleEq(1.0));
}

TEST(BorisParticleAcceleration3D, axpy)
{
  shared_ptr<Acceleration3D> accel_2 = make_shared<Acceleration3D>(0.5, 1.0, -0.5);
  shared_ptr<const Acceleration3D> accel2 = make_shared<const Acceleration3D>(0.5, 1.0, -0.5);
  accel_2->saxpy(1.0, accel2);
  EXPECT_THAT(accel_2->a, DoubleEq(1.0));
  EXPECT_THAT(accel_2->b, DoubleEq(2.0));

  shared_ptr<Acceleration3D> accel_1 = make_shared<Acceleration3D>(0.5, 1.0, -0.5);
  accel_1->saxpy(0.0, accel2);
  EXPECT_THAT(accel_1->a, DoubleEq(0.5));
  EXPECT_THAT(accel_1->b, DoubleEq(1.0));

  shared_ptr<Acceleration3D> accel_0 = make_shared<Acceleration3D>(0.5, 1.0, -0.5);
  accel_0->saxpy(-1.0, accel2);
  EXPECT_THAT(accel_0->a, DoubleEq(0.0));
  EXPECT_THAT(accel_0->b, DoubleEq(0.0));
}

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

TEST(BorisParticleAcceleration3D, operators)
{
  Acceleration3D accel1(0.1, 0.2, -0.1);
  Acceleration3D accel2(-0.1, -0.2, 0.1);

  long double one = 1.0;

  Acceleration3D twice = accel1 * 2.0;
  EXPECT_THAT(twice.a, DoubleEq(0.2));
  EXPECT_THAT(twice.b, DoubleEq(0.4));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));

  Acceleration3D twice2 = 2.0 * accel1;
  EXPECT_THAT(twice2.a, DoubleEq(0.2));
  EXPECT_THAT(twice2.b, DoubleEq(0.4));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));

  Acceleration3D add1 = accel1 + accel2;
  EXPECT_THAT(add1.a, DoubleEq(0.0));
  EXPECT_THAT(add1.b, DoubleEq(0.0));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));
  EXPECT_THAT(accel2.a, DoubleEq(-0.1));
  EXPECT_THAT(accel2.b, DoubleEq(-0.2));

  Acceleration3D add2 = accel1 + one;
  EXPECT_THAT(add2.a, DoubleEq(1.1));
  EXPECT_THAT(add2.b, DoubleEq(1.2));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));

  Acceleration3D add3 = one + accel1;
  EXPECT_THAT(add3.a, DoubleEq(1.1));
  EXPECT_THAT(add3.b, DoubleEq(1.2));
  EXPECT_THAT(accel1.a, DoubleEq(0.1));
  EXPECT_THAT(accel1.b, DoubleEq(0.2));
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
  EXPECT_THAT(default_ctor.pos->x, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.pos->y, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel->u, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel->v, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel->a, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel->b, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.mass, DoubleEq(1.0));
  EXPECT_THAT(default_ctor.charge, DoubleEq(1.0));

  Particle3DT special_ctor = Particle3DT(0.5, 1.0);
  EXPECT_THAT(special_ctor.pos->x, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.pos->y, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel->u, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel->v, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel->a, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel->b, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.mass, DoubleEq(0.5));
  EXPECT_THAT(special_ctor.charge, DoubleEq(1.0));
}

TEST(BorisParticle3D, Copyable)
{
  shared_ptr<Particle3DT> original_ptr = make_shared<Particle3DT>(0.1, 0.2);
  original_ptr->pos->x = 0.3;
  EXPECT_THAT(original_ptr->pos->x, DoubleEq(0.3));
  original_ptr->pos->y = 0.4;
  EXPECT_THAT(original_ptr->pos->y, DoubleEq(0.4));

  original_ptr->vel->u = 0.5;
  EXPECT_THAT(original_ptr->vel->u, DoubleEq(0.5));
  original_ptr->vel->v = 0.6;
  EXPECT_THAT(original_ptr->vel->v, DoubleEq(0.6));

  original_ptr->accel->a = 0.7;
  EXPECT_THAT(original_ptr->accel->a, DoubleEq(0.7));
  original_ptr->accel->b = 0.8;
  EXPECT_THAT(original_ptr->accel->b, DoubleEq(0.8));
  Particle3DT expl_copy;

  expl_copy.copy(const_pointer_cast<const Particle3DT>(original_ptr));
  EXPECT_THAT(expl_copy.pos->x, DoubleEq(0.3));
  EXPECT_THAT(expl_copy.pos->y, DoubleEq(0.4));
  EXPECT_THAT(expl_copy.vel->u, DoubleEq(0.5));
  EXPECT_THAT(expl_copy.vel->v, DoubleEq(0.6));
  EXPECT_THAT(expl_copy.accel->a, DoubleEq(0.7));
  EXPECT_THAT(expl_copy.accel->b, DoubleEq(0.8));
  EXPECT_THAT(expl_copy.mass, DoubleEq(0.1));
  EXPECT_THAT(expl_copy.charge, DoubleEq(0.2));

  Particle3DT original = Particle3DT(0.1, 0.2);
  original.pos->x = 0.3;
  original.pos->y = 0.4;
  original.vel->u = 0.5;
  original.vel->v = 0.6;
  original.accel->a = 0.7;
  original.accel->b = 0.8;
  Particle3DT impl_copy = original;
  EXPECT_THAT(impl_copy.pos->x, DoubleEq(0.3));
  EXPECT_THAT(impl_copy.pos->y, DoubleEq(0.4));
  EXPECT_THAT(impl_copy.vel->u, DoubleEq(0.5));
  EXPECT_THAT(impl_copy.vel->v, DoubleEq(0.6));
  EXPECT_THAT(impl_copy.accel->a, DoubleEq(0.7));
  EXPECT_THAT(impl_copy.accel->b, DoubleEq(0.8));
  EXPECT_THAT(impl_copy.mass, DoubleEq(0.1));
  EXPECT_THAT(impl_copy.charge, DoubleEq(0.2));
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
