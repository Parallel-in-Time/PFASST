#include <memory>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../examples/boris/particle.hpp"
#include "../examples/boris/particle_2d.hpp"

using namespace ::testing;

typedef ParticleEncapsulation<double, double> ParticleT;

typedef Position2DEncapsulation<double, double> Position2D;
typedef Velocity2DEncapsulation<double, double> Velocity2D;
typedef Acceleration2DEncapsulation<double, double> Acceleration2D;
typedef Particle2DEncapsulation<double, double> Particle2DT;


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


TEST(BorisParticleComponents2D, Position)
{
  Position2D default_pos;
  EXPECT_THAT(default_pos.DIM, 2);
  EXPECT_THAT(default_pos.x, DoubleEq(0.0));
  EXPECT_THAT(default_pos.y, DoubleEq(0.0));

  shared_ptr<const Position2D> pos2 = make_shared<const Position2D>(0.5, 1.0);
  EXPECT_THAT(pos2->x, DoubleEq(0.5));
  EXPECT_THAT(pos2->y, DoubleEq(1.0));

  shared_ptr<Position2D> pos_copy = make_shared<Position2D>();
  pos_copy->copy(pos2);
  EXPECT_THAT(pos_copy->x, DoubleEq(0.5));
  EXPECT_THAT(pos_copy->y, DoubleEq(1.0));
}

TEST(BorisParticleComponents2D, Velocity)
{
  Velocity2D default_vel;
  EXPECT_THAT(default_vel.DIM, 2);
  EXPECT_THAT(default_vel.u, DoubleEq(0.0));
  EXPECT_THAT(default_vel.v, DoubleEq(0.0));

  shared_ptr<const Velocity2D> vel2 = make_shared<const Velocity2D>(0.5, 1.0);
  EXPECT_THAT(vel2->u, DoubleEq(0.5));
  EXPECT_THAT(vel2->v, DoubleEq(1.0));

  shared_ptr<Velocity2D> vel_copy = make_shared<Velocity2D>();
  vel_copy->copy(vel2);
  EXPECT_THAT(vel_copy->u, DoubleEq(0.5));
  EXPECT_THAT(vel_copy->v, DoubleEq(1.0));
}

TEST(BorisParticleComponents2D, Acceleration)
{
  Acceleration2D default_accel;
  EXPECT_THAT(default_accel.DIM, 2);
  EXPECT_THAT(default_accel.a, DoubleEq(0.0));
  EXPECT_THAT(default_accel.b, DoubleEq(0.0));

  shared_ptr<const Acceleration2D> accel2 = make_shared<const Acceleration2D>(0.5, 1.0);
  EXPECT_THAT(accel2->a, DoubleEq(0.5));
  EXPECT_THAT(accel2->b, DoubleEq(1.0));

  shared_ptr<Acceleration2D> accel_copy = make_shared<Acceleration2D>();
  accel_copy->copy(accel2);
  EXPECT_THAT(accel_copy->a, DoubleEq(0.5));
  EXPECT_THAT(accel_copy->b, DoubleEq(1.0));
}


TEST(BorisParticle2D, Instantiation)
{
  Particle2DT default_ctor;
  EXPECT_THAT(default_ctor.DIM, 2);
  EXPECT_THAT(default_ctor.pos->x, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.pos->y, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel->u, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.vel->v, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel->a, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.accel->b, DoubleEq(0.0));
  EXPECT_THAT(default_ctor.mass, DoubleEq(1.0));
  EXPECT_THAT(default_ctor.charge, DoubleEq(1.0));

  Particle2DT special_ctor = Particle2DT(0.5, 1.0);
  EXPECT_THAT(special_ctor.pos->x, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.pos->y, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel->u, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.vel->v, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel->a, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.accel->b, DoubleEq(0.0));
  EXPECT_THAT(special_ctor.mass, DoubleEq(0.5));
  EXPECT_THAT(special_ctor.charge, DoubleEq(1.0));
}

TEST(BorisParticle2D, Copyable)
{
  shared_ptr<Particle2DT> original_ptr = make_shared<Particle2DT>(0.1, 0.2);
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
  Particle2DT expl_copy;

  expl_copy.copy(const_pointer_cast<const Particle2DT>(original_ptr));
  EXPECT_THAT(expl_copy.pos->x, DoubleEq(0.3));
  EXPECT_THAT(expl_copy.pos->y, DoubleEq(0.4));
  EXPECT_THAT(expl_copy.vel->u, DoubleEq(0.5));
  EXPECT_THAT(expl_copy.vel->v, DoubleEq(0.6));
  EXPECT_THAT(expl_copy.accel->a, DoubleEq(0.7));
  EXPECT_THAT(expl_copy.accel->b, DoubleEq(0.8));
  EXPECT_THAT(expl_copy.mass, DoubleEq(0.1));
  EXPECT_THAT(expl_copy.charge, DoubleEq(0.2));

  Particle2DT original = Particle2DT(0.1, 0.2);
  original.pos->x = 0.3;
  original.pos->y = 0.4;
  original.vel->u = 0.5;
  original.vel->v = 0.6;
  original.accel->a = 0.7;
  original.accel->b = 0.8;
  Particle2DT impl_copy = original;
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
