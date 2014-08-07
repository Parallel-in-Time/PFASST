#ifndef _TESTS__EXAMPLES__BORIS__MOCKS__HPP_
#define _TESTS__EXAMPLES__BORIS__MOCKS__HPP_

#include <memory>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../examples/boris/physics.hpp"
#include "../examples/boris/particle.hpp"


template<typename scalar, typename time>
class MockPositionEncap
  : public PositionEncapsulation<scalar, time>
{};


template<typename scalar, typename time>
class MockVelocityEncap
  : public VelocityEncapsulation<scalar, time>
{};


template<typename scalar, typename time>
class MockAccelerationEncap
  : public AccelerationEncapsulation<scalar, time>
{};


template<
  typename scalar,
  typename time
>
class MockParticle
  : ParticleEncapsulation<scalar, time,
                          MockPositionEncap,
                          MockVelocityEncap,
                          MockAccelerationEncap>
{
  private:
    typedef ParticleEncapsulation<scalar,
                                  time,
                                  MockPositionEncap,
                                  MockVelocityEncap,
                                  MockAccelerationEncap> parent_type;
};


template<
  typename scalar,
  typename time,
  template <typename, typename> class ParticleT
>
class MockEField
  : public ElectricField<scalar, time, ParticleT>
{
  private:
    typedef ElectricField<scalar, time, MockParticle> parent_type;

  public:
    using parent_type::evaluate;
    MOCK_METHOD2_T(evaluate, scalar(shared_ptr<typename parent_type::particle_type>, time));
};


template<
  typename scalar,
  typename time,
  template <typename, typename> class ParticleT
>
class MockBField
  : public MagneticField<scalar, time, MockParticle>
{
  private:
    typedef MagneticField<scalar, time, MockParticle> parent_type;

  public:
    using parent_type::evaluate;
    MOCK_METHOD2_T(evaluate, scalar(shared_ptr<typename parent_type::particle_type>, time));
};


template<
  typename scalar,
  typename time
>
class MockEOperator
  : public EnergyOperator<scalar, time,
                          MockParticle<scalar, time>,
                          MockEField<scalar, time, MockParticle>,
                          MockBField<scalar, time, MockParticle>
                         >
{
  private:
    typedef EnergyOperator<scalar, time,
                           MockParticle<scalar, time>,
                           MockEField<scalar, time, MockParticle>,
                           MockBField<scalar, time, MockParticle>
                          > parent_type;

  public:
    using parent_type::evaluate;
    MOCK_METHOD4_T(evaluate, scalar(shared_ptr<const typename parent_type::particle_type>,
                                    time,
                                    shared_ptr<const typename parent_type::e_field_type>,
                                    shared_ptr<const typename parent_type::b_field_type>));
};

#endif  // _TESTS__EXAMPLES__BORIS__MOCKS__HPP_
