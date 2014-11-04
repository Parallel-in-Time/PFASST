#ifndef _TESTS__EXAMPLES__BORIS__MOCKS__HPP_
#define _TESTS__EXAMPLES__BORIS__MOCKS__HPP_

#include <memory>

#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "../examples/boris/physics.hpp"
#include "../examples/boris/particle_3d.hpp"

using namespace pfasst::examples::boris;

template<
  typename scalar,
  typename time
>
using MockPositionEncap = Position3DEncapsulation<scalar, time>;


template<
  typename scalar,
  typename time
>
using MockVelocityEncap = Velocity3DEncapsulation<scalar, time>;


template<
  typename scalar,
  typename time
>
using MockAccelerationEncap = Acceleration3DEncapsulation<scalar, time>;


template<
  typename scalar,
  typename time
>
using MockParticle = Particle3DEncapsulation<scalar, time>;


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
    MOCK_METHOD3_T(evaluate,
                   typename parent_type::particle_type::acceleration_type(vector<shared_ptr<typename parent_type::particle_type>>,
                                                                          size_t,
                                                                          time));
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
    MOCK_METHOD3_T(evaluate,
                   typename parent_type::particle_type::acceleration_type(vector<shared_ptr<typename parent_type::particle_type>>,
                                                                          size_t,
                                                                          time));
};


template<
  typename scalar,
  typename time
>
class MockEOperator
  : public EnergyOperator<scalar, time,
                          MockParticle,
                          MockEField,
                          MockBField
                         >
{
  private:
    typedef EnergyOperator<scalar, time,
                           MockParticle,
                           MockEField,
                           MockBField
                          > parent_type;

  public:
    using parent_type::evaluate;
    MOCK_METHOD2_T(evaluate, scalar(vector<shared_ptr<typename parent_type::particle_type>>,
                                    time));
};

#endif  // _TESTS__EXAMPLES__BORIS__MOCKS__HPP_
