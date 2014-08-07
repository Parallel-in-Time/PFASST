#ifndef _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
#define _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_

#include <cstdlib>
#include <map>
#include <cassert>

#include <pfasst/encap/encap_sweeper.hpp>

#include "particle_2d.hpp"
#include "physics.hpp"

using namespace std;
using namespace pfasst;
using namespace pfasst::encap;

typedef map<pair<size_t, size_t>, double> error_map;


template<
  typename scalar,
  typename time = time_precision
>
class BorisSweeper
  : public EncapSweeper<time>
{
  public:
    typedef Particle2DEncapsulation<scalar, time> encap_type;
    typedef Particle2DFactory<scalar, time> factory_type;
    typedef ElectricField<scalar, time, Particle2DEncapsulation> e_field_type;
    typedef MagneticField<scalar, time, Particle2DEncapsulation> b_field_type;

  private:
    shared_ptr<e_field_type> e_field;
    shared_ptr<b_field_type> b_field;
    error_map errors;

  protected:
    vector<shared_ptr<encap_type>> particles;

  public:
    //! @{
    BorisSweeper()
      :   e_field(nullptr)
        , b_field(nullptr)
        , errors()
    {}

    BorisSweeper(const BorisSweeper<scalar, time>& other) = delete;
    BorisSweeper(BorisSweeper<scalar, time>&& other) = delete;

    virtual ~BorisSweeper()
    {}
    //! @}

    //! @{
    virtual void set_state(shared_ptr<const Encapsulation<time>> u0, size_t m) override
    {
      shared_ptr<const encap_type> u0_cast = dynamic_pointer_cast<const encap_type>(u0);
      assert(u0_cast);
      this->set_state(u0_cast, m);
    }
    virtual void set_state(shared_ptr<const encap_type> u0, size_t m)
    {
      UNUSED(u0); UNUSED(m);
      // TODO: implement set_state for BorisSweeper
    }
    virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const override
    {
      UNUSED(m);
      // TODO: implement get_state for BorisSweeper
      return NULL;
    }

    virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const override
    {
      UNUSED(m);
      // TODO: implement get_tau for BorisSweeper
      return NULL;
    }

    virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const override
    {
      UNUSED(m);
      // TODO: implement get_saved_state for BorisSweeper
      return NULL;
    }
    //! @}

    //! @{
    virtual void set_e_field(shared_ptr<const e_field_type> e_field)
    {
      this->e_field = const_pointer_cast<e_field_type>(e_field);
    }
    virtual shared_ptr<e_field_type> get_e_field()
    {
      return this->e_field;
    }
    virtual void set_b_field(shared_ptr<const b_field_type> b_field)
    {
      this->b_field = const_pointer_cast<b_field_type>(b_field);
    }
    virtual shared_ptr<b_field_type> get_b_field()
    {
      return this->b_field;
    }
    //! @}

    //! @{
    virtual void exact(shared_ptr<Encapsulation<time>> q, time t)
    {
      shared_ptr<encap_type> q_cast = dynamic_pointer_cast<encap_type>(q);
      assert(q_cast);
      this->exact(q_cast, t);
    }
    virtual void exact(shared_ptr<encap_type> q, time t)
    {
      UNUSED(q); UNUSED(t);
      // TODO: implement exact solution for Boris
    }

    virtual void echo_error(time t, bool predict = false) const
    {
      UNUSED(t); UNUSED(predict);
      // TODO: implement echo_error
    }

    virtual error_map get_errors() const
    {
      return this->errors;
    }
    //! @}

    //! @{
    virtual void setup(bool coarse = false) override
    {
      UNUSED(coarse);
      // TODO: implement setup for BorisSweeper
    }
    virtual void advance() override
    {
      // TODO: implement advance for BorisSweeper
    }
    virtual void evaluate(size_t m) override
    {
      UNUSED(m);
      // TODO: implement evaluate for BorisSweeper
    }
    virtual void predict(bool initial) override
    {
      UNUSED(initial);
      // TODO: implement predict for BorisSweeper
    }

    virtual void sweep() override
    {
      // TODO: implement sweep for BorisSweeper
    }
    
    virtual void save(bool initial_only=false) override
    {
      UNUSED(initial_only);
      // TODO: implement save for BorisSweeper
    }
    virtual void spread() override
    {
      // TODO: implement spread for BorisSweeper
    }
    //! @}

    //! @{
      virtual void post(ICommunicator* comm, int tag) override
      {
        UNUSED(comm); UNUSED(tag);
      };

      virtual void send(ICommunicator* comm, int tag, bool blocking) override
      {
        UNUSED(comm); UNUSED(tag); UNUSED(blocking);
        NotImplementedYet("pfasst");
      }

      virtual void recv(ICommunicator* comm, int tag, bool blocking) override
      {
        UNUSED(comm); UNUSED(tag); UNUSED(blocking);
        NotImplementedYet("pfasst");
      }

      virtual void broadcast(ICommunicator* comm) override
      {
        UNUSED(comm);
        NotImplementedYet("pfasst");
      }
    //! @}
};

#endif  // _EXAMPLES__BORIS__BORIS_SWEEPER__HPP_
