/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAP_ENCAP_SWEEPER_HPP_
#define _PFASST_ENCAP_ENCAP_SWEEPER_HPP_

#include <cstdlib>
#include <vector>
#include <memory>

#include "../globals.hpp"
#include "../interfaces.hpp"
#include "../quadrature.hpp"
#include "encapsulation.hpp"

namespace pfasst
{
  namespace encap
  {
    template<typename time = time_precision>
    class EncapSweeper
      : public ISweeper<time>
    {
      protected:
        //! @{
        quadrature::IQuadrature<time>* quad;
        shared_ptr<EncapFactory<time>> factory;
        shared_ptr<Encapsulation<time>> u_start;
        shared_ptr<Encapsulation<time>> u_end;
        shared_ptr<Encapsulation<time>> u_end_old;
        //! @}

      public:
        //! @{
        EncapSweeper()
          : quad(nullptr)
        {}

        EncapSweeper(const EncapSweeper<time>& other)
          :   quad(other.quad)
            , factory(other.factory)
            , u_start(other.u_start)
            , u_end(other.u_end)
            , u_end_old(other.u_end_old)
        {}

        EncapSweeper(EncapSweeper<time>&& other)
          : EncapSweeper<time>()
        {
          swap(*this, other);
        }

        virtual ~EncapSweeper()
        {
          if (this->quad) delete this->quad;
        }
        //! @}

        //! @{
        EncapSweeper<time>& operator=(EncapSweeper<time> other)
        {
          swap(*this, other);
          return *this;
        }

        friend void swap(EncapSweeper<time>& first, EncapSweeper<time>& second)
        {
          pfasst::ISweeper<time>::swap(first, second);
          using std::swap;
          swap(first.quad, second.quad);
          swap(first.factory, second.factory);
          swap(first.u_start, second.u_start);
          swap(first.u_end, second.u_end);
          swap(first.u_end_old, second.u_end_old);
        }
        //! @}

        //! @{
        virtual void spread() override
        {
          for (size_t m = 1; m < this->quad->get_num_nodes(); m++) {
            this->get_state(m)->copy(this->u_start);
          }
        }
        //! @}

        //! @{
        void set_quadrature(quadrature::IQuadrature<time>* quad)
        {
          this->quad = quad;
        }

        const quadrature::IQuadrature<time>* get_quadrature()
        {
          return this->quad;
        }

        shared_ptr<Encapsulation<time>> get_u_start() const
        {
          return this->u_start;
        }

        shared_ptr<Encapsulation<time>> get_u_end() const
        {
          return this->u_end;
        }

        shared_ptr<const Encapsulation<time>> get_u_end_old() const
        {
          return this->u_end_old;
        }

        const vector<time> get_nodes() const
        {
          return this->quad->get_nodes();
        }

        void set_factory(shared_ptr<EncapFactory<time>> factory)
        {
          this->factory = factory;
        }

        virtual shared_ptr<EncapFactory<time>> get_factory() const
        {
          return factory;
        }

        /**
         * sets solution values at time node index `m`
         *
         * @param[in] u0 values to set
         * @param[in] m 0-based node index
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void set_state(shared_ptr<const Encapsulation<time>> u0, size_t m)
        {
          UNUSED(u0); UNUSED(m);
          throw NotImplementedYet("sweeper");
        }

        void set_u_start(shared_ptr<const Encapsulation<time>> u_start)
        {
          this->u_start = u_start;
        }

        void set_u_end(shared_ptr<const Encapsulation<time>> u_end)
        {
          this->u_end = u_end;
        }

        /**
         * retrieve solution values of current iteration at time node index `m`
         *
         * @param[in] m 0-based index of time node
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const
        {
          UNUSED(m);
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        /**
         * retrieves FAS correction of current iteration at time node index `m`
         *
         * @param[in] m 0-based index of time node
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const
        {
          UNUSED(m);
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        /**
         * retrieves solution values of previous iteration at time node index `m`
         *
         * @param[in] m 0-based index of time node
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const
        {
          UNUSED(m);
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        virtual shared_ptr<Encapsulation<time>> get_end_state()
        __attribute__((deprecated("use u_end")))
        {
          return this->get_u_end();
        }
        //! @}

        //! @{
        /**
         * @copydoc ISweeper::advance()
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void advance() override
        {
          throw NotImplementedYet("sweeper");
        }

        virtual void do_inner_nodes(const size_t m, const time t, const time ds) override
        {
          UNUSED(m); UNUSED(t); UNUSED(ds);
          throw NotImplementedYet("sweeper");
        }

        virtual void do_first_node(const time t, const time ds) override
        {
          UNUSED(t); UNUSED(ds);
          throw NotImplementedYet("sweeper");
        }

        virtual void do_last_point(const time t, const time dt) override
        {
          UNUSED(t); UNUSED(dt);
          throw NotImplementedYet("sweeper");
        }

        /**
         * evaluates the right hand side at given time node
         *
         * This evaluates the right hand side at the given time node with index `m` as returned by
         * pfasst::encap::EncapSweeper::get_nodes:
         *
         * @param[in] m index of the time node to evaluate at
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void evaluate(size_t m)
        {
          UNUSED(m);
          throw NotImplementedYet("sweeper");
        }

        /**
         * integrates values of right hand side at all time nodes \\( t \\in [0,M-1] \\)
         * simultaneously
         *
         * @param[in] dt width of the time interval to integrate
         * @param[in,out] dst integrated values
         *
         * @note This method must be implemented in derived sweepers.
         */
        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
        {
          UNUSED(dt); UNUSED(dst);
          throw NotImplementedYet("sweeper");
        }
        //! @}

        //! @{
        virtual void post(ICommunicator* comm, int tag) override
        {
          this->get_state(0)->post(comm, tag);
        }

        virtual void send(ICommunicator* comm, int tag, bool blocking) override
        {
          this->get_state(this->get_nodes().size() - 1)->send(comm, tag, blocking);
        }

        virtual void recv(ICommunicator* comm, int tag, bool blocking) override
        {
          this->get_state(0)->recv(comm, tag, blocking);
        }

        virtual void broadcast(ICommunicator* comm) override
        {
          if (comm->rank() == comm->size() - 1) {
            this->get_state(0)->copy(this->get_state(this->get_nodes().size() - 1));
          }
          this->get_state(0)->broadcast(comm);
        }
        //! @}
    };


    template<typename time>
    EncapSweeper<time>& as_encap_sweeper(shared_ptr<ISweeper<time>> x)
    {
      shared_ptr<EncapSweeper<time>> y = dynamic_pointer_cast<EncapSweeper<time>>(x);
      assert(y);
      return *y.get();
    }


    template<typename time>
    const EncapSweeper<time>& as_encap_sweeper(shared_ptr<const ISweeper<time>> x)
    {
      shared_ptr<const EncapSweeper<time>> y = dynamic_pointer_cast<const EncapSweeper<time>>(x);
      assert(y);
      return *y.get();
    }

  }  // ::pfasst::encap
}  // ::pfasst

#endif
