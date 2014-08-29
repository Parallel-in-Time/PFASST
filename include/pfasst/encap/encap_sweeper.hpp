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
      public:
        //! @{
        typedef Encapsulation<time> encap_type;
        typedef EncapFactory<time> factory_type;
        //! @}

      private:
        //! @{
        vector<time> nodes;
        vector<bool> is_proper;
        shared_ptr<factory_type> factory;
        //! @}

      public:
        //! @{
        virtual ~EncapSweeper()
        {}
        //! @}

        //! @{
        virtual void set_nodes(vector<time> nodes)
        {
          auto augmented = pfasst::augment_nodes(nodes);
          this->nodes = get<0>(augmented);
          this->is_proper = get<1>(augmented);
        }

        virtual const vector<time> get_nodes() const
        {
          return nodes;
        }

        const vector<bool> get_is_proper() const
        {
          return is_proper;
        }

        virtual void set_factory(shared_ptr<factory_type> factory)
        {
          this->factory = factory;
        }

        virtual shared_ptr<factory_type> get_factory() const
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
        virtual void set_state(shared_ptr<const encap_type> u0, size_t m)
        {
          UNUSED(u0); UNUSED(m);
          throw NotImplementedYet("sweeper");
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
        {
          return this->get_state(this->get_nodes().size() - 1);
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

        virtual void spread() override
        {
          for (size_t m = 1; m < nodes.size(); m++) {
            this->get_state(m)->copy(this->get_state(0));
          }
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
        virtual void integrate(time dt, vector<shared_ptr<encap_type>> dst) const
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
