/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAP_ENCAP_SWEEPER_HPP_
#define _PFASST_ENCAP_ENCAP_SWEEPER_HPP_

#include <vector>

#include "../interfaces.hpp"
#include "../quadrature.hpp"
#include "encapsulation.hpp"

namespace pfasst
{
  namespace encap
  {

    template<typename time = time_precision>
    class EncapSweeper : public ISweeper<time>
    {
        vector<time> nodes;
        shared_ptr<EncapFactory<time>> factory;

      public:

        void set_nodes(vector<time> nodes)
        {
          this->nodes = nodes;
        }

        const vector<time> get_nodes() const
        {
          return nodes;
        }

        void set_factory(shared_ptr<EncapFactory<time>> factory)
        {
          this->factory = shared_ptr<EncapFactory<time>>(factory);
        }

        shared_ptr<EncapFactory<time>> get_factory() const
        {
          return factory;
        }

        virtual void set_state(shared_ptr<const Encapsulation<time>> q0, size_t m)
        {
          throw NotImplementedYet("sweeper");
        }

        virtual shared_ptr<Encapsulation<time>> get_state(size_t m) const
        {
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        virtual shared_ptr<Encapsulation<time>> get_tau(size_t m) const
        {
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        virtual shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const
        {
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        virtual shared_ptr<Encapsulation<time>> get_end_state()
        {
          return this->get_state(this->get_nodes().size() - 1);
        }

        virtual void evaluate(size_t m)
        {
          throw NotImplementedYet("sweeper");
        }

        virtual void advance()
        {
          throw NotImplementedYet("sweeper");
        }

        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
        {
          throw NotImplementedYet("sweeper");
        }
    };

  }

}

#endif
