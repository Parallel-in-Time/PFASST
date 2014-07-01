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

    template<typename ScalarT, typename time>
    class EncapSweeper : public ISweeper
    {
        vector<ScalarT> nodes;
        shared_ptr<EncapFactory<ScalarT, time>> factory;

      public:

        void set_nodes(vector<time> nodes)
        {
          this->nodes = nodes;
        }

        const vector<time> get_nodes() const
        {
          return nodes;
        }

        void set_factory(EncapFactory<ScalarT, time>* factory)
        {
          this->factory = shared_ptr<EncapFactory<ScalarT, time>>(factory);
        }

        EncapFactory<ScalarT, time>* get_factory() const
        {
          return factory.get();
        }

        virtual void set_state(const Encapsulation<ScalarT, time>* q0, unsigned int m)
        {
          throw NotImplementedYet("sweeper");
        }

        virtual Encapsulation<ScalarT, time>* get_state(unsigned int m) const
        {
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        virtual Encapsulation<ScalarT, time>* get_tau(unsigned int m) const
        {
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        virtual Encapsulation<ScalarT, time>* get_saved_state(unsigned int m) const
        {
          throw NotImplementedYet("sweeper");
          return NULL;
        }

        virtual Encapsulation<ScalarT, time>* get_end_state()
        {
          return this->get_state(this->get_nodes().size() - 1);
        }

        virtual void evaluate(int m)
        {
          throw NotImplementedYet("sweeper");
        }

        virtual void advance()
        {
          throw NotImplementedYet("sweeper");
        }

        virtual void integrate(time dt, vector<Encapsulation<ScalarT, time>*> dst) const
        {
          throw NotImplementedYet("sweeper");
        }
    };

  }

}

#endif
