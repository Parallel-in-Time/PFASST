
#ifndef _PFASST_ENCAP_IMEX_SWEEPER_HPP_
#define _PFASST_ENCAP_IMEX_SWEEPER_HPP_

#include <vector>
#include <memory>

#include <boost/numeric/ublas/matrix.hpp>

#include "encapsulation.hpp"
#include "encap_sweeper.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {
    using pfasst::encap::Encapsulation;

    template<typename time = time_precision>
    class IMEXSweeper
      : public pfasst::encap::EncapSweeper<time>
    {
        vector<shared_ptr<Encapsulation<time>>> Q, pQ, S, T, Fe, Fi;
        matrix<time> Smat, SEmat, SImat;

      public:
        ~IMEXSweeper()
        {}

        void set_state(shared_ptr<const Encapsulation<time>> q0, size_t m)
        {
          Q[m]->copy(q0);
        }

        shared_ptr<Encapsulation<time>> get_state(size_t m) const
        {
          return Q[m];
        }

        shared_ptr<Encapsulation<time>> get_tau(size_t m) const
        {
          return T[m];
        }

        shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const
        {
          return pQ[m];
        }

        virtual void advance()
        {
          Q[0]->copy(Q.back());
          Fe[0]->copy(Fe.back());
          Fi[0]->copy(Fi.back());
        }

        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
        {
          dst[0]->mat_apply(dst, dt, Smat, Fe, true);
          dst[0]->mat_apply(dst, dt, Smat, Fi, false);
        }

        void setup(bool coarse)
        {
          auto nodes = this->get_nodes();
          assert(nodes.size() >= 1);

          Smat = compute_quadrature(nodes, nodes, 's');

          SEmat = Smat;
          SImat = Smat;
          for (size_t m = 0; m < nodes.size() - 1; m++) {
            time ds = nodes[m + 1] - nodes[m];
            SEmat(m, m)     -= ds;
            SImat(m, m + 1) -= ds;
          }

          for (size_t m = 0; m < nodes.size(); m++) {
            Q.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              pQ.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
            Fe.push_back(this->get_factory()->create(pfasst::encap::function));
            Fi.push_back(this->get_factory()->create(pfasst::encap::function));
          }

          for (size_t m = 0; m < nodes.size() - 1; m++) {
            S.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              T.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
          }
        }

        virtual void sweep(time t0, time dt)
        {
          const auto   nodes  = this->get_nodes();
          const size_t nnodes = nodes.size();
          assert(nnodes >= 1);

          // integrate
          S[0]->mat_apply(S, dt, SEmat, Fe, true);
          S[0]->mat_apply(S, dt, SImat, Fi, false);
          if (T.size() > 0) {
            for (size_t m = 0; m < nnodes - 1; m++) {
              S[m]->saxpy(1.0, T[m]);
            }
          }

          // sweep
          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          time t = t0;
          for (size_t m = 0; m < nnodes - 1; m++) {
            time ds = dt * (nodes[m + 1] - nodes[m]);

            rhs->copy(Q[m]);
            rhs->saxpy(ds, Fe[m]);
            rhs->saxpy(1.0, S[m]);
            f2comp(Fi[m + 1], Q[m + 1], t, ds, rhs);
            f1eval(Fe[m + 1], Q[m + 1], t + ds);

            t += ds;
          }
        }

        virtual void predict(time t0, time dt, bool initial)
        {
          const auto   nodes  = this->get_nodes();
          const size_t nnodes = nodes.size();
          assert(nnodes >= 1);

          if (initial) {
            f1eval(Fe[0], Q[0], t0);
            f2eval(Fi[0], Q[0], t0);
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          time t = t0;
          for (size_t m = 0; m < nnodes - 1; m++) {
            time ds = dt * (nodes[m + 1] - nodes[m]);
            rhs->copy(Q[m]);
            rhs->saxpy(ds, Fe[m]);
            f2comp(Fi[m + 1], Q[m + 1], t, ds, rhs);
            f1eval(Fe[m + 1], Q[m + 1], t + ds);

            t += ds;
          }
        }

        virtual void save(bool initial_only)
        {
	  if (initial_only) {
	    pQ[0]->copy(Q[0]);
	  } else {
	    for (size_t m = 0; m < pQ.size(); m++) {
	      pQ[m]->copy(Q[m]);
	    }
	  }
        }

        virtual void evaluate(size_t m)
        {
          time t = this->get_nodes()[m]; // XXX
          f1eval(Fe[m], Q[m], t);
          f2eval(Fi[m], Q[m], t);
        }

        virtual void f1eval(shared_ptr<Encapsulation<time>> F, shared_ptr<Encapsulation<time>> Q,
                            time t)
        {
          throw NotImplementedYet("imex (f1eval)");
        }

        virtual void f2eval(shared_ptr<Encapsulation<time>> F, shared_ptr<Encapsulation<time>> Q,
                            time t)
        {
          throw NotImplementedYet("imex (f2eval)");
        }

        virtual void f2comp(shared_ptr<Encapsulation<time>> F, shared_ptr<Encapsulation<time>> Q,
                            time t, time dt,
                            shared_ptr<Encapsulation<time>> rhs)
        {
          throw NotImplementedYet("imex (f2comp)");
        }

    };

  }
}

#endif
