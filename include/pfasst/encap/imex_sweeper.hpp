
#ifndef _PFASST_ENCAP_IMEX_SWEEPER_HPP_
#define _PFASST_ENCAP_IMEX_SWEEPER_HPP_

#include <cstdlib>
#include <cassert>
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
          this->Q[m]->copy(q0);
        }

        shared_ptr<Encapsulation<time>> get_state(size_t m) const
        {
          return this->Q[m];
        }

        shared_ptr<Encapsulation<time>> get_tau(size_t m) const
        {
          return this->T[m];
        }

        shared_ptr<Encapsulation<time>> get_saved_state(size_t m) const
        {
          return this->pQ[m];
        }

        virtual void advance()
        {
          this->Q[0]->copy(this->Q.back());
          this->Fe[0]->copy(this->Fe.back());
          this->Fi[0]->copy(this->Fi.back());
        }

        virtual void integrate(time dt, vector<shared_ptr<Encapsulation<time>>> dst) const
        {
          dst[0]->mat_apply(dst, dt, this->Smat, this->Fe, true);
          dst[0]->mat_apply(dst, dt, this->Smat, this->Fi, false);
        }

        void setup(bool coarse)
        {
          auto nodes = this->get_nodes();
          assert(nodes.size() >= 1);

          this->Smat = compute_quadrature(nodes, nodes, 's');

          this->SEmat = this->Smat;
          this->SImat = this->Smat;
          for (size_t m = 0; m < nodes.size() - 1; m++) {
            time ds = nodes[m + 1] - nodes[m];
            this->SEmat(m, m)     -= ds;
            this->SImat(m, m + 1) -= ds;
          }

          for (size_t m = 0; m < nodes.size(); m++) {
            this->Q.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->pQ.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
            this->Fe.push_back(this->get_factory()->create(pfasst::encap::function));
            this->Fi.push_back(this->get_factory()->create(pfasst::encap::function));
          }

          for (size_t m = 0; m < nodes.size() - 1; m++) {
            this->S.push_back(this->get_factory()->create(pfasst::encap::solution));
            if (coarse) {
              this->T.push_back(this->get_factory()->create(pfasst::encap::solution));
            }
          }
        }

        virtual void sweep()
        {
          const auto   nodes  = this->get_nodes();
          const size_t nnodes = nodes.size();
          assert(nnodes >= 1);

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          // integrate
          this->S[0]->mat_apply(this->S, dt, this->SEmat, this->Fe, true);
          this->S[0]->mat_apply(this->S, dt, this->SImat, this->Fi, false);
          if (this->T.size() > 0) {
            for (size_t m = 0; m < nnodes - 1; m++) {
              this->S[m]->saxpy(1.0, this->T[m]);
            }
          }

          // sweep
          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          for (size_t m = 0; m < nnodes - 1; m++) {
            time ds = dt * (nodes[m + 1] - nodes[m]);

            rhs->copy(this->Q[m]);
            rhs->saxpy(ds, this->Fe[m]);
            rhs->saxpy(1.0, this->S[m]);
            this->f_solve(this->Fi[m + 1], this->Q[m + 1], t, ds, rhs);
            this->f_eval_expl(this->Fe[m + 1], this->Q[m + 1], t + ds);

            t += ds;
          }
        }

        virtual void predict(bool initial)
        {
          const auto   nodes  = this->get_nodes();
          const size_t nnodes = nodes.size();
          assert(nnodes >= 1);

          time dt = this->get_controller()->get_time_step();
          time t  = this->get_controller()->get_time();

          if (initial) {
            this->f_eval_expl(this->Fe[0], this->Q[0], t);
            this->f_eval_impl(this->Fi[0], this->Q[0], t);
          }

          shared_ptr<Encapsulation<time>> rhs = this->get_factory()->create(pfasst::encap::solution);

          for (size_t m = 0; m < nnodes - 1; m++) {
            time ds = dt * (nodes[m + 1] - nodes[m]);
            rhs->copy(this->Q[m]);
            rhs->saxpy(ds, this->Fe[m]);
            this->f_solve(this->Fi[m + 1], this->Q[m + 1], t, ds, rhs);
            this->f_eval_expl(this->Fe[m + 1], this->Q[m + 1], t + ds);

            t += ds;
          }
        }

        virtual void save()
        {
          for (size_t m = 0; m < this->pQ.size(); m++) {
            this->pQ[m]->copy(Q[m]);
          }
        }

        virtual void evaluate(size_t m)
        {
          time t = this->get_nodes()[m]; // XXX
          this->f_eval_expl(this->Fe[m], this->Q[m], t);
          this->f_eval_impl(this->Fi[m], this->Q[m], t);
        }

        virtual void f_eval_expl(shared_ptr<Encapsulation<time>> /*f*/,
                                 shared_ptr<Encapsulation<time>> /*q*/,
                                 time /*t*/)
        {
          throw NotImplementedYet("imex (f1eval)");
        }

        virtual void f_eval_impl(shared_ptr<Encapsulation<time>> /*f*/,
                                 shared_ptr<Encapsulation<time>> /*q*/,
                                 time /*t*/)
        {
          throw NotImplementedYet("imex (f2eval)");
        }

        virtual void f_solve(shared_ptr<Encapsulation<time>> /*f*/,
                             shared_ptr<Encapsulation<time>> /*q*/,
                             time /*t*/, time /*dt*/,
                             shared_ptr<Encapsulation<time>> /*rhs*/)
        {
          throw NotImplementedYet("imex (f2comp)");
        }

    };

  }
}

#endif
