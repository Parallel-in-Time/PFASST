
#ifndef _PFASST_ENCAP_IMEX_SWEEPER_HPP_
#define _PFASST_ENCAP_IMEX_SWEEPER_HPP_

#include <iostream>

#include "encap_sweeper.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {

    using pfasst::encap::Encapsulation;

    template<typename ScalarT, typename timeT>
    class IMEXSweeper : public pfasst::encap::EncapSweeper<ScalarT, timeT>
    {
        vector<Encapsulation<ScalarT, timeT>*> Q, pQ, S, T, Fe, Fi;
        matrix<timeT> Smat, SEmat, SImat;

      public:

        ~IMEXSweeper()
        {
          for (int m = 0; m < Q.size(); m++)  { delete Q[m]; }

          for (int m = 0; m < S.size(); m++)  { delete S[m]; }

          for (int m = 0; m < T.size(); m++)  { delete T[m]; }

          for (int m = 0; m < pQ.size(); m++) { delete pQ[m]; }

          for (int m = 0; m < Fe.size(); m++) { delete Fe[m]; }

          for (int m = 0; m < Fi.size(); m++) { delete Fi[m]; }
        }

        void set_state(const Encapsulation<ScalarT, timeT>* q0, unsigned int m)
        {
          Q[m]->copy(q0);

        }

        Encapsulation<ScalarT, timeT>* get_state(unsigned int m) const
        {
          return Q[m];
        }

        Encapsulation<ScalarT, timeT>* get_tau(unsigned int m) const
        {
          return T[m];
        }

        Encapsulation<ScalarT, timeT>* get_saved_state(unsigned int m) const
        {
          return pQ[m];
        }

        virtual void advance()
        {
          Q[0]->copy(Q[Q.size() - 1]);
          Fe[0]->copy(Fe[Fe.size() - 1]);
          Fi[0]->copy(Fi[Fi.size() - 1]);
        }


        virtual void integrate(timeT dt, vector<Encapsulation<ScalarT, timeT>*> dst) const
        {
          auto* encap = dst[0];
          encap->mat_apply(dst, dt, Smat, Fe, true);
          encap->mat_apply(dst, dt, Smat, Fi, false);
        }

        void setup(bool coarse)
        {
          auto nodes = this->get_nodes();

          Smat = compute_quadrature(nodes, nodes, 's');

          SEmat = Smat;
          SImat = Smat;

          for (int m = 0; m < nodes.size() - 1; m++) {
            ScalarT ds = nodes[m + 1] - nodes[m];
            SEmat(m, m)   -= ds;
            SImat(m, m + 1) -= ds;
          }

          for (int m = 0; m < nodes.size(); m++) {
            Q.push_back(this->get_factory()->create(pfasst::encap::solution));

            if (coarse)
            { pQ.push_back(this->get_factory()->create(pfasst::encap::solution)); }

            Fe.push_back(this->get_factory()->create(pfasst::encap::function));
            Fi.push_back(this->get_factory()->create(pfasst::encap::function));
          }

          for (int m = 0; m < nodes.size() - 1; m++) {
            S.push_back(this->get_factory()->create(pfasst::encap::solution));

            if (coarse)
            { T.push_back(this->get_factory()->create(pfasst::encap::solution)); }
          }
        }

        virtual void sweep(timeT t0, timeT dt)
        {
          const auto nodes  = this->get_nodes();
          const int  nnodes = nodes.size();

          // integrate
          S[0]->mat_apply(S, dt, SEmat, Fe, true);
          S[0]->mat_apply(S, dt, SImat, Fi, false);

          if (T.size() > 0)
            for (int m = 0; m < nnodes - 1; m++)
            { S[m]->saxpy(1.0, T[m]); }


          // sweep
          Encapsulation<ScalarT, timeT>* rhs = this->get_factory()->create(pfasst::encap::solution);

          timeT t = t0;

          for (int m = 0; m < nnodes - 1; m++) {
            timeT ds = dt * (nodes[m + 1] - nodes[m]);

            rhs->copy(Q[m]);
            rhs->saxpy(ds, Fe[m]);
            rhs->saxpy(1.0, S[m]);
            f2comp(Fi[m + 1], Q[m + 1], t, ds, rhs);
            f1eval(Fe[m + 1], Q[m + 1], t + ds);

            t += ds;
          }

          delete rhs;
        }

        virtual void predict(timeT t0, timeT dt, bool initial)
        {
          const auto nodes  = this->get_nodes();
          const int  nnodes = nodes.size();

          if (initial) {
            f1eval(Fe[0], Q[0], t0);
            f2eval(Fi[0], Q[0], t0);
          }

          Encapsulation<ScalarT, timeT>* rhs = this->get_factory()->create(pfasst::encap::solution);

          timeT t = t0;

          for (int m = 0; m < nnodes - 1; m++) {
            timeT ds = dt * (nodes[m + 1] - nodes[m]);
            rhs->copy(Q[m]);
            rhs->saxpy(ds, Fe[m]);
            f2comp(Fi[m + 1], Q[m + 1], t, ds, rhs);
            f1eval(Fe[m + 1], Q[m + 1], t + ds);

            t += ds;
          }

          delete rhs;
        }

        virtual void save()
        {
          for (int m = 0; m < pQ.size(); m++)
          { pQ[m]->copy(Q[m]); }
        }

        virtual void evaluate(int m)
        {
          timeT t = this->get_nodes()[m]; // XXX
          f1eval(Fe[m], Q[m], t);
          f2eval(Fi[m], Q[m], t);
        }


        virtual void f1eval(Encapsulation<ScalarT, timeT>* F, Encapsulation<ScalarT, timeT>* Q, timeT t)
        {
          throw NotImplementedYet("imex (f1eval)");
        }

        virtual void f2eval(Encapsulation<ScalarT, timeT>* F, Encapsulation<ScalarT, timeT>* Q, timeT t)
        {
          throw NotImplementedYet("imex (f2eval)");
        }

        virtual void f2comp(Encapsulation<ScalarT, timeT>* F, Encapsulation<ScalarT, timeT>* Q, timeT t, timeT dt,
                            Encapsulation<ScalarT, timeT>* rhs)
        {
          throw NotImplementedYet("imex (f2comp)");
        }

    };

  }
}

#endif