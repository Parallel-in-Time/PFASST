
#ifndef _PFASST_IMEX_HPP_
#define _PFASST_IMEX_HPP_

#include <iostream>

#include "encapsulation.hpp"

using namespace std;

namespace pfasst {
  namespace imex {

    using pfasst::encap::Encapsulation;

    template<typename scalar, typename time>
    class IMEX : public pfasst::encap::EncapsulatedSweeperMixin<scalar,time> {
      vector<Encapsulation<scalar,time>*> Q, pQ, S, T, Fe, Fi;
      matrix<time> Smat, SEmat, SImat;

    public:

      ~IMEX()
      {
	for (int m=0; m<Q.size(); m++)  delete Q[m];
	for (int m=0; m<S.size(); m++)  delete S[m];
	for (int m=0; m<T.size(); m++)  delete T[m];
	for (int m=0; m<pQ.size(); m++) delete pQ[m];
	for (int m=0; m<Fe.size(); m++) delete Fe[m];
	for (int m=0; m<Fi.size(); m++) delete Fi[m];
      }

      void set_state(const Encapsulation<scalar,time> *q0, unsigned int m)
      {
	Q[m]->copy(q0);
      }

      Encapsulation<scalar,time>* get_state(unsigned int m) const
      {
	return Q[m];
      }

      Encapsulation<scalar,time>* get_tau(unsigned int m) const
      {
	return T[m];
      }

      Encapsulation<scalar,time>* get_saved_state(unsigned int m) const
      {
	return pQ[m];
      }

      virtual void integrate(time dt, vector<Encapsulation<scalar,time>*> dst) const
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
	for (int m=0; m<nodes.size()-1; m++) {
	  scalar ds = nodes[m+1] - nodes[m];
	  SEmat(m, m)   -= ds;
	  SImat(m, m+1) -= ds;
	}

	for (int m=0; m<nodes.size(); m++) {
	  Q.push_back(this->get_factory()->create(pfasst::encap::solution));
	  if (coarse)
	    pQ.push_back(this->get_factory()->create(pfasst::encap::solution));
	  Fe.push_back(this->get_factory()->create(pfasst::encap::function));
	  Fi.push_back(this->get_factory()->create(pfasst::encap::function));
	}

	for (int m=0; m<nodes.size()-1; m++) {
	  S.push_back(this->get_factory()->create(pfasst::encap::solution));
	  if (coarse)
	    T.push_back(this->get_factory()->create(pfasst::encap::solution));
	}
      }

      virtual void sweep(time t0, time dt)
      {
	const auto nodes  = this->get_nodes();
	const int  nnodes = nodes.size();

	// integrate
	S[0]->mat_apply(S, dt, SEmat, Fe, true);
	S[0]->mat_apply(S, dt, SImat, Fi, false);
	if (T.size() > 0)
	  for (int m=0; m<nnodes-1; m++)
	    S[m]->saxpy(1.0, T[m]);


	// sweep
	Encapsulation<scalar,time> *rhs = this->get_factory()->create(pfasst::encap::solution);

	time t = t0;
	for (int m=0; m<nnodes-1; m++) {
	  time ds = dt * ( nodes[m+1] - nodes[m] );

	  rhs->copy(Q[m]);
	  rhs->saxpy(ds, Fe[m]);
	  rhs->saxpy(1.0, S[m]);
	  f2comp(Fi[m+1], Q[m+1], t, ds, rhs);
	  f1eval(Fe[m+1], Q[m+1], t + ds);

	  t += ds;
	}

	delete rhs;
      }

      virtual void predict(time t0, time dt)
      {
	const auto nodes  = this->get_nodes();
	const int  nnodes = nodes.size();

	f1eval(Fe[0], Q[0], t0);
	f2eval(Fi[0], Q[0], t0);

	Encapsulation<scalar,time> *rhs = this->get_factory()->create(pfasst::encap::solution);

	time t = t0;
	for (int m=0; m<nnodes-1; m++) {
	  time ds = dt * ( nodes[m+1] - nodes[m] );
	  rhs->copy(Q[m]);
	  rhs->saxpy(ds, Fe[m]);
	  f2comp(Fi[m+1], Q[m+1], t, ds, rhs);
	  f1eval(Fe[m+1], Q[m+1], t + ds);

	  t += ds;
	}

	delete rhs;
      }

      virtual void save()
      {
	for (int m=0; m<pQ.size(); m++)
	  pQ[m]->copy(Q[m]);
      }

      virtual void evaluate(int m)
      {
	// XXX: time
	f1eval(Fe[m], Q[m], 0.0);
	f2eval(Fi[m], Q[m], 0.0);
      }


      virtual void f1eval(Encapsulation<scalar,time> *F, Encapsulation<scalar,time> *Q, time t)
      {
	throw NotImplementedYet("imex (f1eval)");
      }

      virtual void f2eval(Encapsulation<scalar,time> *F, Encapsulation<scalar,time> *Q, time t)
      {
	throw NotImplementedYet("imex (f2eval)");
      }

      virtual void f2comp(Encapsulation<scalar,time> *F, Encapsulation<scalar,time> *Q, time t, time dt, Encapsulation<scalar,time>* rhs)
      {
	throw NotImplementedYet("imex (f2comp)");
      }

    };

  }
}

#endif
