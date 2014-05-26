
#ifndef _PFASST_IMEX_HPP_
#define _PFASST_IMEX_HPP_

#include <iostream>

#include "pfasst-encapsulated.hpp"
#include "pfasst-quadrature.hpp"

using namespace std;

namespace pfasst {
  namespace imex {

    using pfasst::encap::Encapsulation;

    template<typename time>
    class IMEX : public pfasst::encap::EncapsulatedSweeperMixin<time> {
      vector<Encapsulation*> Q, S, Fe, Fi;
      matrix<time> Smat, SEmat, SImat;

    public:

      ~IMEX() {
	for (int m=0; m<Q.size(); m++) delete Q[m];
	for (int m=0; m<S.size(); m++) delete S[m];
	for (int m=0; m<Fe.size(); m++) delete Fe[m];
	for (int m=0; m<Fi.size(); m++) delete Fi[m];
      }

      void set_q0(Encapsulation *q0)
      {
	Q[0]->copy(q0);
      }

      Encapsulation* get_qend() {
	return Q[Q.size()-1];
      }

      void advance() {
	set_q0(get_qend());
      }

      void setup() {
	auto nodes = this->get_nodes();

	Smat = compute_quadrature(nodes, nodes, 's');

	SEmat = Smat;
	SImat = Smat;
	for (int m=0; m<nodes.size()-1; m++) {
	  time ds = nodes[m+1] - nodes[m];
	  SEmat(m, m)   -= ds;
	  SImat(m, m+1) -= ds;
	}

	for (int m=0; m<nodes.size(); m++) {
	  Q.push_back(this->get_factory()->create(pfasst::encap::solution));
	  Fe.push_back(this->get_factory()->create(pfasst::encap::function));
	  Fi.push_back(this->get_factory()->create(pfasst::encap::function));
	}

	for (int m=0; m<nodes.size()-1; m++) {
	  S.push_back(this->get_factory()->create(pfasst::encap::solution));
	}
      }

      virtual void integrate(time t0, time dt)
      {
	throw NotImplementedYet("imex integrate");
      }

      virtual void residual(time t0, time dt)
      {
	throw NotImplementedYet("imex residual");
      }

      virtual void sweep(time t0, time dt)
      {
	const auto nodes  = this->get_nodes();
	const int  nnodes = nodes.size();

	// integrate
	for (int n=0; n<nnodes-1; n++) {
	  S[n]->setval(0.0);
	  for (int m=0; m<nnodes; m++) {
	    S[n]->saxpy(dt * SEmat(n,m), Fe[m]);
	    S[n]->saxpy(dt * SImat(n,m), Fi[m]);
	  }
	}

	// sweep
	Encapsulation *rhs = this->get_factory()->create(pfasst::encap::solution);

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

      virtual void predict(time t0, time dt) {
	const auto nodes  = this->get_nodes();
	const int  nnodes = nodes.size();

	f1eval(Fe[0], Q[0], t0);
	f2eval(Fi[0], Q[0], t0);

	Encapsulation *rhs = this->get_factory()->create(pfasst::encap::solution);

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

      virtual void f1eval(Encapsulation *F, Encapsulation *Q, time t)
      {
	throw NotImplementedYet("imex (f1eval)");
      }

      virtual void f2eval(Encapsulation *F, Encapsulation *Q, time t)
      {
	throw NotImplementedYet("imex (f2eval)");
      }

      virtual void f2comp(Encapsulation *F, Encapsulation *Q, time t, time dt, Encapsulation* rhs)
      {
	throw NotImplementedYet("imex (f2comp)");
      }

    };

  }
}

#endif
