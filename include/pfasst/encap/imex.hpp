
#ifndef _PFASST_IMEX_HPP_
#define _PFASST_IMEX_HPP_

#include <iostream>

#include "encapsulation.hpp"

using namespace std;

namespace pfasst {
  namespace imex {

    using pfasst::encap::Encapsulation;

    template<typename scalar>
    class IMEX : public pfasst::encap::EncapsulatedSweeperMixin<scalar> {
      vector<Encapsulation<scalar>*> Q, pQ, S, Fe, Fi;
      matrix<scalar> Smat, SEmat, SImat;

    public:

      ~IMEX() {
	for (int m=0; m<Q.size(); m++) delete Q[m];
	for (int m=0; m<pQ.size(); m++) delete pQ[m];
	for (int m=0; m<S.size(); m++) delete S[m];
	for (int m=0; m<Fe.size(); m++) delete Fe[m];
	for (int m=0; m<Fi.size(); m++) delete Fi[m];
      }

      void set_q(const Encapsulation<scalar> *q0, unsigned int m)
      {
	Q[m]->copy(q0);
      }

      Encapsulation<scalar>* get_q(unsigned int m) const
      {
	return Q[m];
      }

      Encapsulation<scalar>* get_pq(unsigned int m) const
      {
	return pQ[m];
      }

      void advance() const {
	set_q(get_q(Q.size()-1), 0);
      }

      void setup() {
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
	  pQ.push_back(this->get_factory()->create(pfasst::encap::solution));
	  Fe.push_back(this->get_factory()->create(pfasst::encap::function));
	  Fi.push_back(this->get_factory()->create(pfasst::encap::function));
	}

	for (int m=0; m<nodes.size()-1; m++) {
	  S.push_back(this->get_factory()->create(pfasst::encap::solution));
	}
      }

      virtual void integrate(scalar t0, scalar dt)
      {
	throw NotImplementedYet("imex integrate");
      }

      virtual void residual(scalar t0, scalar dt)
      {
	throw NotImplementedYet("imex residual");
      }

      virtual void sweep(scalar t0, scalar dt)
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
	Encapsulation<scalar> *rhs = this->get_factory()->create(pfasst::encap::solution);

	scalar t = t0;
	for (int m=0; m<nnodes-1; m++) {
	  scalar ds = dt * ( nodes[m+1] - nodes[m] );

	  rhs->copy(Q[m]);
	  rhs->saxpy(ds, Fe[m]);
	  rhs->saxpy(1.0, S[m]);
	  f2comp(Fi[m+1], Q[m+1], t, ds, rhs);
	  f1eval(Fe[m+1], Q[m+1], t + ds);

	  t += ds;
	}

	delete rhs;
      }

      virtual void predict(scalar t0, scalar dt) {
	const auto nodes  = this->get_nodes();
	const int  nnodes = nodes.size();

	f1eval(Fe[0], Q[0], t0);
	f2eval(Fi[0], Q[0], t0);

	Encapsulation<scalar> *rhs = this->get_factory()->create(pfasst::encap::solution);

	scalar t = t0;
	for (int m=0; m<nnodes-1; m++) {
	  scalar ds = dt * ( nodes[m+1] - nodes[m] );
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

      virtual void evaluate(int m) {
	// XXX: time
	f1eval(Fe[m], Q[m], 0.0);
	f2eval(Fi[m], Q[m], 0.0);
      }


      virtual void f1eval(Encapsulation<scalar> *F, Encapsulation<scalar> *Q, scalar t)
      {
	throw NotImplementedYet("imex (f1eval)");
      }

      virtual void f2eval(Encapsulation<scalar> *F, Encapsulation<scalar> *Q, scalar t)
      {
	throw NotImplementedYet("imex (f2eval)");
      }

      virtual void f2comp(Encapsulation<scalar> *F, Encapsulation<scalar> *Q, scalar t, scalar dt, Encapsulation<scalar>* rhs)
      {
	throw NotImplementedYet("imex (f2comp)");
      }

    };

  }
}

#endif
