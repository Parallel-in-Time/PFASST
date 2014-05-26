
#ifndef _PFASST_IMEX_HPP_
#define _PFASST_IMEX_HPP_

#include <iostream>

#include "pfasst-encapsulated.hpp"
#include "pfasst-quadrature.hpp"

using namespace std;

namespace pfasst {

  template<typename timeT>
  class IMEX : public EncapsulatedSweeperMixin<timeT> {
    vector<Encapsulation*> Q, S, Fe, Fi;
    matrix<timeT> Smat, SEmat, SImat;

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

      SEmat = Smat; SImat = Smat;
      for (int m=0; m<nodes.size()-1; m++) {
	timeT dt = nodes[m+1] - nodes[m];
	SEmat(m, m)   -= dt;
	SImat(m, m+1) -= dt;
      }

      for (int m=0; m<nodes.size(); m++) {
      	this->Q.push_back(this->get_factory()->create(solution));
      	this->Fe.push_back(this->get_factory()->create(function));
      	this->Fi.push_back(this->get_factory()->create(function));
      }

      for (int m=0; m<nodes.size()-1; m++) {
      	S.push_back(this->get_factory()->create(solution));
      }
    }

    virtual void integrate(timeT t0, timeT dt0)
    {
      throw NotImplementedYet("imex integrate");
    }
    virtual void residual(timeT t0, timeT dt0)
    {
      throw NotImplementedYet("imex residual");
    }

    virtual void sweep(timeT t0, timeT dt0)
    {
      const auto nodes  = this->get_nodes();
      const int  nnodes = nodes.size();

      // integrate
      for (int n=0; n<nnodes-1; n++) {
	this->S[n]->setval(0.0);
	for (int m=0; m<nnodes; m++) {
	  this->S[n]->saxpy(dt0 * this->SEmat(n,m), this->Fe[m]);
	  this->S[n]->saxpy(dt0 * this->SImat(n,m), this->Fi[m]);
	}
      }

      // sweep
      Encapsulation *rhs = this->get_factory()->create(solution);

      timeT t = t0;
      for (int m=0; m<nnodes-1; m++) {
	timeT dt = dt0 * ( nodes[m+1] - nodes[m] );

      	rhs->copy(this->Q[m]);
      	rhs->saxpy(dt, this->Fe[m]);
      	rhs->saxpy(1.0, this->S[m]);
      	this->f2comp(this->Fi[m+1], this->Q[m+1], t, dt, rhs);
      	this->f1eval(this->Fe[m+1], this->Q[m+1], t+dt);

	t += dt;
      }

      delete rhs;
    }

    virtual void predict(timeT t0, timeT dt0) {
      const auto nodes  = this->get_nodes();
      const int  nnodes = nodes.size();

      this->f1eval(this->Fe[0], this->Q[0], t0);
      this->f2eval(this->Fi[0], this->Q[0], t0);

      Encapsulation *rhs = this->get_factory()->create(solution);

      timeT t = t0;
      for (int m=0; m<nnodes-1; m++) {
	timeT dt = dt0 * ( nodes[m+1] - nodes[m] );
	rhs->copy(this->Q[m]);
	rhs->saxpy(dt, this->Fe[m]);
	this->f2comp(this->Fi[m+1], this->Q[m+1], t, dt, rhs);
	this->f1eval(this->Fe[m+1], this->Q[m+1], t+dt);

	t += dt;
      }

      delete rhs;
    }

    virtual void f1eval(Encapsulation *F, Encapsulation *Q, timeT t)
    {
      throw NotImplementedYet("imex (f1eval)");
    }

    virtual void f2eval(Encapsulation *F, Encapsulation *Q, timeT t)
    {
      throw NotImplementedYet("imex (f2eval)");
    }

    virtual void f2comp(Encapsulation *F, Encapsulation *Q, timeT t, timeT dt, Encapsulation* rhs)
    {
      throw NotImplementedYet("imex (f2comp)");
    }

  };

}

#endif
