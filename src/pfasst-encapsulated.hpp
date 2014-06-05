/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>

#include "pfasst-interfaces.hpp"
#include "pfasst-quadrature.hpp"

using namespace std;

namespace pfasst {
  namespace encap {

    typedef enum EncapType { solution, function } EncapType;

    //
    // encapsulation
    //
    template<typename scalar>
    class Encapsulation {
    public:
      virtual ~Encapsulation() { }

      // required for time-parallel communications
      virtual void send() {
	throw NotImplementedYet("pfasst");
      }
      virtual void recv() {
	throw NotImplementedYet("pfasst");
      }

      // required for host based encap helpers
      virtual void setval(scalar) {
	throw NotImplementedYet("encap");
      }
      virtual void copy(const Encapsulation<scalar> *) {
	throw NotImplementedYet("encap");
      }
      virtual void saxpy(scalar a, const Encapsulation<scalar> *) {
	throw NotImplementedYet("encap");
      }
      virtual void mat_apply(vector<Encapsulation<scalar>*> dst, scalar a, matrix<scalar> m,
			     vector<Encapsulation<scalar>*> src, bool zero=true) {
        throw NotImplementedYet("encap");
      }
    };

    template<typename scalar>
    class EncapsulationFactory {
    public:
      virtual Encapsulation<scalar>* create(const EncapType) = 0;
    };

    template<typename scalar>
    class EncapsulatedSweeperMixin : public ISweeper {
      vector<scalar> nodes;
      shared_ptr<EncapsulationFactory<scalar>> factory;

    public:
      void set_nodes(vector<scalar> nodes) {
	this->nodes = nodes;      //
      }

      const vector<scalar> get_nodes() const {
	return nodes;
      }

      void set_factory(EncapsulationFactory<scalar>* factory) {
	this->factory = shared_ptr<EncapsulationFactory<scalar>>(factory);
      }

      EncapsulationFactory<scalar>* get_factory() const {
	return factory.get();
      }

      virtual void set_q(const Encapsulation<scalar>* q0, unsigned int m) {
	throw NotImplementedYet("sweeper");
      }

      virtual Encapsulation<scalar>* get_q(unsigned int m) const {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar>* get_pq(unsigned int m) const {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar>* get_qend() {
	return this->get_q(this->get_nodes().size()-1);
      }

      virtual void evaluate(int m) {
	throw NotImplementedYet("sweeper");
      }


      virtual void advance() {
	this->set_q(this->get_qend(), 0);
      }
    };

    template<typename scalar>
    class PolyInterpMixin : public pfasst::ITransfer {
      matrix<scalar> tmat;
    public:

      virtual void interpolate(ISweeper *DST, const ISweeper *SRC, bool initial) {
	auto* dst = dynamic_cast<EncapsulatedSweeperMixin<scalar>*>(DST);
	auto* src = dynamic_cast<const EncapsulatedSweeperMixin<scalar>*>(SRC);

	if (tmat.size() == 0)
	  tmat = pfasst::compute_interp<scalar>(dst->get_nodes(), src->get_nodes());

	int ndst = dst->get_nodes().size();
	int nsrc = src->get_nodes().size();

	auto* crse_factory = src->get_factory();
	auto* fine_factory = dst->get_factory();

	vector<Encapsulation<scalar>*> fine_q(ndst), fine_tmp(ndst);

	for (int m=0; m<ndst; m++) fine_q[m]   = dst->get_q(m);
	for (int m=0; m<ndst; m++) fine_tmp[m] = fine_factory->create(solution);

	if (initial)
	  for (int m=1; m<ndst; m++)
	    fine_q[m]->copy(fine_q[0]);

	auto* crse_tmp = crse_factory->create(solution);
	for (int m=0; m<nsrc; m++) {
	  crse_tmp->copy(src->get_q(m));
	  if (initial)
	    crse_tmp->saxpy(-1.0, src->get_q(0));
	  else
	    crse_tmp->saxpy(-1.0, src->get_pq(m));
	  interpolate(fine_tmp[m], crse_tmp);
	}
	delete crse_tmp;

	dst->get_q(0)->mat_apply(fine_q, 1.0, tmat, fine_tmp, false);

	for (int m=0; m<ndst; m++) delete fine_tmp[m];
	for (int m=0; m<ndst; m++) dst->evaluate(m);


      }

      virtual void restrict(ISweeper *DST, const ISweeper *SRC) {
	auto* dst = dynamic_cast<EncapsulatedSweeperMixin<scalar>*>(DST);
	auto* src = dynamic_cast<const EncapsulatedSweeperMixin<scalar>*>(SRC);

	auto dnodes = dst->get_nodes();
	auto snodes = src->get_nodes();

	int ndst = dst->get_nodes().size();
	int nsrc = src->get_nodes().size();

	int trat = (nsrc - 1) / (ndst - 1);

	for (int m=0; m<ndst; m++) {
	  if (dnodes[m] != snodes[m*trat])
	    throw NotImplementedYet("coarse nodes must be nested");
	  this->restrict(dst->get_q(m), src->get_q(m*trat));
	}

	for (int m=0; m<ndst; m++) dst->evaluate(m);
      }

      // required for interp/restrict helpers
      virtual void interpolate(Encapsulation<scalar> *dst, const Encapsulation<scalar> *src) {
	throw NotImplementedYet("mlsdc/pfasst");
      }

      virtual void restrict(Encapsulation<scalar> *dst, const Encapsulation<scalar> *src) {
	throw NotImplementedYet("mlsdc/pfasst");
      }

    };

  }
}

#endif
