/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>

#include "../interfaces.hpp"
#include "../quadrature.hpp"

using namespace std;

namespace pfasst {
  namespace encap {

    typedef enum EncapType { solution, function } EncapType;

    //
    // encapsulation
    //
    template<typename scalar, typename time>
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
      virtual void copy(const Encapsulation<scalar,time> *) {
	throw NotImplementedYet("encap");
      }
      virtual void saxpy(time a, const Encapsulation<scalar,time> *) {
	throw NotImplementedYet("encap");
      }
      virtual void mat_apply(vector<Encapsulation<scalar,time>*> dst, time a, matrix<time> m,
			     vector<Encapsulation<scalar,time>*> src, bool zero=true) {
        throw NotImplementedYet("encap");
      }
    };

    template<typename scalar, typename time>
    class EncapsulationFactory {
    public:
      virtual Encapsulation<scalar,time>* create(const EncapType) = 0;
    };

    template<typename scalar, typename time>
    class EncapsulatedSweeperMixin : public ISweeper {
      vector<scalar> nodes;
      shared_ptr<EncapsulationFactory<scalar,time>> factory;

    public:
      void set_nodes(vector<time> nodes) {
	this->nodes = nodes;
      }

      const vector<time> get_nodes() const {
	return nodes;
      }

      void set_factory(EncapsulationFactory<scalar,time>* factory) {
	this->factory = shared_ptr<EncapsulationFactory<scalar,time>>(factory);
      }

      EncapsulationFactory<scalar,time>* get_factory() const {
	return factory.get();
      }

      virtual void set_q(const Encapsulation<scalar,time>* q0, unsigned int m) {
	throw NotImplementedYet("sweeper");
      }

      virtual Encapsulation<scalar,time>* get_q(unsigned int m) const {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar,time>* get_tau(unsigned int m) const {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar,time>* get_pq(unsigned int m) const {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar,time>* get_qend() {
	return this->get_q(this->get_nodes().size()-1);
      }

      virtual void evaluate(int m) {
	throw NotImplementedYet("sweeper");
      }

      virtual void advance() {
	this->set_q(this->get_qend(), 0);
      }

      virtual void integrate(Encapsulation<scalar,time>* dst, time dt) {
	throw NotImplementedYet("sweeper");
      }
    };

    template<typename scalar, typename time>
    class PolyInterpMixin : public pfasst::ITransfer {
      matrix<time> tmat;
    public:

      virtual void interpolate(ISweeper *DST, const ISweeper *SRC, bool initial)
      {
	auto* dst = dynamic_cast<EncapsulatedSweeperMixin<scalar,time>*>(DST);
	auto* src = dynamic_cast<const EncapsulatedSweeperMixin<scalar,time>*>(SRC);

	if (tmat.size1() == 0)
	  tmat = pfasst::compute_interp<time>(dst->get_nodes(), src->get_nodes());

	int ndst = dst->get_nodes().size();
	int nsrc = src->get_nodes().size();

	auto* crse_factory = src->get_factory();
	auto* fine_factory = dst->get_factory();

	vector<Encapsulation<scalar,time>*> fine_q(ndst), fine_tmp(nsrc);

	for (int m=0; m<ndst; m++) fine_q[m]   = dst->get_q(m);
	for (int m=0; m<nsrc; m++) fine_tmp[m] = fine_factory->create(solution);

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

	for (int m=0; m<nsrc; m++) delete fine_tmp[m];
	for (int m=0; m<ndst; m++) dst->evaluate(m);
      }

      virtual void restrict(ISweeper *DST, const ISweeper *SRC)
      {
	auto* dst = dynamic_cast<EncapsulatedSweeperMixin<scalar,time>*>(DST);
	auto* src = dynamic_cast<const EncapsulatedSweeperMixin<scalar,time>*>(SRC);

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

      virtual void fas(ISweeper *DST, const ISweeper *SRC)
      {
	auto* dst = dynamic_cast<EncapsulatedSweeperMixin<scalar,time>*>(DST);
	auto* src = dynamic_cast<const EncapsulatedSweeperMixin<scalar,time>*>(SRC);

	int ndst = dst->get_nodes().size();
	int nsrc = src->get_nodes().size();

	auto* crse_factory = src->get_factory();
	auto* fine_factory = dst->get_factory();

	vector<Encapsulation<scalar,time>*> fine_q(ndst), fine_tmp(nsrc);

	for (int m=0; m<ndst; m++) fine_q[m]   = dst->get_q(m);
	for (int m=0; m<nsrc; m++) fine_tmp[m] = fine_factory->create(solution);

      }

      // required for interp/restrict helpers
      virtual void interpolate(Encapsulation<scalar,time> *dst, const Encapsulation<scalar,time> *src) {
	throw NotImplementedYet("mlsdc/pfasst");
      }

      virtual void restrict(Encapsulation<scalar,time> *dst, const Encapsulation<scalar,time> *src) {
	throw NotImplementedYet("mlsdc/pfasst");
      }

    };

  }
}

#endif
