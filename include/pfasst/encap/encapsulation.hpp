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
    class EncapFactory {
    public:
      virtual Encapsulation<scalar,time>* create(const EncapType) = 0;
    };

    template<typename scalar, typename time>
    class EncapSweeper : public ISweeper {
      vector<scalar> nodes;
      shared_ptr<EncapFactory<scalar,time>> factory;

    public:

      void set_nodes(vector<time> nodes)
      {
	this->nodes = nodes;
      }

      const vector<time> get_nodes() const
      {
	return nodes;
      }

      void set_factory(EncapFactory<scalar,time>* factory)
      {
	this->factory = shared_ptr<EncapFactory<scalar,time>>(factory);
      }

      EncapFactory<scalar,time>* get_factory() const
      {
	return factory.get();
      }

      virtual void set_state(const Encapsulation<scalar,time>* q0, unsigned int m)
      {
      	throw NotImplementedYet("sweeper");
      }

      virtual Encapsulation<scalar,time>* get_state(unsigned int m) const
      {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar,time>* get_tau(unsigned int m) const
      {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar,time>* get_saved_state(unsigned int m) const
      {
	throw NotImplementedYet("sweeper");
	return NULL;
      }

      virtual Encapsulation<scalar,time>* get_end_state()
      {
	return this->get_state(this->get_nodes().size()-1);
      }

      virtual void evaluate(int m)
      {
	throw NotImplementedYet("sweeper");
      }

      virtual void advance()
      {
	this->set_state(this->get_end_state(), 0);
      }

      virtual void integrate(time dt, vector<Encapsulation<scalar,time>*> dst) const
      {
	throw NotImplementedYet("sweeper");
      }
    };

    template<typename scalar, typename time>
    class PolyInterpMixin : public pfasst::ITransfer {
      matrix<time> tmat, fmat;

    public:

      virtual ~PolyInterpMixin() { }

      virtual void interpolate(ISweeper *DST, const ISweeper *SRC, bool initial)
      {
	auto* dst = dynamic_cast<EncapSweeper<scalar,time>*>(DST);
	auto* src = dynamic_cast<const EncapSweeper<scalar,time>*>(SRC);

	if (tmat.size1() == 0)
	  tmat = pfasst::compute_interp<time>(dst->get_nodes(), src->get_nodes());

	int ndst = dst->get_nodes().size();
	int nsrc = src->get_nodes().size();

	auto* crse_factory = src->get_factory();
	auto* fine_factory = dst->get_factory();

	vector<Encapsulation<scalar,time>*> fine_q(ndst), fine_tmp(nsrc);

	for (int m=0; m<ndst; m++) fine_q[m]   = dst->get_state(m);
	for (int m=0; m<nsrc; m++) fine_tmp[m] = fine_factory->create(solution);

	if (initial)
	  for (int m=1; m<ndst; m++)
	    fine_q[m]->copy(fine_q[0]);

	auto* crse_tmp = crse_factory->create(solution);
	for (int m=0; m<nsrc; m++) {
	  crse_tmp->copy(src->get_state(m));
	  if (initial)
	    crse_tmp->saxpy(-1.0, src->get_state(0));
	  else
	    crse_tmp->saxpy(-1.0, src->get_saved_state(m));
	  interpolate(fine_tmp[m], crse_tmp);
	}
	delete crse_tmp;

	dst->get_state(0)->mat_apply(fine_q, 1.0, tmat, fine_tmp, false);

	for (int m=0; m<nsrc; m++) delete fine_tmp[m];
	for (int m=0; m<ndst; m++) dst->evaluate(m);
      }

      virtual void restrict(ISweeper *DST, const ISweeper *SRC)
      {
	auto* dst = dynamic_cast<EncapSweeper<scalar,time>*>(DST);
	auto* src = dynamic_cast<const EncapSweeper<scalar,time>*>(SRC);

	auto dnodes = dst->get_nodes();
	auto snodes = src->get_nodes();

	int ndst = dst->get_nodes().size();
	int nsrc = src->get_nodes().size();

	int trat = (nsrc - 1) / (ndst - 1);

	for (int m=0; m<ndst; m++) {
	  if (dnodes[m] != snodes[m*trat])
	    throw NotImplementedYet("coarse nodes must be nested");
	  this->restrict(dst->get_state(m), src->get_state(m*trat));
	}

	for (int m=0; m<ndst; m++) dst->evaluate(m);
      }

      virtual void fas(time dt, ISweeper *dst, const ISweeper *src)
      {
	auto* crse = dynamic_cast<EncapSweeper<scalar,time>*>(dst);
	auto* fine = dynamic_cast<const EncapSweeper<scalar,time>*>(src);

	int ncrse = crse->get_nodes().size();
	int nfine = fine->get_nodes().size();

	auto* crse_factory = crse->get_factory();
	auto* fine_factory = fine->get_factory();

	vector<Encapsulation<scalar,time>*> crse_z2n(ncrse-1), fine_z2n(nfine-1), rstr_z2n(ncrse-1);
	for (int m=0; m<ncrse-1; m++) crse_z2n[m] = crse_factory->create(solution);
	for (int m=0; m<ncrse-1; m++) rstr_z2n[m] = crse_factory->create(solution);
	for (int m=0; m<nfine-1; m++) fine_z2n[m] = fine_factory->create(solution);

	// compute '0 to node' integral on the coarse level
	crse->integrate(dt, crse_z2n);
	for (int m=1; m<ncrse-1; m++)
	  crse_z2n[m]->saxpy(1.0, crse_z2n[m-1]);

	// compute '0 to node' integral on the fine level
	fine->integrate(dt, fine_z2n);
	for (int m=1; m<nfine-1; m++)
	  fine_z2n[m]->saxpy(1.0, fine_z2n[m-1]);

	// restrict '0 to node' fine integral
	int trat = (nfine - 1) / (ncrse - 1);
	for (int m=1; m<ncrse; m++)
	  this->restrict(rstr_z2n[m-1], fine_z2n[m*trat-1]);

	// compute 'node to node' tau correction
	vector<Encapsulation<scalar,time>*> tau(ncrse-1);
	for (int m=0; m<ncrse-1; m++) tau[m] = crse->get_tau(m);

	tau[0]->copy(rstr_z2n[0]);
	tau[0]->saxpy(-1.0, crse_z2n[0]);

	for (int m=1; m<ncrse-1; m++) {
	  tau[m]->copy(rstr_z2n[m]);
	  tau[m]->saxpy(-1.0, rstr_z2n[m-1]);

	  tau[m]->saxpy(-1.0, crse_z2n[m]);
	  tau[m]->saxpy(1.0, crse_z2n[m-1]);
	}

	for (int m=0; m<ncrse-1; m++) delete crse_z2n[m];
	for (int m=0; m<ncrse-1; m++) delete rstr_z2n[m];
	for (int m=0; m<nfine-1; m++) delete fine_z2n[m];
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
