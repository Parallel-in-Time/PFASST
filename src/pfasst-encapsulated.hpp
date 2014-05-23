/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_ENCAPSULATED_HPP_
#define _PFASST_ENCAPSULATED_HPP_

#include <vector>

#include "pfasst.hpp"

using namespace std;

namespace pfasst {

  template<typename T>
  vector<T> compute_nodes(unsigned int nnodes, string qtype) {
    vector<T> nodes(nnodes);

    // ...

    return nodes;
  }

  typedef enum encaptype { solution, function } encaptype;

  template<typename T>
  class matrix : public vector<T> {

  public:
    unsigned int n, m;
    matrix() { }
    matrix(unsigned int n, unsigned int m) {
      zeros(n, m);
    }
    void zeros(unsigned int n, unsigned int m) {
      this->n = n; this->m = m;
      this->resize(n*m);
      // ...
    }
    T& operator()(unsigned int i, unsigned int j) {
      return (*this)[i*m+j];
    }
  };

  //
  // encapsulation
  //

  struct encapsulation {
    virtual ~encapsulation() { }

    // required for time-parallel communications
    virtual unsigned int nbytes() { }
    virtual void pack(char *buf) { }
    virtual void unpack(char *buf) { }

    // required for interp/restrict helpers
    virtual void interpolate(const encapsulation *) { }
    virtual void restrict(const encapsulation *) { }

    // required for host based encap helpers
    virtual void setval(double) { }
    virtual void copy(const encapsulation *) { }
    virtual void mat_apply(encapsulation dst[], double a, matrix m, const encapsulation src[]) { }
  };

  struct encapsulation_factory {
    virtual encapsulation* create(const encaptype) = 0;
  };


  template<typename T>
  class encapsulated_sweeper_mixin : public isweeper {
    shared_ptr<vector<T>>             nodes;
    shared_ptr<encapsulation_factory> encap;

  public:
    vector<encapsulation*> q;
    vector<T>* get_nodes() { return nodes.get(); }

    virtual void set_q0(const encapsulation* q0) { }
    virtual encapsulation* get_qend() { }
  };

  template<class T>
  class poly_interp_mixin : public T {
    virtual void interpolate(const isweeper*) { }
    virtual void restrict(const isweeper*) { }
  };

  template<typename T>
  struct vector_encapsulation : public vector<T>, public encapsulation {
    vector_encapsulation(int size) : vector<T>(size) { }
    virtual unsigned int nbytes() const {
      return sizeof(T) * this->size();
    }
    void setval(double v) {
      for (int i=0; i<this->size(); i++)
	this->data()[i] = v;
    }
    // ...
  };

  template<typename T>
  class vector_factory : public pfasst::encapsulation_factory {
    int size;
  public:
    vector_factory(const int size) : size(size) { }
    encapsulation* create(const pfasst::encap_type) {
      return new vector_encapsulation<T>(size);
    }
  };

}

#endif
