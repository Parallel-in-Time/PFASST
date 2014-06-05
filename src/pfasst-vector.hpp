/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_VECTOR_HPP_
#define _PFASST_VECTOR_HPP_

#include <algorithm>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#define GPCHKERR(err, msg) if ((err) == -1) { perror(msg); return; }

#include "pfasst-encapsulated.hpp"

using namespace std;

namespace pfasst {
  namespace encap {

    template<typename scalar, typename time>
    class VectorEncapsulation : public vector<scalar>, public Encapsulation<scalar> {

      FILE *gp_to, *gp_fr;

    public:
      VectorEncapsulation(int size) : vector<scalar>(size) { gp_to = NULL; }

      void setval(scalar v) {
	std::fill(this->begin(), this->end(), v);
      }

      void copy(const Encapsulation<scalar>* X) {
	const auto* x = dynamic_cast<const VectorEncapsulation*>(X);
	std::copy(x->begin(), x->end(), this->begin());
      }

      void saxpy(time a, const Encapsulation<scalar> *X) {
	const auto& x = *dynamic_cast<const VectorEncapsulation*>(X);
	auto&       y = *this;

	for (int i=0; i<y.size(); i++)
	  y[i] += a * x[i];
      }

      void mat_apply(vector<Encapsulation<scalar>*> DST, time a, matrix<time> mat,
		     vector<Encapsulation<scalar>*> SRC, bool zero=true) {

	vector<VectorEncapsulation<scalar,time>*> dst(mat.n), src(mat.m);
	for (int n=0; n<mat.n; n++)
	  dst[n] = dynamic_cast<VectorEncapsulation*>(DST[n]);
	for (int m=0; m<mat.m; m++)
	  src[m] = dynamic_cast<VectorEncapsulation*>(SRC[m]);

	if (zero)
	  for (int n=0; n<mat.n; n++)
	    dst[n]->setval(0.0);

	for (int i=0; i<(*dst[0]).size(); i++)
	  for (int n=0; n<mat.n; n++)
	    for (int m=0; m<mat.m; m++)
	      (*dst[n])[i] += a * mat(n,m) * (*src[m])[i];
      }

      scalar norm0() const
      {
	scalar max = 0.0;
	for (int i=0; i<this->size(); i++) {
	  scalar v = abs((*this)[i]);
	  if (v > max)
	    max = v;
	}
	return max;
      }

#ifdef PFASST_ENABLE_GNUPLOT
      void plot(int window, bool wait)
      {
	if (gp_to == NULL) {

	  /*
	   * Fork a new gnuplot process and connect pipes.
	   */

	  int err, pid, to[2], fr[2];

	  err = pipe(to); GPCHKERR(err, "Unable to create pipe.");
	  err = pipe(fr); GPCHKERR(err, "Unable to create pipe.");
	  pid = fork();   GPCHKERR(pid, "Unable to fork gnuplot.");

	  if (pid) {
	    close(to[0]); close(fr[1]);
	  } else {
	    close(to[1]); close(fr[0]);
	    dup2(to[0], 0); close(to[0]);
	    dup2(fr[1], 1); close(fr[1]);
	    execl("/usr/bin/gnuplot", "gnuplot", NULL);
	  }

	  gp_to = fdopen(to[1], "w");
	  gp_fr = fdopen(fr[0], "r");

	}

	/*
	 * Send vector to gnuplot window.
	 */

	fprintf(gp_to, "set term wxt %d\n", window);
	fprintf(gp_to, "plot '-'\n");
	for (int i=0; i<this->size(); i++)
	  fprintf(gp_to, "%lg\n", (*this)[i]);
	fprintf(gp_to, "e\n");

	if (wait) {
	  fprintf(gp_to, "set print\n");
	  fprintf(gp_to, "pause mouse button2 \"===> paused - button 2 in window %d to release\"\n", window);
	  fprintf(gp_to, "print ''\n");
	  fprintf(gp_to, "set print '-'\n");
	  fprintf(gp_to, "print 'hoser!'\n");
	}
	fflush(gp_to);

	if (wait) {
	  char buf[8];
	  char *s = fgets(buf, 8, gp_fr);
	}
      }
#else
      void plot(int, bool) const { }
#endif

    };

    template<typename scalar, typename time>
    class VectorFactory : public EncapsulationFactory<scalar> {
      int size;
    public:
      int dofs() { return size; }
      VectorFactory(const int size) : size(size) { }
      Encapsulation<scalar>* create(const EncapType) {
	return new VectorEncapsulation<scalar,time>(size);
      }
    };

  }
}

#endif
