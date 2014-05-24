
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;

typedef unsigned int uint;

namespace pfasst {

  template<typename ptype>
  class polynomial {
    vector<ptype> c;

  public:

    polynomial(uint n) : c(n) {
      fill(c.begin(), c.end(), 0.0);
    }

    uint order() const {
      return c.size()-1;
    }

    ptype& operator[](const unsigned int i) { return c.at(i); }

    polynomial<ptype> differentiate() const {
      polynomial<ptype> p(c.size()-1);
      for (int j=1; j<c.size(); j++)
	p[j-1] = j * c[j];
      return p;
    }

    polynomial<ptype> integrate() const {
      polynomial<ptype> p(c.size()+1);
      for (int j=0; j<c.size(); j++)
    	p[j+1] = c[j] / (j+1);
      return p;
    }

    template<typename xtype>
    xtype evaluate(const xtype x) const {
      int   n = c.size()-1;
      xtype v = c[n];
      for (int j=n-1; j>=0; j--)
    	v = x * v + c[j];
      return v;
    }

    polynomial<ptype> normalize() const {
      polynomial<ptype> p(c.size());
      for (int j=0; j<c.size(); j++)
	p[j] = c[j] / c[c.size()-1];
      return p;
    }

    vector<ptype> roots() const {
      uint n = c.size()-1;

      // initial guess
      polynomial<complex<ptype>> z0(n), z1(n);
      for (int j=0; j<n; j++) {
    	z0[j] = pow((0.4, 0.9), j);
    	z1[j] = z0[j];
      }

      // durand-kerner-weierstrass iterations
      polynomial<ptype> p = normalize();
      for (int k=0; k<100; k++) {
	complex<ptype> num, den;
    	for (int i=0; i<n; i++) {
    	  num = p.evaluate(z0[i]);
    	  den = 1.0;
    	  for (int j=0; j<n; j++) {
    	    if (j == i) continue;
    	    den = den * (z0[i] - z0[j]);
    	  }
    	  z0[i] = z0[i] - num / den;
    	}

    	// converged?
    	ptype acc = 0.0;
    	for (int j=0; j<n; j++)
    	  acc += abs(z0[j] - z1[j]);
    	if (acc < 1e-24)
	  break;

	z1 = z0;
      }

      vector<ptype> roots(n);
      for (int j=0; j<n; j++)
    	roots[j] = abs(z0[j]) < 1e-12 ? 0.0 : real(z0[j]);

      sort(roots.begin(), roots.end());
      return roots;
    }

    static polynomial<ptype> legendre(const uint order)
    {
      if (order == 0) {
	polynomial<ptype> p(1);
        p[0] = 1.0;
        return p;
      }

      if (order == 1) {
	polynomial<ptype> p(2);
        p[0] = 0.0;
        p[1] = 1.0;
        return p;
      }

      polynomial<ptype> p0(order+1), p1(order+1), p2(order+1);
      p0[0] = 1.0; p1[1] = 1.0;

      // (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
      for (int m=1; m<order; m++) {
        for (int j=1; j<order+1; j++)
	  p2[j] = ( (2*m + 1) * p1[j-1] - m * p0[j] ) / (m + 1);
        p2[0] = - m * p0[0] / (m + 1);

        for (int j=0; j<order+1; j++) {
	  p0[j] = p1[j];
	  p1[j] = p2[j];
        }
      }

      return p2;
    }
  };

  //#define pi 3.1415926535897932384626433832795028841971693993751

  template<typename ntype>
  vector<ntype> compute_nodes(int nnodes, string qtype)
  {
    vector<ntype> nodes(nnodes);

    if (qtype == "gauss-legendre") {
      auto roots = polynomial<ntype>::legendre(nnodes).roots();
      for (int j=0; j<nnodes; j++)
      	nodes[j] = 0.5 * (1.0 + roots[j]);
    } else if (qtype == "gauss-lobatto") {
      auto roots = polynomial<ntype>::legendre(nnodes-1).differentiate().roots();
      for (int j=0; j<nnodes-2; j++)
	nodes[j+1] = 0.5 * (1.0 + roots[j]);
      nodes[0] = 0.0; nodes[nnodes-1] = 1.0;
    }

    return nodes;
  }


// void sdc_smat(sdc_mat *smat,
// 	      const int n, const int m, const int sign,
// 	      const sdc_dtype *dst, const sdc_dtype *src,
// 	      const int *flags, int ndst, int nsrc)
// {
//   /* for (int n=0; n<(ndst-1)*nsrc; n++) */
//   /*   smat[n] = 0.0; */

//   sdc_dtype p[nsrc+1], p1[nsrc+1];
//   for (int i=0; i<nsrc; i++) {
//     if ((flags[i] & SDC_NODE_PROPER) == 0) continue;

//     // construct interpolating polynomial coefficients
//     p[0] = 1.0; for (int j=1; j<nsrc+1; j++) p[j] = 0.0;
//     for (int m=0; m<nsrc; m++) {
//       if (((flags[m] & SDC_NODE_PROPER) == 0) || (m == i)) continue;
//       // p_{m+1}(x) = (x - x_j) * p_m(x)
//       p1[0] = 0.0;
//       for (int j=0; j<nsrc;   j++) p1[j+1]  = p[j];
//       for (int j=0; j<nsrc+1; j++) p1[j]   -= p[j] * src[m];
//       for (int j=0; j<nsrc+1; j++) p[j] = p1[j];
//     }

//     // evaluate integrals
//     sdc_dtype den = poly_eval(p, nsrc, src[i]);
//     poly_int(p, nsrc+1);
//     for (int j=1; j<ndst; j++) {
//       sdc_dtype s = poly_eval(p, nsrc, dst[j]) - poly_eval(p, nsrc, dst[j-1]);
//       sdc_mat_setvalue(smat, n+j-1, m+i, s / den, sign);
//     }
//   }
// }



}
