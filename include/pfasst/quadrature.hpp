
#ifndef _PFASST_QUADRATURE_HPP_
#define _PFASST_QUADRATURE_HPP_

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using std::complex;
using std::string;
using std::vector;

using boost::numeric::ublas::matrix;

#include "interfaces.hpp"

#define PI 3.1415926535897932384626433832795028841971693993751

namespace pfasst
{
  template<typename CoeffT>
  class Polynomial
  {
      vector<CoeffT> c;

    public:
      Polynomial(size_t n) : c(n)
      {
        fill(c.begin(), c.end(), 0.0);
      }

      size_t order() const
      {
        return c.size() - 1;
      }

      CoeffT& operator[](const size_t i)
      {
        return c.at(i);
      }

      Polynomial<CoeffT> differentiate() const
      {
        Polynomial<CoeffT> p(c.size() - 1);
        for (size_t j = 1; j < c.size(); j++) {
          p[j - 1] = j * c[j];
        }
        return p;
      }

      Polynomial<CoeffT> integrate() const
      {
        Polynomial<CoeffT> p(c.size() + 1);
        for (size_t j = 0; j < c.size(); j++) {
          p[j + 1] = c[j] / (j + 1);
        }
        return p;
      }

      template<typename xtype>
      xtype evaluate(const xtype x) const
      {
        int n = c.size() - 1;
        xtype v = c[n];
        for (int j = n - 1; j >= 0; j--) {
          v = x * v + c[j];
        }
        return v;
      }

      Polynomial<CoeffT> normalize() const
      {
        Polynomial<CoeffT> p(c.size());
        for (size_t j = 0; j < c.size(); j++) {
          p[j] = c[j] / c.back();
        }
        return p;
      }

      vector<CoeffT> roots() const
      {
        assert(c.size() >= 1);
        size_t n = c.size() - 1;

        // initial guess
        Polynomial<complex<CoeffT>> z0(n), z1(n);
        for (size_t j = 0; j < n; j++) {
          z0[j] = pow(complex<double>(0.4, 0.9), j);
          z1[j] = z0[j];
        }

        // durand-kerner-weierstrass iterations
        Polynomial<CoeffT> p = normalize();
        for (size_t k = 0; k < 100; k++) {
          complex<CoeffT> num, den;
          for (size_t i = 0; i < n; i++) {
            num = p.evaluate(z0[i]);
            den = 1.0;
            for (size_t j = 0; j < n; j++) {
              if (j == i) { continue; }
              den = den * (z0[i] - z0[j]);
            }
            z0[i] = z0[i] - num / den;
          }

          // converged?
          CoeffT acc = 0.0;
          for (size_t j = 0; j < n; j++) { acc += abs(z0[j] - z1[j]); }
          if (acc < 2 * std::numeric_limits<CoeffT>::epsilon()) { break; }

          z1 = z0;
        }

        vector<CoeffT> roots(n);
        for (size_t j = 0; j < n; j++) {
          roots[j] = (abs(z0[j]) < 4 * std::numeric_limits<CoeffT>::epsilon()) ? 0.0 : real(z0[j]);
        }

        sort(roots.begin(), roots.end());
        return roots;
      }

      static Polynomial<CoeffT> legendre(const size_t order)
      {
        if (order == 0) {
          Polynomial<CoeffT> p(1);
          p[0] = 1.0;
          return p;
        }

        if (order == 1) {
          Polynomial<CoeffT> p(2);
          p[0] = 0.0;
          p[1] = 1.0;
          return p;
        }

        Polynomial<CoeffT> p0(order + 1), p1(order + 1), p2(order + 1);
        p0[0] = 1.0; p1[1] = 1.0;

        // (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
        for (size_t m = 1; m < order; m++) {
          for (size_t j = 1; j < order + 1; j++) {
            p2[j] = ((2 * m + 1) * p1[j - 1] - m * p0[j]) / (m + 1);
          }
          p2[0] = - int(m) * p0[0] / (m + 1);

          for (size_t j = 0; j < order + 1; j++) {
            p0[j] = p1[j];
            p1[j] = p2[j];
          }
        }

        return p2;
      }
  };

  template<typename node = time_precision>
  vector<node> compute_nodes(size_t nnodes, string qtype)
  {
    vector<node> nodes(nnodes);

    if (qtype == "gauss-legendre") {

      auto roots = Polynomial<node>::legendre(nnodes).roots();
      for (size_t j = 0; j < nnodes; j++) {
        nodes[j] = 0.5 * (1.0 + roots[j]);
      }

    } else if (qtype == "gauss-lobatto") {

      auto roots = Polynomial<node>::legendre(nnodes - 1).differentiate().roots();
      assert(nnodes >= 2);
      for (size_t j = 0; j < nnodes - 2; j++) {
        nodes[j + 1] = 0.5 * (1.0 + roots[j]);
      }
      nodes.front() = 0.0;
      nodes.back() = 1.0;

    } else if (qtype == "gauss-radau") {

      auto l   = Polynomial<node>::legendre(nnodes);
      auto lm1 = Polynomial<node>::legendre(nnodes - 1);
      for (size_t i = 0; i < nnodes; i++) {
        l[i] += lm1[i];
      }
      auto roots = l.roots();
      for (size_t j = 1; j < nnodes; j++) {
        nodes[j - 1] = 0.5 * (1.0 - roots[nnodes - j]);
      }
      nodes.back() = 1.0;

    } else if (qtype == "clenshaw-curtis") {

      for (size_t j = 0; j < nnodes; j++) {
        nodes[j] = 0.5 * (1.0 - cos(j * PI / (nnodes-1)));
      }

    } else if (qtype == "uniform") {

      for (size_t j = 0; j < nnodes; j++) {
        nodes[j] = node(j) / (nnodes-1);
      }

    } else {

      throw ValueError("invalid node type passed to compute_nodes.");

    }

    return nodes;
  }

  template<typename node = time_precision>
  matrix<node> compute_quadrature(vector<node> dst, vector<node> src, char type)
  {
    const size_t ndst = dst.size();
    const size_t nsrc = src.size();

    assert(ndst >= 1);
    matrix<node> mat(ndst - 1, nsrc);

    Polynomial<node> p(nsrc + 1), p1(nsrc + 1);

    for (size_t i = 0; i < nsrc; i++) {

      // construct interpolating polynomial coefficients
      p[0] = 1.0;
      for (size_t j = 1; j < nsrc + 1; j++) { p[j] = 0.0; }
      for (size_t m = 0; m < nsrc; m++) {
        if (m == i) { continue; }

        // p_{m+1}(x) = (x - x_j) * p_m(x)
        p1[0] = 0.0;
        for (size_t j = 0; j < nsrc;   j++) { p1[j + 1]  = p[j]; }
        for (size_t j = 0; j < nsrc + 1; j++) { p1[j]   -= p[j] * src[m]; }
        for (size_t j = 0; j < nsrc + 1; j++) { p[j] = p1[j]; }
      }

      // evaluate integrals
      auto den = p.evaluate(src[i]);
      auto P = p.integrate();
      for (size_t j = 1; j < ndst; j++) {
        node q = 0.0;
        if (type == 's') {
          q = P.evaluate(dst[j]) - P.evaluate(dst[j - 1]);
        } else {
          q = P.evaluate(dst[j]) - P.evaluate(0.0);
        }

        mat(j - 1, i) = q / den;
      }
    }

    return mat;
  }

  template<typename node = time_precision>
  matrix<node> compute_interp(vector<node> dst, vector<node> src)
  {
    const size_t ndst = dst.size();
    const size_t nsrc = src.size();

    matrix<node> mat(ndst, nsrc);

    for (size_t i = 0; i < ndst; i++) {
      for (size_t j = 0; j < nsrc; j++) {
        node den = 1.0;
        node num = 1.0;

        for (size_t k = 0; k < nsrc; k++) {
          if (k == j) { continue; }
          den *= src[j] - src[k];
          num *= dst[i] - src[k];
        }

        if (abs(num) > 1e-32) {
          mat(i, j) = num / den;
        } else {
	  mat(i, j) = 0.0;
	}
      }
    }

    return mat;
  }

}

#endif
