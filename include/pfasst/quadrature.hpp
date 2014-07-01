
#ifndef _PFASST_QUADRATURE_HPP_
#define _PFASST_QUADRATURE_HPP_

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using std::complex;
using std::string;
using std::vector;

using boost::numeric::ublas::matrix;

namespace pfasst
{

  typedef unsigned int uint;

  template<typename coeffT>
  class Polynomial
  {
      vector<coeffT> c;

    public:

      Polynomial(uint n) : c(n)
      {
        fill(c.begin(), c.end(), 0.0);
      }

      uint order() const
      {
        return c.size() - 1;
      }

      coeffT& operator[](const unsigned int i) { return c.at(i); }

      Polynomial<coeffT> differentiate() const
      {
        Polynomial<coeffT> p(c.size() - 1);

        for (int j = 1; j < c.size(); j++)
        { p[j - 1] = j * c[j]; }

        return p;
      }

      Polynomial<coeffT> integrate() const
      {
        Polynomial<coeffT> p(c.size() + 1);

        for (int j = 0; j < c.size(); j++)
        { p[j + 1] = c[j] / (j + 1); }

        return p;
      }

      template<typename xtype>
      xtype evaluate(const xtype x) const
      {
        int   n = c.size() - 1;
        xtype v = c[n];

        for (int j = n - 1; j >= 0; j--)
        { v = x * v + c[j]; }

        return v;
      }

      Polynomial<coeffT> normalize() const
      {
        Polynomial<coeffT> p(c.size());

        for (int j = 0; j < c.size(); j++)
        { p[j] = c[j] / c[c.size() - 1]; }

        return p;
      }

      vector<coeffT> roots() const
      {
        uint n = c.size() - 1;

        // initial guess
        Polynomial<complex<coeffT>> z0(n), z1(n);

        for (int j = 0; j < n; j++) {
          z0[j] = pow(complex<double>(0.4, 0.9), j);
          z1[j] = z0[j];
        }

        // durand-kerner-weierstrass iterations
        Polynomial<coeffT> p = normalize();

        for (int k = 0; k < 100; k++) {
          complex<coeffT> num, den;

          for (int i = 0; i < n; i++) {
            num = p.evaluate(z0[i]);
            den = 1.0;

            for (int j = 0; j < n; j++) {
              if (j == i) { continue; }

              den = den * (z0[i] - z0[j]);
            }

            z0[i] = z0[i] - num / den;
          }

          // converged?
          coeffT acc = 0.0;

          for (int j = 0; j < n; j++)
          { acc += abs(z0[j] - z1[j]); }

          if (acc < 2 * std::numeric_limits<coeffT>::epsilon())
          { break; }

          z1 = z0;
        }

        vector<coeffT> roots(n);

        for (int j = 0; j < n; j++)
        { roots[j] = abs(z0[j]) < 4 * std::numeric_limits<coeffT>::epsilon() ? 0.0 : real(z0[j]); }

        sort(roots.begin(), roots.end());
        return roots;
      }

      static Polynomial<coeffT> legendre(const uint order)
      {
        if (order == 0) {
          Polynomial<coeffT> p(1);
          p[0] = 1.0;
          return p;
        }

        if (order == 1) {
          Polynomial<coeffT> p(2);
          p[0] = 0.0;
          p[1] = 1.0;
          return p;
        }

        Polynomial<coeffT> p0(order + 1), p1(order + 1), p2(order + 1);
        p0[0] = 1.0; p1[1] = 1.0;

        // (n + 1) P_{n+1} = (2n + 1) x P_{n} - n P_{n-1}
        for (int m = 1; m < order; m++) {
          for (int j = 1; j < order + 1; j++)
          { p2[j] = ((2 * m + 1) * p1[j - 1] - m * p0[j]) / (m + 1); }

          p2[0] = - m * p0[0] / (m + 1);

          for (int j = 0; j < order + 1; j++) {
            p0[j] = p1[j];
            p1[j] = p2[j];
          }
        }

        return p2;
      }
  };

  //#define pi 3.1415926535897932384626433832795028841971693993751

  template<typename timeT>
  vector<timeT> compute_nodes(int nnodes, string qtype)
  {
    vector<timeT> nodes(nnodes);

    if (qtype == "gauss-legendre") {
      auto roots = Polynomial<timeT>::legendre(nnodes).roots();

      for (int j = 0; j < nnodes; j++)
      { nodes[j] = 0.5 * (1.0 + roots[j]); }
    } else if (qtype == "gauss-lobatto") {
      auto roots = Polynomial<timeT>::legendre(nnodes - 1).differentiate().roots();

      for (int j = 0; j < nnodes - 2; j++)
      { nodes[j + 1] = 0.5 * (1.0 + roots[j]); }

      nodes[0] = 0.0; nodes[nnodes - 1] = 1.0;
    } else if (qtype == "gauss-radau") {
      auto l   = Polynomial<timeT>::legendre(nnodes);
      auto lm1 = Polynomial<timeT>::legendre(nnodes - 1);

      for (int i = 0; i < nnodes; i++)
      { l[i] += lm1[i]; }

      auto roots = l.roots();

      for (int j = 1; j < nnodes; j++)
      { nodes[j - 1] = 0.5 * (1.0 - roots[nnodes - j]); }

      nodes[nnodes - 1] = 1.0;
    }

    return nodes;
  }

  template<typename timeT>
  matrix<timeT> compute_quadrature(vector<timeT> dst, vector<timeT> src, char type)
  {
    const int ndst = dst.size();
    const int nsrc = src.size();

    matrix<timeT> mat(ndst - 1, nsrc);

    //   /* for (int n=0; n<(ndst-1)*nsrc; n++) */
    //   /*   smat[n] = 0.0; */

    Polynomial<timeT> p(nsrc + 1), p1(nsrc + 1);

    for (int i = 0; i < nsrc; i++) {
      //      if ((flags[i] & SDC_NODE_PROPER) == 0) continue;

      // construct interpolating polynomial coefficients
      p[0] = 1.0;

      for (int j = 1; j < nsrc + 1; j++) { p[j] = 0.0; }

      for (int m = 0; m < nsrc; m++) {
        //if (((flags[m] & SDC_NODE_PROPER) == 0) || (m == i)) continue;
        if (m == i) { continue; }

        // p_{m+1}(x) = (x - x_j) * p_m(x)
        p1[0] = 0.0;

        for (int j = 0; j < nsrc;   j++) { p1[j + 1]  = p[j]; }

        for (int j = 0; j < nsrc + 1; j++) { p1[j]   -= p[j] * src[m]; }

        for (int j = 0; j < nsrc + 1; j++) { p[j] = p1[j]; }
      }

      // evaluate integrals
      auto den = p.evaluate(src[i]);
      auto P = p.integrate();

      for (int j = 1; j < ndst; j++) {
        timeT q = 0.0;

        if (type == 's')
        { q = P.evaluate(dst[j]) - P.evaluate(dst[j - 1]); }
        else
        { q = P.evaluate(dst[j]) - P.evaluate(0.0); }

        mat(j - 1, i) = q / den;
      }
    }

    return mat;
  }

  template<typename timeT>
  matrix<timeT> compute_interp(vector<timeT> dst, vector<timeT> src)
  {
    const int ndst = dst.size();
    const int nsrc = src.size();

    matrix<timeT> mat(ndst, nsrc);

    for (int i = 0; i < ndst; i++) {
      for (int j = 0; j < nsrc; j++) {
        timeT den = 1.0;
        timeT num = 1.0;

        for (int k = 0; k < nsrc; k++) {
          if (k == j) { continue; }

          den *= src[j] - src[k];
          num *= dst[i] - src[k];
        }

        if (abs(num) > 1e-32)
        { mat(i, j) = num / den; }
      }
    }

    return mat;
  }

}

#endif
