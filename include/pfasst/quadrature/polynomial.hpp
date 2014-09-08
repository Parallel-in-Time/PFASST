#ifndef _PFASST__POLYNOMIAL_HPP_
#define _PFASST__POLYNOMIAL_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>


using namespace std;


namespace pfasst
{
  template<typename CoeffT>
  class Polynomial
  {
      vector<CoeffT> c;

    public:
      Polynomial(size_t n)
        : c(n, CoeffT(0.0))
      {}

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
          if (acc < 2 * numeric_limits<CoeffT>::epsilon()) { break; }

          z1 = z0;
        }

        vector<CoeffT> roots(n);
        for (size_t j = 0; j < n; j++) {
          roots[j] = (abs(z0[j]) < 4 * numeric_limits<CoeffT>::epsilon()) ? 0.0 : real(z0[j]);
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

}  // ::pfasst

#endif  // _PFASST__POLYNOMIAL_HPP_
