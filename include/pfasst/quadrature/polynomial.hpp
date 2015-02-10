#ifndef _PFASST__POLYNOMIAL_HPP_
#define _PFASST__POLYNOMIAL_HPP_

#include <vector>
using namespace std;


namespace pfasst
{
  template<typename CoeffT>
  class Polynomial
  {
    protected:
      vector<CoeffT> c;

    public:
      Polynomial(size_t n);

      size_t order() const;
      CoeffT& operator[](const size_t i);
      Polynomial<CoeffT> differentiate() const;
      Polynomial<CoeffT> integrate() const;

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

      Polynomial<CoeffT> normalize() const;
      vector<CoeffT> roots() const;
      static Polynomial<CoeffT> legendre(const size_t order);
  };
}  // ::pfasst

#include "pfasst/quadrature/polynomial_impl.hpp"

#endif  // _PFASST__POLYNOMIAL_HPP_
