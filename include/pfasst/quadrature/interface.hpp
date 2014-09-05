#ifndef _PFASST__QUADRATURE__INTERFACE_HPP_
#define _PFASST__QUADRATURE__INTERFACE_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "../globals.hpp"
#include "../interfaces.hpp"
#include "polynomial.hpp"
#include "traits.hpp"

template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename scalar>
using Index = typename Matrix<scalar>::Index;

using namespace std;


namespace pfasst
{
  namespace quadrature
  {
    template<typename scalar>
    static Polynomial<scalar> build_polynomial(const size_t node, const vector<scalar>& from)
    {
      const size_t from_size = from.size();
      Polynomial<scalar> p(from_size + 1), p1(from_size + 1);
      p[0] = 1.0;

      for (size_t m = 0; m < from_size; ++m) {
        if (m == node) { continue; }

        // p_{m+1}(x) = (x - x_j) * p_m(x)
        p1[0] = scalar(0.0);
        for (size_t j = 0; j < from_size;     ++j) { p1[j + 1]  = p[j]; }
        for (size_t j = 0; j < from_size + 1; ++j) { p1[j]     -= p[j] * from[m]; }
        for (size_t j = 0; j < from_size + 1; ++j) { p[j]       = p1[j]; }
      }

      return p;
    }


    template<typename scalar>
    static Polynomial<scalar> build_polynomial(const size_t node, const vector<scalar>& nodes)
    {
      return build_polynomial(node, nodes, nodes);
    }


    template<typename scalar>
    static Matrix<scalar> compute_q_matrix(const vector<scalar>& from, const vector<scalar>& to)
    {
      const size_t to_size = to.size();
      const size_t from_size = from.size();
      assert(to_size >= 1 && from_size >= 1);

      Matrix<scalar> q_mat = Matrix<scalar>::Zero(to_size, from_size);

      for (size_t m = 0; m < from_size; ++m) {
        Polynomial<scalar> p = build_polynomial(m, from, to);
        // evaluate integrals
        auto den = p.evaluate(from[m]);
        auto P = p.integrate();
        for (size_t j = 0; j < to_size; ++j) {
          q_mat(j, m) = (P.evaluate(to[j]) - P.evaluate(0.0)) / den;
        }
      }

      return q_mat;
    }


    template<typename scalar>
    static Matrix<scalar> compute_q_matrix(const vector<scalar>& nodes)
    {
      return compute_q_matrix<scalar>(nodes, nodes);
    }


    template<typename scalar>
    static Matrix<scalar> compute_q_matrix(const Matrix<scalar>& s_mat)
    {
      Matrix<scalar> q_mat = Matrix<scalar>::Zero(s_mat.rows(), s_mat.cols());
      q_mat.col(0) = s_mat.col(0);
      for (Index<scalar> q_mat_col = 1; q_mat_col < q_mat.cols(); ++q_mat_col) {
        q_mat.col(q_mat_col) = q_mat.col(q_mat_col - 1) + s_mat.col(q_mat_col);
      }
      return q_mat;
    }


    template<typename scalar>
    static Matrix<scalar> compute_s_matrix(const Matrix<scalar>& q_mat)
    {
      Matrix<scalar> s_mat = Matrix<scalar>::Zero(q_mat.rows(), q_mat.cols());
      s_mat.row(0) = q_mat.row(0);
      for (Index<scalar> row = 1; row < s_mat.rows(); ++row) {
        s_mat.row(row) = q_mat.row(row) - q_mat.row(row - 1);
      }
      return s_mat;
    }


    template<typename scalar>
    static Matrix<scalar> compute_s_matrix(const vector<scalar>& from, const vector<scalar>& to)
    {
      return compute_s_matrix(compute_q_matrix(from, to));
    }


    template<typename scalar>
    static vector<scalar> compute_q_vec(const vector<scalar>& nodes)
    {
      const size_t num_nodes = nodes.size();
      assert(num_nodes >= 1);

      vector<scalar> q_vec = vector<scalar>(num_nodes, scalar(0.0));

      for (size_t m = 0; m < num_nodes; ++m) {
        Polynomial<scalar> p = build_polynomial(m, nodes, nodes);
        // evaluate integrals
        auto den = p.evaluate(nodes[m]);
        auto P = p.integrate();
        q_vec[m] = (P.evaluate(scalar(1.0)) - P.evaluate(scalar(0.0))) / den;
      }

      return q_vec;
    }

    template<typename precision = pfasst::time_precision>
    class IQuadrature
    {
      protected:
        //! @{
        size_t m_num_nodes;
        Matrix<precision> m_q_mat;
        Matrix<precision> m_s_mat;
        vector<precision> m_q_vec;
        vector<precision> m_nodes;
        vector<precision> m_delta_nodes;
        //! @}

      public:
        //! @{
        IQuadrature(const size_t num_nodes)
          : m_num_nodes(num_nodes)
        {
          if (this->m_num_nodes == 0) {
            throw invalid_argument("Any quadrature requires at least one quadrature nodes.");
          }
        }

        IQuadrature()
          : m_num_nodes(0)
        {}

        IQuadrature(const IQuadrature<precision>& other)
          :   m_num_nodes(other.m_num_nodes)
            , m_q_mat(other.m_q_mat)
            , m_s_mat(other.m_s_mat)
            , m_q_vec(other.m_q_vec)
            , m_nodes(other.m_nodes)
            , m_delta_nodes(other.m_delta_nodes)
        {}

        IQuadrature(IQuadrature<precision>&& other)
          : IQuadrature<precision>()
        {
          swap(*this, other);
        }

        virtual ~IQuadrature()
        {}
        //! @}

        //! @{
        virtual const Matrix<precision>& q_mat() const
        { return this->m_q_mat; }

        virtual const Matrix<precision>& s_mat() const
        { return this->m_s_mat; }

        virtual const vector<precision>& q_vec() const
        { return this->m_q_vec; }

        virtual const vector<precision>& nodes() const
        { return this->m_nodes; }

        virtual const vector<precision>& delta_nodes() const
        { return this->m_delta_nodes; }

        virtual size_t num_nodes() const
        { return this->m_num_nodes; }

        virtual bool left_is_node() const
        { return quadrature_traits<IQuadrature<precision>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<IQuadrature<precision>>::right_is_node; }
        //! @}

        //! @{
        //! copy-and-nothrow-swap idiom
        IQuadrature<precision>& operator=(IQuadrature<precision> other)
        {
          // pass-by-value to implicitly call copy-constructor by compiler
          swap(*this, other);
          return *this;
        }

        //! copy-and-nothrow-swap idiom
        friend void swap(IQuadrature<precision>& first, IQuadrature<precision>& second) noexcept
        {
          using std::swap;
          swap(first.m_num_nodes, second.m_num_nodes);
          swap(first.m_q_mat, second.m_q_mat);
          swap(first.m_s_mat, second.m_s_mat);
          swap(first.m_q_vec, second.m_q_vec);
          swap(first.m_nodes, second.m_nodes);
          swap(first.m_delta_nodes, second.m_delta_nodes);
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes()
        {
          throw NotImplementedYet("Implemented by derived classes.");
        }

        virtual void compute_weights()
        {
          this->m_q_mat = compute_q_matrix(this->m_nodes);
          this->m_s_mat = compute_s_matrix(this->m_q_mat);
          this->m_q_vec = compute_q_vec(this->m_nodes);
        }

        virtual void compute_delta_nodes()
        {
          this->m_delta_nodes.resize(this->m_num_nodes);
          this->m_delta_nodes[0] = this->m_nodes[0] - precision(0.0);
          for (size_t m = 1; m < this->m_num_nodes; ++m) {
            this->m_delta_nodes[m] = this->m_nodes[m] - this->m_nodes[m - 1];
          }
        }
        //! @}
    };

  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__INTERFACE_HPP_
