#ifndef _PFASST__QUADRATURE__INTERFACE_HPP_
#define _PFASST__QUADRATURE__INTERFACE_HPP_

#include <cassert>
#include <vector>

#include <Eigen/Dense>

#include "../globals.hpp"
#include "../interfaces.hpp"
#include "polynomial.hpp"

template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename scalar>
using Index = typename Matrix<scalar>::Index;

using namespace std;


namespace pfasst
{
  namespace quadrature
  {
    enum class QuadratureType : int {
        GaussLegendre   =  0
      , GaussLobatto    =  1
      , GaussRadau      =  2
      , ClenshawCurtis  =  3
      , Uniform         =  4
      , UNDEFINED       = -1
    };


    template<typename scalar>
    static Polynomial<scalar> build_polynomial(const size_t node, const vector<scalar>& nodes)
    {
      const size_t num_nodes = nodes.size();
      Polynomial<scalar> p(num_nodes + 1), p1(num_nodes + 1);
      p[0] = 1.0;

      for (size_t m = 0; m < num_nodes; ++m) {
        if (m == node) { continue; }

        // p_{m+1}(x) = (x - x_j) * p_m(x)
        p1[0] = scalar(0.0);
        for (size_t j = 0; j < num_nodes;     ++j) { p1[j + 1]  = p[j]; }
        for (size_t j = 0; j < num_nodes + 1; ++j) { p1[j]     -= p[j] * nodes[m]; }
        for (size_t j = 0; j < num_nodes + 1; ++j) { p[j]       = p1[j]; }
      }

      return p;
    }


    template<typename scalar>
    static Matrix<scalar> compute_q_matrix(const vector<scalar>& from, const vector<scalar>& to)
    {
      const size_t to_size = to.size();
      const size_t from_size = from.size();
      assert(to_size >= 1 && from_size >= 1);

      Matrix<scalar> q_mat = Matrix<scalar>::Zero(to_size, from_size);

      for (size_t m = 0; m < from_size; ++m) {
        Polynomial<scalar> p = build_polynomial(m, from);
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
        Polynomial<scalar> p = build_polynomial(m, nodes);
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
        static const bool LEFT_IS_NODE = false;
        static const bool RIGHT_IS_NODE = false;

        size_t num_nodes;
        Matrix<precision> q_mat;
        Matrix<precision> s_mat;
        vector<precision> q_vec; // XXX: black spot
        Matrix<precision> b_mat;
        vector<precision> nodes;
        //! @}

      public:
        //! @{
        explicit IQuadrature(const size_t num_nodes)
          : num_nodes(num_nodes)
        {
          if (this->num_nodes == 0) {
            throw invalid_argument("Any quadrature requires at least one quadrature nodes.");
          }
        }

        IQuadrature()
          : num_nodes(0)
        {}

        virtual ~IQuadrature() = default;
        //! @}

        //! @{
        virtual const Matrix<precision>& get_q_mat() const { return this->q_mat; }

        virtual const Matrix<precision>& get_s_mat() const { return this->s_mat; }

        virtual const Matrix<precision>& get_b_mat() const { return this->b_mat; }

        virtual const vector<precision>& get_q_vec() const { return this->q_vec; }

        virtual const vector<precision>& get_nodes() const { return this->nodes; }

        virtual size_t get_num_nodes() const { return this->num_nodes; }

        virtual bool left_is_node() const
        {
          throw NotImplementedYet("Quadrature");
          return LEFT_IS_NODE;
        }

        virtual bool right_is_node() const
        {
          throw NotImplementedYet("Quadrature");
          return RIGHT_IS_NODE;
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
          this->q_mat = compute_q_matrix(this->nodes);
          this->s_mat = compute_s_matrix(this->q_mat);
          this->q_vec = compute_q_vec(this->nodes);
          this->b_mat = Matrix<precision>::Zero(1, this->num_nodes);
          for (size_t i = 0; i < this->num_nodes; i++){
            this->b_mat(0,i) = this->q_vec[i];
          }
        }
        //! @}
    };

  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE__INTERFACE_HPP_
