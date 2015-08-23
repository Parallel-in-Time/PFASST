/**
 * @file pfasst/quadrature/interface.hpp
 * @since v0.3.0
 */
#ifndef _PFASST__QUADRATURE__INTERFACE_HPP_
#define _PFASST__QUADRATURE__INTERFACE_HPP_

#include <cassert>
#include <vector>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <Eigen/Dense>
#include <leathers/pop>
template<typename precision>
using Matrix = Eigen::Matrix<precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename precision>
using Index = typename Matrix<precision>::Index;

#include "pfasst/globals.hpp"
#include "pfasst/exceptions.hpp"
#include "pfasst/quadrature/polynomial.hpp"


namespace pfasst
{
  namespace quadrature
  {
    /**
     * Quadrature type descriptors.
     * @since v0.3.0
     */
    enum class QuadratureType : int {
        GaussLegendre   =  0  //!< @ref pfasst::quadrature::GaussLegendre "Gauss-Legendre" quadrature
      , GaussLobatto    =  1  //!< @ref pfasst::quadrature::GaussLobatto "Gauss-Lobatto" quadrature
      , GaussRadau      =  2  //!< @ref pfasst::quadrature::GaussRadau "Gauss-Radau" quadrature
      , ClenshawCurtis  =  3  //!< @ref pfasst::quadrature::ClenshawCurtis "Clenshaw-Curtis" quadrature
      , Uniform         =  4  //!< Uniform quadrature
      , UNDEFINED       = -1
    };


    // TODO: get the grasp of this stuff
    template<typename precision>
    static Polynomial<precision> build_polynomial(const size_t node, const vector<precision>& nodes)
    {
      const size_t num_nodes = nodes.size();
      Polynomial<precision> p(num_nodes + 1), p1(num_nodes + 1);
      p[0] = 1.0;

      for (size_t m = 0; m < num_nodes; ++m) {
        if (m == node) { continue; }

        // p_{m+1}(x) = (x - x_j) * p_m(x)
        p1[0] = precision(0.0);
        for (size_t j = 0; j < num_nodes;     ++j) { p1[j + 1]  = p[j]; }
        for (size_t j = 0; j < num_nodes + 1; ++j) { p1[j]     -= p[j] * nodes[m]; }
        for (size_t j = 0; j < num_nodes + 1; ++j) { p[j]       = p1[j]; }
      }

      return p;
    }


    /**
     * Compute quadrature matrix \\( Q \\) between two sets of nodes.
     *
     * Computing the quadrature matrix \\( Q \\) for polynomial-based integration from one set of
     * quadrature nodes (@p from) to another set of quadrature nodes (@p to).
     *
     * @tparam precision precision of quadrature (i.e. `double`)
     * @param[in] from first set of quadrature nodes
     * @param[in] to second set of quadrature nodes
     * @returns quadrature matrix \\( Q \\) with `to.size()` rows and `from.size()` colums
     *
     * @pre For correctness of the algorithm it is assumed, that both sets of nodes are in the range
     *   \\( [0, 1] \\).
     *
     * @since v0.3.0
     */
    template<typename precision>
    static Matrix<precision> compute_q_matrix(const vector<precision>& from, const vector<precision>& to)
    {
      const size_t to_size = to.size();
      const size_t from_size = from.size();
      assert(to_size >= 1 && from_size >= 1);

      Matrix<precision> q_mat = Matrix<precision>::Zero(to_size, from_size);

      for (size_t m = 0; m < from_size; ++m) {
        Polynomial<precision> p = build_polynomial(m, from);
        // evaluate integrals
        auto den = p.evaluate(from[m]);
        auto P = p.integrate();
        for (size_t j = 0; j < to_size; ++j) {
          q_mat(j, m) = (P.evaluate(to[j]) - P.evaluate(0.0)) / den;
        }
      }

      return q_mat;
    }


    /**
     * Compute quadrature matrix \\( Q \\) for one set of nodes.
     *
     * @tparam precision precision of quadrature (i.e. `double`)
     * @param[in] nodes quadrature nodes to compute \\( Q \\) matrix for
     *
     * @since v0.3.0
     *
     * @overload
     */
    template<typename precision>
    static Matrix<precision> compute_q_matrix(const vector<precision>& nodes)
    {
      return compute_q_matrix<precision>(nodes, nodes);
    }


    /**
     * Compute quadrature matrix \\( Q \\) from a given node-to-node quadrature matrix \\( S \\).
     *
     * @tparam precision precision of quadrature (i.e. `double`)
     * @param[in] s_mat \\( S \\) matrix to compute \\( Q \\) from
     * @see pfasst::quadrature::compute_s_matrix
     *
     * @since v0.3.0
     *
     * @overload
     */
    template<typename precision>
    static Matrix<precision> compute_q_matrix(const Matrix<precision>& s_mat)
    {
      Matrix<precision> q_mat = Matrix<precision>::Zero(s_mat.rows(), s_mat.cols());
      q_mat.col(0) = s_mat.col(0);
      for (Index<precision> q_mat_col = 1; q_mat_col < q_mat.cols(); ++q_mat_col) {
        q_mat.col(q_mat_col) = q_mat.col(q_mat_col - 1) + s_mat.col(q_mat_col);
      }
      return q_mat;
    }


    /**
     * Compute node-to-node quadrature matrix \\( S \\) from a given quadrature matrix \\( Q \\).
     *
     * The \\( S \\) matrix provides a node-to-node quadrature where the \\( i \\)-th row of
     * \\( S \\) represents a quadrature from the \\( i-1 \\)-th node to the \\( i \\)-th node.
     *
     * The procedure is simply subtracting the \\( i-1 \\)-th row of \\( Q \\) from the
     * \\( i \\)-th row of \\( Q \\).
     *
     * @tparam precision precision of quadrature (i.e. `double`)
     * @param[in] q_mat \\( Q \\) matrix to compute \\( S \\) of
     * @returns \\( S \\) matrix
     *
     * @since v0.3.0
     */
    template<typename precision>
    static Matrix<precision> compute_s_matrix(const Matrix<precision>& q_mat)
    {
      Matrix<precision> s_mat = Matrix<precision>::Zero(q_mat.rows(), q_mat.cols());
      s_mat.row(0) = q_mat.row(0);
      for (Index<precision> row = 1; row < s_mat.rows(); ++row) {
        s_mat.row(row) = q_mat.row(row) - q_mat.row(row - 1);
      }
      return s_mat;
    }


    /**
     * Compute node-to-node quadrature matrix \\( S \\) from two given sets of nodes
     *
     * @tparam precision precision of quadrature (i.e. `double`)
     * @param[in] from first set of quadrature nodes
     * @param[in] to second set of quadrature nodes
     *
     * @since v0.3.0
     *
     * @overload
     */
    template<typename precision>
    static Matrix<precision> compute_s_matrix(const vector<precision>& from, const vector<precision>& to)
    {
      return compute_s_matrix(compute_q_matrix(from, to));
    }


    /**
     * Compute vector \\( q \\) for integration from \\( 0 \\) to \\( 1 \\) for given set of nodes.
     *
     * This equals to the last row of the quadrature matrix \\( Q \\) for the given set of nodes.
     *
     * @tparam precision precision of quadrature (i.e. `double`)
     * @param[in] nodes quadrature nodes to compute \\( Q \\) matrix for
     * @pre For correctness of the algorithm it is assumed, that the nodes are in the range
     *   \\( [0, 1] \\).
     *
     * @since v0.3.0
     */
    template<typename precision>
    static vector<precision> compute_q_vec(const vector<precision>& nodes)
    {
      const size_t num_nodes = nodes.size();
      assert(num_nodes >= 1);

      vector<precision> q_vec = vector<precision>(num_nodes, precision(0.0));

      for (size_t m = 0; m < num_nodes; ++m) {
        Polynomial<precision> p = build_polynomial(m, nodes);
        // evaluate integrals
        auto den = p.evaluate(nodes[m]);
        auto P = p.integrate();
        q_vec[m] = (P.evaluate(precision(1.0)) - P.evaluate(precision(0.0))) / den;
      }

      return q_vec;
    }

    /**
     * Interface for quadrature handlers.
     *
     * Quadrature handlers provide \\( Q \\), \\( S \\) and \\( B \\) matrices respecting the left
     * and right nodes, i.e. whether \\( 0 \\) and \\( 1 \\) are part of the nodes or not.
     *
     * Computation of the quadrature nodes and matrices (i.e. quadrature weights) is done on
     * initialization.
     *
     * @tparam precision precision of quadrature (i.e. `double`)
     *
     * @since v0.3.0
     */
    template<typename precision>
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
        /**
         * @throws invalid_argument if number of nodes is invalid for quadrature type
         */
        explicit IQuadrature(const size_t num_nodes);
        /**
         * @throws invalid_argument if number of nodes is invalid for quadrature type
         */
        IQuadrature();
        IQuadrature(const IQuadrature<precision>& other) = default;
        IQuadrature(IQuadrature<precision>&& other) = default;
        virtual ~IQuadrature() = default;
        //! @}

        //! @{
        IQuadrature<precision>& operator=(const IQuadrature<precision>& other) = default;
        IQuadrature<precision>& operator=(IQuadrature<precision>&& other) = default;
        //! @}

        //! @{
        virtual const Matrix<precision>& get_q_mat() const;
        virtual const Matrix<precision>& get_s_mat() const;
        virtual const Matrix<precision>& get_b_mat() const;
        virtual const vector<precision>& get_q_vec() const;
        virtual const vector<precision>& get_nodes() const;
        virtual size_t get_num_nodes() const;

        /**
         * @throws pfasst::NotImplementedYet if not overwritten by implementation;
         *   required for quadrature of any kind
         */
        virtual bool left_is_node() const;
        /**
         * @throws pfasst::NotImplementedYet if not overwritten by implementation;
         *   required for quadrature of any kind
         */
        virtual bool right_is_node() const;
        //! @}

        /**
         * Compute a rough estimate of the numerical error... XXX
         */
        precision expected_error() const;

      protected:
        //! @{
        /**
         * @throws pfasst::NotImplementedYet if not overwritten by implementation;
         *   required for quadrature of any kind
         */
        virtual void compute_nodes();
        virtual void compute_weights();
        //! @}
    };
  }  // ::pfasst::quadrature
}  // ::pfasst

#include "pfasst/quadrature/quadrature_impl.hpp"

#endif  // _PFASST__QUADRATURE__INTERFACE_HPP_
