/**
 * @file pfasst/quadrature.hpp
 * @since v0.1.0
 */
#ifndef _PFASST__QUADRATURE_HPP_
#define _PFASST__QUADRATURE_HPP_

#include <cmath>
#include <exception>
#include <vector>
using namespace std;

#include <Eigen/Dense>
template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

#include <boost/math/constants/constants.hpp>

#include "pfasst/config.hpp"
#include "pfasst/interfaces.hpp"
#include "pfasst/quadrature/polynomial.hpp"
#include "pfasst/quadrature/interface.hpp"
#include "pfasst/quadrature/gauss_lobatto.hpp"
#include "pfasst/quadrature/gauss_legendre.hpp"
#include "pfasst/quadrature/gauss_radau.hpp"
#include "pfasst/quadrature/clenshaw_curtis.hpp"
#include "pfasst/quadrature/uniform.hpp"

template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


namespace pfasst
{
  /**
   * Functionality related to computing quadrature nodes and weights.
   *
   * @note Please note, that all quadrature nodes are in the range \\( [0, 1] \\).
   */
  namespace quadrature
  {
    /**
     * Instantiates quadrature handler for given number of nodes and type descriptor.
     *
     * @tparam precision numerical type of the nodes (e.g. `double`)
     * @param[in] nnodes number of quadrature nodes
     * @param[in] qtype type descriptor of the quadrature
     * @returns instance of pfasst::quadrature::IQuadrature of specified type with desired number
     *   of nodes
     * @throws pfasst::ValueError if @p qtype is not a valid quadrature type descriptor
     */
    template<typename precision = pfasst::time_precision>
    shared_ptr<IQuadrature<precision>> quadrature_factory(const size_t nnodes,
                                                          const QuadratureType qtype)
    {
      if (qtype == QuadratureType::GaussLegendre) {
        return make_shared<GaussLegendre<precision>>(nnodes);
      } else if (qtype == QuadratureType::GaussLobatto) {
        return make_shared<GaussLobatto<precision>>(nnodes);
      } else if (qtype == QuadratureType::GaussRadau) {
        return make_shared<GaussRadau<precision>>(nnodes);
      } else if (qtype == QuadratureType::ClenshawCurtis) {
        return make_shared<ClenshawCurtis<precision>>(nnodes);
      } else if (qtype == QuadratureType::Uniform) {
        return make_shared<Uniform<precision>>(nnodes);
      } else {
        throw ValueError("invalid node type passed to compute_nodes.");
        return nullptr;
      }
    }

    /**
     * Compute quadrature nodes for given quadrature type descriptor.
     *
     * @tparam precision numerical type of the nodes (e.g. `double`)
     * @param[in] nnodes number of quadrature nodes to compute
     * @param[in] qtype type descriptor of the quadrature nodes
     * @returns std::vector of quadrature nodes of given type
     *
     * @see pfasst::quadrature::QuadratureType for valid types
     * @see pfasst::quadrature::quadrature_factory for further details
     */
    template<typename precision = pfasst::time_precision>
    vector<precision> compute_nodes(size_t nnodes, QuadratureType qtype)
    {
      return quadrature_factory<precision>(nnodes, qtype)->get_nodes();
    }

    /**
     * @todo write documentation
     *
     * @tparam precision numerical type of the interpolation (e.g. `double`)
     */
    template<typename precision = time_precision>
    Matrix<precision> compute_interp(vector<precision> dst, vector<precision> src)
    {
      const size_t ndst = dst.size();
      const size_t nsrc = src.size();

      Matrix<precision> mat(ndst, nsrc);

      for (size_t i = 0; i < ndst; i++) {
        for (size_t j = 0; j < nsrc; j++) {
          precision den = 1.0;
          precision num = 1.0;

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
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif  // _PFASST__QUADRATURE_HPP_
