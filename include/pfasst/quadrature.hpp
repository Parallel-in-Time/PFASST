#ifndef _PFASST__QUADRATURE_HPP_
#define _PFASST__QUADRATURE_HPP_

#include <cmath>
#include <exception>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>

#include "config.hpp"
#include "interfaces.hpp"
#include "quadrature/polynomial.hpp"
#include "quadrature/interface.hpp"
#include "quadrature/gauss_lobatto.hpp"
#include "quadrature/gauss_legendre.hpp"
#include "quadrature/gauss_radau.hpp"
#include "quadrature/clenshaw_curtis.hpp"
#include "quadrature/uniform.hpp"

template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using namespace std;

namespace pfasst
{
  namespace quadrature
  {
    template<typename precision = pfasst::time_precision>
    IQuadrature<precision>* quadrature_factory(const size_t nnodes,
                                               const QuadratureType qtype)
    {
      if (qtype == QuadratureType::GaussLegendre) {
        return new GaussLegendre<precision>(nnodes);
      } else if (qtype == QuadratureType::GaussLobatto) {
        return new GaussLobatto<precision>(nnodes);
      } else if (qtype == QuadratureType::GaussRadau) {
        return new GaussRadau<precision>(nnodes);
      } else if (qtype == QuadratureType::ClenshawCurtis) {
        return new ClenshawCurtis<precision>(nnodes);
      } else if (qtype == QuadratureType::Uniform) {
        return new Uniform<precision>(nnodes);
      } else {
        throw ValueError("invalid node type passed to compute_nodes.");
        return nullptr;
      }
    }


    template<typename precision = pfasst::time_precision>
    vector<precision> compute_nodes(size_t nnodes, QuadratureType qtype)
    {
      return quadrature_factory<precision>(nnodes, qtype)->get_nodes();
    }


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

    class Quadrature
    {
      private:
        static void init_config_options(po::options_description& opts)
        {
          opts.add_options()
            ("nodes_type", po::value<string>(), "type of quadrature nodes")
            ("num_nodes", po::value<size_t>(), "number of quadrature nodes");
        }

      public:
        static void enable_config_options(size_t index = -1)
        {
          pfasst::config::Options::get_instance()
            .register_init_function("Quadrature",
                                    function<void(po::options_description&)>(init_config_options),
                                    index);
        }
    };

    namespace config
    {
      // note: GCC fails with "error: explicit template specialization cannot have a storage class"
      //       if this template specialization is also declared 'static'; Clang does not care.
      //      template<>
      pfasst::quadrature::QuadratureType get_value(const string& name)
      {
        const string type = pfasst::config::Options::get_instance().get_variables_map()[name].as<string>();
        if (type == "gauss-lobatto") {
          return pfasst::quadrature::QuadratureType::GaussLobatto;
        } else if (type == "gauss-legendre") {
          return pfasst::quadrature::QuadratureType::GaussLegendre;
        } else if (type == "gauss-radau") {
          return pfasst::quadrature::QuadratureType::GaussRadau;
        } else if (type == "clenshaw-curtis") {
          return pfasst::quadrature::QuadratureType::ClenshawCurtis;
        } else if (type == "uniform") {
          return pfasst::quadrature::QuadratureType::Uniform;
        } else {
          throw invalid_argument("Quadrature type '" + type + "' not known.");
        }
      }

      // note: GCC fails with "error: explicit template specialization cannot have a storage class"
      //       if this template specialization is also declared 'static'; Clang does not care.
      //      template<>
      pfasst::quadrature::QuadratureType get_value(const string& name,
                                                   const pfasst::quadrature::QuadratureType& default_value)
      {
        if (pfasst::config::Options::get_instance().get_variables_map().count(name) == 1) {
          return pfasst::config::get_value<pfasst::quadrature::QuadratureType>(name);
        } else {
          return default_value;
        }
      }

    } // ::pfasst::quadrature::config

  }  // ::pfasst::quadrature

}  // ::pfasst

#endif  // _PFASST__QUADRATURE_HPP_
