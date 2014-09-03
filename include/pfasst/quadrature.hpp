#ifndef _PFASST_QUADRATURE_HPP_
#define _PFASST_QUADRATURE_HPP_

#include <cmath>
#include <exception>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>

#include "interfaces.hpp"
#include "polynomial.hpp"


#include <boost/math/constants/constants.hpp>
using namespace boost::math::constants;
template<typename scalar>
using Matrix = Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

template<typename scalar>
using IndexT = typename Matrix<scalar>::Index;


using namespace std;
using namespace boost::math::constants;


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

    template<typename precision = time_precision>
    class IQuadrature;
    template<typename precision = time_precision>
    class GaussLegendre;
    template<typename precision = time_precision>
    class GaussLobatto;
    template<typename precision = time_precision>
    class GaussRadau;
    template<typename precision = time_precision>
    class ClenshawCurtis;
    template<typename precision = time_precision>
    class Uniform;

    typedef integral_constant<QuadratureType, QuadratureType::GaussLegendre> gauss_legendre;
    typedef integral_constant<QuadratureType, QuadratureType::GaussLobatto> gauss_lobatto;
    typedef integral_constant<QuadratureType, QuadratureType::GaussRadau> gauss_radau;
    typedef integral_constant<QuadratureType, QuadratureType::ClenshawCurtis> clenshaw_curtis;
    typedef integral_constant<QuadratureType, QuadratureType::Uniform> uniform;
    typedef integral_constant<QuadratureType, QuadratureType::UNDEFINED> undefined;


    template<typename QuadratureT>
    struct quadrature_traits
    {
      typedef pfasst::quadrature::undefined integral_constant;
      static const bool left_is_node = false;
      static const bool right_is_node = false;
    };

    template<>
    struct quadrature_traits<GaussLegendre<>>
    {
      typedef pfasst::quadrature::gauss_legendre integral_constant;
      static const bool left_is_node = false;
      static const bool right_is_node = false;
    };

    template<>
    struct quadrature_traits<GaussLobatto<>>
    {
      typedef pfasst::quadrature::gauss_lobatto integral_constant;
      static const bool left_is_node = true;
      static const bool right_is_node = true;
    };

    template<>
    struct quadrature_traits<GaussRadau<>>
    {
      typedef pfasst::quadrature::gauss_radau integral_constant;
      static const bool left_is_node = false;
      static const bool right_is_node = true;
    };

    template<>
    struct quadrature_traits<ClenshawCurtis<>>
    {
      typedef pfasst::quadrature::clenshaw_curtis integral_constant;
      static const bool left_is_node = true;
      static const bool right_is_node = true;
    };

    template<>
    struct quadrature_traits<Uniform<>>
    {
      typedef pfasst::quadrature::uniform integral_constant;
      static const bool left_is_node = true;
      static const bool right_is_node = true;
    };


    template<typename scalar>
    static Polynomial<scalar> build_polynomial(const size_t node, const vector<scalar>& from, const vector<scalar>& to)
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
      for (IndexT<scalar> q_mat_col = 1; q_mat_col < q_mat.cols(); ++q_mat_col) {
        q_mat.col(q_mat_col) = q_mat.col(q_mat_col - 1) + s_mat.col(q_mat_col);
      }
      return q_mat;
    }


    template<typename scalar>
    static Matrix<scalar> compute_s_matrix(const Matrix<scalar>& q_mat)
    {
      Matrix<scalar> s_mat = Matrix<scalar>::Zero(q_mat.rows(), q_mat.cols());
      s_mat.row(0) = q_mat.row(0);
      for (IndexT<scalar> row = 1; row < s_mat.rows(); ++row) {
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


  //  enum class QuadratureMatrix { S, Q, QQ }; // returning QQ might be cool for 2nd-order stuff
    enum class QuadratureMatrix { S, Q };


    template<typename precision /*= time_precision*/>
    class IQuadrature
    {
      protected:
        //! @{
        size_t m_num_nodes;
        Matrix<precision> m_q_mat;
        Matrix<precision> m_s_mat;
        vector<precision> m_q_vec;
        vector<precision> m_nodes;
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

        virtual size_t num_nodes() const
        { return this->m_num_nodes; }

        virtual bool left_is_node() const
        { return quadrature_traits<IQuadrature<>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<IQuadrature<>>::right_is_node; }
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
        //! @}
    };


    template<typename precision /*= time_precision*/>
    class GaussLobatto
      : public IQuadrature<precision>
    {
      public:
        //! @{
        GaussLobatto(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->m_num_nodes < 2) {
            throw invalid_argument("Gauss-Lobatto quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
        }

        GaussLobatto()
          : IQuadrature<precision>()
        {}

        GaussLobatto(const GaussLobatto<precision>& other)
          : IQuadrature<precision>(other)
        {}

        GaussLobatto(GaussLobatto<precision>&& other)
          : GaussLobatto<precision>()
        {
          swap(*this, other);
        }

        virtual ~GaussLobatto()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return quadrature_traits<GaussLobatto<>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<GaussLobatto<>>::right_is_node; }
        //! @}

        //! @{
        GaussLobatto<precision>& operator=(GaussLobatto<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override
        {
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          auto roots = Polynomial<precision>::legendre(this->m_num_nodes - 1).differentiate().roots();

          for (size_t j = 0; j < this->m_num_nodes - 2; j++) {
            this->m_nodes[j + 1] = 0.5 * (1.0 + roots[j]);
          }
          this->m_nodes.front() = 0.0;
          this->m_nodes.back() = 1.0;
        }
        //! @}
    };


    template<typename precision /*= time_precision*/>
    class GaussLegendre
      : public IQuadrature<precision>
    {
      public:
        //! @{
        GaussLegendre(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          this->compute_nodes();
          this->compute_weights();
        }

        GaussLegendre()
          : IQuadrature<precision>()
        {}

        GaussLegendre(const GaussLegendre<precision>& other)
          : IQuadrature<precision>(other)
        {}

        GaussLegendre(GaussLegendre<precision>&& other)
          : GaussLegendre<precision>()
        {
          swap(*this, other);
        }

        virtual ~GaussLegendre()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return quadrature_traits<GaussLegendre<>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<GaussLegendre<>>::right_is_node; }
        //! @}

        //! @{
        GaussLegendre<precision>& operator=(GaussLegendre<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override
        {
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          auto roots = Polynomial<precision>::legendre(this->m_num_nodes).roots();

          for (size_t j = 0; j < this->m_num_nodes; j++) {
            this->m_nodes[j] = 0.5 * (1.0 + roots[j]);
          }
        }
        //! @}
    };


    template<typename precision /*= time_precision*/>
    class GaussRadau
      : public IQuadrature<precision>
    {
      public:
        //! @{
        GaussRadau(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->m_num_nodes < 2) {
            throw invalid_argument("Gauss-Radau quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
        }

        GaussRadau()
          : IQuadrature<precision>()
        {}

        GaussRadau(const GaussRadau<precision>& other)
          : IQuadrature<precision>(other)
        {}

        GaussRadau(GaussRadau<precision>&& other)
          : GaussRadau<precision>()
        {
          swap(*this, other);
        }

        virtual ~GaussRadau()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return quadrature_traits<GaussRadau<>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<GaussRadau<>>::right_is_node; }
        //! @}

        //! @{
        GaussRadau<precision>& operator=(GaussRadau<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes() override
        {
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          auto l   = Polynomial<precision>::legendre(this->m_num_nodes);
          auto lm1 = Polynomial<precision>::legendre(this->m_num_nodes - 1);

          for (size_t i = 0; i < this->m_num_nodes; i++) {
            l[i] += lm1[i];
          }
          auto roots = l.roots();
          for (size_t j = 1; j < this->m_num_nodes; j++) {
            this->m_nodes[j - 1] = 0.5 * (1.0 - roots[this->m_num_nodes - j]);
          }
          this->m_nodes.back() = 1.0;
        }
        //! @}
    };


    template<typename precision /*= time_precision*/>
    class ClenshawCurtis
      : public IQuadrature<precision>
    {
      public:
        //! @{
        ClenshawCurtis(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->m_num_nodes < 2) {
            throw invalid_argument("Clenshaw-Curtis quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
        }

        ClenshawCurtis()
          : IQuadrature<precision>()
        {}

        ClenshawCurtis(const ClenshawCurtis<precision>& other)
          : IQuadrature<precision>(other)
        {}

        ClenshawCurtis(ClenshawCurtis<precision>&& other)
          : ClenshawCurtis<precision>()
        {
          swap(*this, other);
        }

        virtual ~ClenshawCurtis()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return quadrature_traits<ClenshawCurtis<>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<ClenshawCurtis<>>::right_is_node; }
        //! @}

        //! @{
        ClenshawCurtis<precision>& operator=(ClenshawCurtis<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes()
        {
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          auto roots = Polynomial<precision>::legendre(this->m_num_nodes).roots();

          for (size_t j = 0; j < this->m_num_nodes; j++) {
            this->m_nodes[j] = 0.5 * (1.0 - cos(j * pi<precision>() / (this->m_num_nodes - 1)));
          }
        }
        //! @}
    };


    template<typename precision /*= time_precision*/>
    class Uniform
      : public IQuadrature<precision>
    {
      public:
        //! @{
        Uniform(const size_t num_nodes)
          : IQuadrature<precision>(num_nodes)
        {
          if (this->m_num_nodes < 2) {
            throw invalid_argument("Uniform quadrature requires at least two quadrature nodes.");
          }
          this->compute_nodes();
          this->compute_weights();
        }

        Uniform()
          : IQuadrature<precision>()
        {}

        Uniform(const Uniform<precision>& other)
          : IQuadrature<precision>(other)
        {}

        Uniform(Uniform<precision>&& other)
          : Uniform<precision>()
        {
          swap(*this, other);
        }

        virtual ~Uniform()
        {}
        //! @}

        //! @{
        virtual bool left_is_node() const
        { return quadrature_traits<Uniform<>>::left_is_node; }

        virtual bool right_is_node() const
        { return quadrature_traits<Uniform<>>::right_is_node; }
        //! @}

        //! @{
        Uniform<precision>& operator=(Uniform<precision> other)
        {
          swap(*this, other);
          return *this;
        }
        //! @}

      protected:
        //! @{
        virtual void compute_nodes()
        {
          this->m_nodes = vector<precision>(this->m_num_nodes, precision(0.0));
          for (size_t j = 0; j < this->m_num_nodes; j++) {
            this->m_nodes[j] = precision(j) / (this->m_num_nodes - 1);
          }
        }
        //! @}
    };


    template<typename precision = time_precision>
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


    template<typename precision = time_precision>
    vector<precision> compute_nodes(size_t nnodes, QuadratureType qtype)
    {
      return quadrature_factory<precision>(nnodes, qtype)->nodes();
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
  }  // ::pfasst::quadrature
}  // ::pfasst

#endif
