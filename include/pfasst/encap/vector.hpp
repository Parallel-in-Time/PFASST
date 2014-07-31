/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_VECTOR_HPP_
#define _PFASST_VECTOR_HPP_

#include <cstdlib>
#include <algorithm>
#include <memory>
#include <vector>
#include <cassert>

#include <boost/numeric/ublas/matrix.hpp>

#include "encapsulation.hpp"

using namespace std;

namespace pfasst
{
  namespace encap
  {
    /**
     * @tparam scalar
     *     precision and numerical type of the data values
     * @tparam time
     *     precision of the time points; defaults to pfasst::time_precision
     */
    template<typename scalar, typename time = time_precision>
    class VectorEncapsulation
      : public vector<scalar>, public Encapsulation<time>
    {
      public:
        //! @{
        /**
         */
        VectorEncapsulation(int size)
          : vector<scalar>(size)
        {
          zero();
        }

        /**
         * @brief copy constuctor
         * @note delegated to sdt::vector<scalar>
         */
        VectorEncapsulation(const VectorEncapsulation<scalar, time>& other)
          : vector<scalar>(other)
        {}
        /**
         * @throws std::bad_cast
         *     if `other` can not be transformed into pfasst::encap::VectorEncapsulation via
         *     `dynamic_cast`
         */
        VectorEncapsulation(const Encapsulation<time>& other)
          : VectorEncapsulation(dynamic_cast<const VectorEncapsulation<scalar, time>>(other))
        {}

        /**
         * @brief move constructor
         * @note delegated to std::vector<scalar>
         */
        VectorEncapsulation(VectorEncapsulation<scalar, time>&& other)
          : vector<scalar>(other)
        {}
        /**
         * @throws std::bad_cast
         *     if `other` can not be transformed into pfasst::encap::VectorEncapsulation via
         *     `dynamic_cast`
         */
        VectorEncapsulation(Encapsulation<time>&& other)
          : VectorEncapsulation(dynamic_cast<VectorEncapsulation<scalar, time>&&>(other))
        {}
        //! @}

        //! @{
        void zero()
        {
          this->assign(this->size(), scalar(0.0));
        }

        void copy(shared_ptr<const Encapsulation<time>> x)
        {
          shared_ptr<const VectorEncapsulation<scalar, time>> x_cast = dynamic_pointer_cast<const VectorEncapsulation<scalar, time>>(x);
          assert(x_cast);
          this->copy(x_cast);
        }

        void copy(shared_ptr<const VectorEncapsulation<scalar, time>> x)
        {
          std::copy(x->cbegin(), x->cend(), this->begin());
        }
        //! @}

        //! @{
        void saxpy(time a, shared_ptr<const Encapsulation<time>> x)
        {
          shared_ptr<const VectorEncapsulation<scalar, time>> x_cast = dynamic_pointer_cast<const VectorEncapsulation<scalar, time>>(x);
          assert(x_cast);

          this->saxpy(a, x_cast);
        }

        void saxpy(time a, shared_ptr<const VectorEncapsulation<scalar, time>> x)
        {
          assert(this->size() == x->size());
          for (size_t i = 0; i < this->size(); i++)
          { this->data()[i] += a * x->data()[i]; }
        }

        /**
         * @note In case any of the elements of `dst` or `src` can not be transformed via
         *     `dynamic_cast` into pfasst::encap::VectorEncapsulation std::abort is called.
         */
        void mat_apply(vector<shared_ptr<Encapsulation<time>>> dst, time a, matrix<time> mat,
                       vector<shared_ptr<Encapsulation<time>>> src, bool zero = true)
        {
          size_t ndst = dst.size();
          size_t nsrc = src.size();

          vector<shared_ptr<VectorEncapsulation<scalar, time>>> dst_cast(ndst), src_cast(nsrc);
          for (size_t n = 0; n < ndst; n++) {
            dst_cast[n] = dynamic_pointer_cast<VectorEncapsulation<scalar, time>>(dst[n]);
            assert(dst_cast[n]);
          }
          for (size_t m = 0; m < nsrc; m++) {
            src_cast[m] = dynamic_pointer_cast<VectorEncapsulation<scalar, time>>(src[m]);
            assert(src_cast[m]);
          }

          dst_cast[0]->mat_apply(dst_cast, a, mat, src_cast, zero);
        }

        void mat_apply(vector<shared_ptr<VectorEncapsulation<scalar, time>>> dst, time a, matrix<time> mat,
                       vector<shared_ptr<VectorEncapsulation<scalar, time>>> src, bool zero = true)
        {
          size_t ndst = dst.size();
          size_t nsrc = src.size();

          if (zero) { for (auto elem : dst) { elem->zero(); } }

          size_t ndofs = dst[0]->size();
          for (size_t i = 0; i < ndofs; i++) {
            for (size_t n = 0; n < ndst; n++) {
              assert(dst[n]->size() == ndofs);
              for (size_t m = 0; m < nsrc; m++) {
                assert(src[m]->size() == ndofs);
                dst[n]->data()[i] += a * mat(n, m) * src[m]->data()[i];
              }
            }
          }
        }

        /**
         * maximum norm of contained elements.
         *
         * This uses std::max with custom comparison function.
         */
        scalar norm0() const
        {
          return std::max(this->cbegin(), this->cend(),
                          [](scalar a, scalar b) {return std::abs(a) < std::abs(b); } );
        }
        //! @}
    };

    /**
     * @tparam scalar
     *     precision and numerical type of the data values
     * @tparam time
     *     precision of the time points; defaults to pfasst::time_precision
     */
    template<typename scalar, typename time = time_precision>
    class VectorFactory : public EncapFactory<time>
    {
        int size;
      public:
        int dofs() { return size; }
        VectorFactory(const int size) : size(size) { }
        shared_ptr<Encapsulation<time>> create(const EncapType)
        {
          return make_shared<VectorEncapsulation<scalar, time>>(size);
        }
    };

    template<typename scalar, typename time = time_precision>
    VectorEncapsulation<scalar,time>& as_vector(shared_ptr<Encapsulation<time>> x)
    {
      typedef VectorEncapsulation<scalar,time> VectorT;
      shared_ptr<VectorT> y = dynamic_pointer_cast<VectorT>(x);
      assert(y);
      return *y.get();
    }

    template<typename scalar, typename time = time_precision>
    const VectorEncapsulation<scalar,time>& as_vector(shared_ptr<const Encapsulation<time>> x)
    {
      typedef VectorEncapsulation<scalar,time> VectorT;
      shared_ptr<const VectorT> y = dynamic_pointer_cast<const VectorT>(x);
      assert(y);
      return *y.get();
    }

  }
}

#endif
