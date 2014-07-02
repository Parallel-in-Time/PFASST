/*
 * Host based encapsulated base sweeper.
 */

#ifndef _PFASST_VECTOR_HPP_
#define _PFASST_VECTOR_HPP_

#include <algorithm>
#include <vector>
#include <cassert>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#define GPCHKERR(err, msg) if ((err) == -1) { perror(msg); return; }

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
    class VectorEncapsulation : public vector<scalar>, public Encapsulation<time>
    {
      public:
        //! @{
        VectorEncapsulation(int size) : vector<scalar>(size)
        {
          zero();
        }
        //! @}

        //! @{
        void zero()
        {
          std::fill(this->begin(), this->end(), 0.0);
        }

        void copy(const Encapsulation<time>* x)
        {
          const VectorEncapsulation<scalar, time>* x_cast = dynamic_cast<const VectorEncapsulation<scalar, time>*>(x);
          assert(x_cast != nullptr);
          this->copy(x_cast);
        }

        void copy(const VectorEncapsulation<scalar, time>* x)
        {
          std::copy(x->begin(), x->end(), this->begin());
        }
        //! @}

        //! @{
        void saxpy(time a, const Encapsulation<time>* x)
        {
          const VectorEncapsulation<scalar, time>* x_cast = dynamic_cast<const VectorEncapsulation<scalar, time>*>(x);
          assert(x_cast != nullptr);

          this->saxpy(a, x_cast);
        }

        void saxpy(time a, const VectorEncapsulation<scalar, time>* x)
        {
          for (size_t i = 0; i < this->size(); i++)
          { this->at(i) += a * x->at(i); }
        }

        void mat_apply(vector<Encapsulation<time>*> dst, time a, matrix<time> mat,
                       vector<Encapsulation<time>*> src, bool zero = true)
        {
          size_t ndst = dst.size();
          size_t nsrc = src.size();

          vector<VectorEncapsulation<scalar, time>*> dst_cast(ndst), src_cast(nsrc);
          for (int n = 0; n < ndst; n++) {
            dst_cast[n] = dynamic_cast<VectorEncapsulation<scalar, time>*>(dst[n]);
            assert(dst_cast[n] != nullptr);
          }
          for (int m = 0; m < nsrc; m++) {
            src_cast[m] = dynamic_cast<VectorEncapsulation<scalar, time>*>(src[m]);
            assert(src_cast[m] != nullptr);
          }

          this->mat_apply(dst_cast, a, mat, src_cast, zero);
        }

        void mat_apply(vector<VectorEncapsulation<scalar, time>*> dst, time a, matrix<time> mat,
                       vector<VectorEncapsulation<scalar, time>*> src, bool zero = true)
        {
          size_t ndst = dst.size();
          size_t nsrc = src.size();

          if (zero) { for (size_t n = 0; n < ndst; n++) { dst[n]->zero(); } }

          size_t ndofs = dst[0]->size();
          for (size_t i = 0; i < ndofs; i++) {
            for (size_t n = 0; n < ndst; n++) {
              for (size_t m = 0; m < nsrc; m++) {
                dst[n]->at(i) += a * mat(n, m) * src[m]->at(i);
              }
            }
          }
        }

        /**
         * maximum norm of contained data.
         */
        scalar norm0() const
        {
          scalar max = 0.0;
          for (size_t i = 0; i < this->size(); i++) {
            scalar v = abs(this->at(i));
            if (v > max) { max = v; }
          }
          return max;
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
        Encapsulation<time>* create(const EncapType)
        {
          return new VectorEncapsulation<scalar, time>(size);
        }
    };

  }
}

#endif
