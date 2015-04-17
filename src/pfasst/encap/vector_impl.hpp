#include "pfasst/encap/vector.hpp"

#include <algorithm>
#include <cassert>
#include <vector>
using namespace std;


namespace pfasst
{
  namespace encap
  {
    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(const size_t size)
      : vector<scalar>(size)
    {
      zero();
    }

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(const VectorEncapsulation<scalar, time>& other)
      : vector<scalar>(other)
    {}

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(const Encapsulation<time>& other)
      : VectorEncapsulation(dynamic_cast<const VectorEncapsulation<scalar, time>>(other))
    {}

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(VectorEncapsulation<scalar, time>&& other)
      : vector<scalar>(other)
    {}

    template<typename scalar, typename time>
    VectorEncapsulation<scalar, time>::VectorEncapsulation(Encapsulation<time>&& other)
      : VectorEncapsulation(dynamic_cast<VectorEncapsulation<scalar, time>&&>(other))
    {}

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::zero()
    {
      this->assign(this->size(), scalar(0.0));
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::copy(shared_ptr<const Encapsulation<time>> x)
    {
      shared_ptr<const VectorEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const VectorEncapsulation<scalar, time>>(x);
      assert(x_cast);
      this->copy(x_cast);
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::copy(shared_ptr<const VectorEncapsulation<scalar, time>> x)
    {
      std::copy(x->cbegin(), x->cend(), this->begin());
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::saxpy(time a, shared_ptr<const Encapsulation<time>> x)
    {
      shared_ptr<const VectorEncapsulation<scalar, time>> x_cast = \
        dynamic_pointer_cast<const VectorEncapsulation<scalar, time>>(x);
      assert(x_cast);

      this->saxpy(a, x_cast);
    }

    template<typename scalar, typename time>
    void VectorEncapsulation<scalar, time>::saxpy(time a, shared_ptr<const VectorEncapsulation<scalar, time>> x)
    {
      assert(this->size() == x->size());
      for (size_t i = 0; i < this->size(); i++)
      { this->data()[i] += a * x->data()[i]; }
    }

    template<typename scalar, typename time>
    void
    VectorEncapsulation<scalar, time>::mat_apply(vector<shared_ptr<Encapsulation<time>>> dst,
                                                 time a, Matrix<time> mat,
                                                 vector<shared_ptr<Encapsulation<time>>> src,
                                                 bool zero)
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

    template<typename scalar, typename time>
    void
    VectorEncapsulation<scalar, time>::mat_apply(vector<shared_ptr<VectorEncapsulation<scalar, time>>> dst,
                                                 time a, Matrix<time> mat,
                                                 vector<shared_ptr<VectorEncapsulation<scalar, time>>> src,
                                                 bool zero)
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

    template<typename scalar, typename time>
    time VectorEncapsulation<scalar, time>::norm0() const
    {
      return std::abs(*std::max_element(this->cbegin(), this->cend(),
                                        [](scalar a, scalar b) {return std::abs(a) < std::abs(b); } ));
    }


    template<typename scalar, typename time>
    VectorFactory<scalar, time>::VectorFactory(const size_t size)
      : size(size)
    {}

    template<typename scalar, typename time>
    size_t VectorFactory<scalar, time>::dofs() const
    {
      return size;
    }

    template<typename scalar, typename time>
    shared_ptr<Encapsulation<time>> VectorFactory<scalar, time>::create(const EncapType)
    {
      return make_shared<VectorEncapsulation<scalar, time>>(this->dofs());
    }


    template<typename scalar, typename time>
    VectorEncapsulation<scalar,time>& as_vector(shared_ptr<Encapsulation<time>> x)
    {
      typedef VectorEncapsulation<scalar,time> VectorT;
      shared_ptr<VectorT> y = dynamic_pointer_cast<VectorT>(x);
      assert(y);
      return *y.get();
    }

    template<typename scalar, typename time>
    const VectorEncapsulation<scalar,time>& as_vector(shared_ptr<const Encapsulation<time>> x)
    {
      typedef VectorEncapsulation<scalar,time> VectorT;
      shared_ptr<const VectorT> y = dynamic_pointer_cast<const VectorT>(x);
      assert(y);
      return *y.get();
    }

  }  // ::pfasst::encap
}  // ::pfasst
