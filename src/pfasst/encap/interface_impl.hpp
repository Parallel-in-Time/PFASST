#include "pfasst/encap/interface.hpp"

#include <cassert>
using namespace std;

#include "pfasst/exceptions.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  namespace encap
  {
    template<class EncapsulationTrait>
    shared_ptr<Encapsulation<EncapsulationTrait>>
    axpy(const typename EncapsulationTrait::time_type& a,
         const shared_ptr<Encapsulation<EncapsulationTrait>> x,
         const shared_ptr<Encapsulation<EncapsulationTrait>> y)
    {
      shared_ptr<Encapsulation<EncapsulationTrait>> result = \
        make_shared<Encapsulation<EncapsulationTrait>>(*y);
      result->scaled_add(a, x);
      return result;
    }

    template<class EncapsulationTrait>
    void
    mat_apply(vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x,
              const typename EncapsulationTrait::time_type& a,
              const Matrix<typename EncapsulationTrait::time_type>& mat,
              const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& y,
              const bool zero_vec_x)
    {
      CLOG_IF(x.size() != (size_t)mat.rows(), WARNING, "ENCAP")
        << "size of result vector (" << x.size()
        << ") does not match result of matrix-vector multiplication (" << mat.rows() << ")";

      if (zero_vec_x) {
        for_each(x.begin(), x.end(),
                 [](shared_ptr<Encapsulation<EncapsulationTrait>> xi) {
                   xi->zero();
                });
      }

      const size_t cols = mat.cols();
      const size_t rows = mat.rows();

      for (size_t n = 0; n < rows; ++n) {
        for (size_t m = 0; m < cols; ++m) {
          x[n]->scaled_add(a * mat(n, m), y[m]);
        }
      }
    }

    template<class EncapsulationTrait>
    vector<shared_ptr<Encapsulation<EncapsulationTrait>>>
    mat_mul_vec(const typename EncapsulationTrait::time_type& a,
                const Matrix<typename EncapsulationTrait::time_type>& mat,
                const vector<shared_ptr<Encapsulation<EncapsulationTrait>>>& x)
    {
      assert((size_t)mat.cols() == x.size());
      const size_t cols = mat.cols();
      const size_t rows = mat.rows();

      // initialize result vector of encaps
      vector<shared_ptr<Encapsulation<EncapsulationTrait>>> result(rows);
      for(auto& ri : result) {
        ri = make_shared<Encapsulation<EncapsulationTrait>>();
        ri->data() = x[0]->get_data();
        ri->zero();
      }

      mat_apply(result, a, mat, x, false);
      return result;
    }

    template<class EncapsulationTrait>
    typename EncapsulationTrait::spacial_type
    norm0(const shared_ptr<Encapsulation<EncapsulationTrait>> x)
    {
      return x->norm0();
    }
  }  // ::pfasst::encap
}  // ::pfasst
