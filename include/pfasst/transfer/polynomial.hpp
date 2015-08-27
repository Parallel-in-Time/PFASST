#ifndef _PFASST__TRANSFER__POLYNOMIAL_HPP_
#define _PFASST__TRANSFER__POLYNOMIAL_HPP_

#include "pfasst/transfer/transfer.hpp"

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/quadrature.hpp"


namespace pfasst
{
  template<
    class TransferTraits,
    typename Enabled = void
  >
  class PolynomialTransfer
    : public Transfer<TransferTraits, Enabled>
  {
    public:
      typedef          TransferTraits              traits;

      typedef typename traits::coarse_encap_traits coarse_encap_traits;
      typedef typename traits::coarse_encap_type   coarse_encap_type;
      typedef typename traits::coarse_time_type    coarse_time_type;
      typedef typename traits::coarse_spacial_type coarse_spacial_type;

      typedef typename traits::fine_encap_traits   fine_encap_traits;
      typedef typename traits::fine_encap_type     fine_encap_type;
      typedef typename traits::fine_time_type      fine_time_type;
      typedef typename traits::fine_spacial_type   fine_spacial_type;

    protected:
      Matrix<fine_time_type> tmat;
      Matrix<fine_time_type> fmat;

      virtual void setup_tmat(const shared_ptr<quadrature::IQuadrature<typename TransferTraits::fine_time_type>> fine_quad,
                              const shared_ptr<quadrature::IQuadrature<typename TransferTraits::coarse_time_type>> coarse_quad);

    public:
      PolynomialTransfer() = default;
      PolynomialTransfer(const PolynomialTransfer<TransferTraits, Enabled>& other) = default;
      PolynomialTransfer(PolynomialTransfer<TransferTraits, Enabled>&& other) = default;
      virtual ~PolynomialTransfer() = default;
      PolynomialTransfer<TransferTraits, Enabled>& operator=(const PolynomialTransfer<TransferTraits, Enabled>& other) = default;
      PolynomialTransfer<TransferTraits, Enabled>& operator=(PolynomialTransfer<TransferTraits, Enabled>&& other) = default;

      virtual void interpolate_initial(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                                       shared_ptr<typename TransferTraits::fine_sweeper_type> fine);
      virtual void interpolate(const shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                               shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                               const bool initial = false);
      virtual void interpolate_data(const shared_ptr<typename TransferTraits::coarse_encap_type> coarse,
                                    shared_ptr<typename TransferTraits::fine_encap_type> fine);

      virtual void restrict_initial(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                                    shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse);
      virtual void restrict(const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                            shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse,
                            const bool initial = false);
      virtual void restrict_data(const shared_ptr<typename TransferTraits::fine_encap_type> fine,
                                 shared_ptr<typename TransferTraits::coarse_encap_type> coarse);

      virtual void fas(const typename TransferTraits::fine_time_type& dt,
                       const shared_ptr<typename TransferTraits::fine_sweeper_type> fine,
                       shared_ptr<typename TransferTraits::coarse_sweeper_type> coarse);
  };
}  // ::pfasst

#include "pfasst/transfer/polynomial_impl.hpp"

#endif  // _PFASST__TRANSFER__POLYNOMIAL_HPP_
