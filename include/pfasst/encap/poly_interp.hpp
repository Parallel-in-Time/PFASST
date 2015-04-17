#ifndef _PFASST_ENCAP_POLYINTERP_HPP_
#define _PFASST_ENCAP_POLYINTERP_HPP_

#include <memory>
#include <vector>
using namespace std;

#include "pfasst/interfaces.hpp"
#include "pfasst/encap/encapsulation.hpp"


namespace pfasst
{
  namespace encap
  {
    /**
     * Polynomial time interpolation mixin.
     */
    template<typename time = time_precision>
    class PolyInterpMixin
      : public pfasst::ITransfer<time>
    {
      protected:
        //! @{
        typedef vector<shared_ptr<Encapsulation<time>>> EncapVecT;
        Matrix<time> tmat;
        Matrix<time> fmat;
        //! @}

      public:
        //! @{
        virtual ~PolyInterpMixin();
        //! @}

        //! @{
        virtual void interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                         shared_ptr<const ISweeper<time>> src) override;
        virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                                 shared_ptr<const ISweeper<time>> src,
                                 bool interp_initial) override;
        virtual void interpolate(shared_ptr<Encapsulation<time>> dst,
                                 shared_ptr<const Encapsulation<time>> src);
        //! @}

        //! @{
        virtual void restrict_initial(shared_ptr<ISweeper<time>> dst,
                                      shared_ptr<const ISweeper<time>> src) override;
        virtual void restrict(shared_ptr<ISweeper<time>> dst,
                              shared_ptr<const ISweeper<time>> src,
                              bool restrict_initial) override;
        virtual void restrict(shared_ptr<Encapsulation<time>> dst,
                              shared_ptr<const Encapsulation<time>> src);
        //! @}

        //! @{
        virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                         shared_ptr<const ISweeper<time>> src) override;
        //! @}
    };
  }  // ::pfasst::encap
}  // ::pfasst

#include "pfasst/encap/poly_interp_impl.hpp"

#endif
