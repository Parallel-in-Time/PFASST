/*
 * Polynomial time interpolation mixin.
 */

#ifndef _PFASST_ENCAP_POLYINTERP_HPP_
#define _PFASST_ENCAP_POLYINTERP_HPP_

#include <vector>
#include <cassert>
#include <memory>

#include "../globals.hpp"
#include "../interfaces.hpp"
#include "../quadrature.hpp"
#include "encap_sweeper.hpp"

namespace pfasst
{
  namespace encap
  {

    template<typename time = time_precision>
    class PolyInterpMixin
      : public pfasst::ITransfer<time>
    {
        //! @{
        typedef vector<shared_ptr<Encapsulation<time>>> EncapVecT;
        Matrix<time> tmat;
        Matrix<time> fmat;
        //! @}

      public:
        //! @{
        virtual ~PolyInterpMixin()
        {}
        //! @}

        //! @{
        virtual void interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                         shared_ptr<const ISweeper<time>> src) override
        {
          auto& fine = as_encap_sweeper(dst);
          auto& crse = as_encap_sweeper(src);

          auto crse_factory = crse.get_factory();
          auto fine_factory = fine.get_factory();

          auto crse_delta = crse_factory->create(solution);
          this->restrict(crse_delta, fine.get_start_state());
          crse_delta->saxpy(-1.0, crse.get_start_state());

          auto fine_delta = fine_factory->create(solution);
          this->interpolate(fine_delta, crse_delta);
          fine.get_start_state()->saxpy(-1.0, fine_delta);

          fine.evaluate(0);
        }


        virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                                 shared_ptr<const ISweeper<time>> src,
                                 bool interp_initial) override
        {
          auto& fine = as_encap_sweeper(dst);
          auto& crse = as_encap_sweeper(src);

          if (tmat.rows() == 0) {
            tmat = pfasst::quadrature::compute_interp<time>(fine.get_nodes(), crse.get_nodes());
          }

          size_t nfine = fine.get_nodes().size();
          size_t ncrse = crse.get_nodes().size();

          auto crse_factory = crse.get_factory();
          auto fine_factory = fine.get_factory();

          EncapVecT fine_state(nfine), fine_delta(ncrse);

          for (size_t m = 0; m < nfine; m++) { fine_state[m] = fine.get_state(m); }
          for (size_t m = 0; m < ncrse; m++) { fine_delta[m] = fine_factory->create(solution); }

          auto crse_delta = crse_factory->create(solution);
          size_t m0 = interp_initial ? 0 : 1;
          for (size_t m = m0; m < ncrse; m++) {
            crse_delta->copy(crse.get_state(m));
            crse_delta->saxpy(-1.0, crse.get_saved_state(m));
            interpolate(fine_delta[m], crse_delta);
          }

          if (!interp_initial) {
            fine_delta[0]->zero();
          }

          fine.get_state(0)->mat_apply(fine_state, 1.0, tmat, fine_delta, false);

          for (size_t m = m0; m < nfine; m++) { fine.evaluate(m); }
        }


        virtual void interpolate(shared_ptr<Encapsulation<time>> dst,
                                 shared_ptr<const Encapsulation<time>> src)
        {
          UNUSED(dst); UNUSED(src);
          throw NotImplementedYet("mlsdc/pfasst");
        }
        //! @}


        //! @{
        virtual void restrict_initial(shared_ptr<ISweeper<time>> dst,
                                      shared_ptr<const ISweeper<time>> src) override
        {
          auto& crse = as_encap_sweeper(dst);
          auto& fine = as_encap_sweeper(src);
          this->restrict(crse.get_start_state(), fine.get_start_state());
        }


        virtual void restrict(shared_ptr<ISweeper<time>> dst,
                              shared_ptr<const ISweeper<time>> src,
                              bool restrict_initial) override
        {
          auto& crse = pfasst::encap::as_encap_sweeper(dst);
          auto& fine = pfasst::encap::as_encap_sweeper(src);

          auto const crse_nodes = crse.get_nodes();
          auto const fine_nodes = fine.get_nodes();
          auto const num_crse = crse_nodes.size();
          auto const num_fine = fine_nodes.size();

          int trat = (int(num_fine) - 1) / (int(num_crse) - 1);

          int m0 = restrict_initial ? 0 : 1;
          for (size_t m = m0; m < num_crse; m++) {
            if (crse_nodes[m] != fine_nodes[m * trat]) {
              throw NotImplementedYet("coarse nodes must be nested");
            }
            this->restrict(crse.get_state(m), fine.get_state(m * trat));
          }

          for (size_t m = m0; m < num_crse; m++) { crse.evaluate(m); }
        }


        virtual void restrict(shared_ptr<Encapsulation<time>> dst,
                              shared_ptr<const Encapsulation<time>> src)
        {
          UNUSED(dst); UNUSED(src);
          throw NotImplementedYet("mlsdc/pfasst");
        }
        //! @}

        virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                         shared_ptr<const ISweeper<time>> src) override
        {
          auto& crse = pfasst::encap::as_encap_sweeper(dst);
          auto& fine = pfasst::encap::as_encap_sweeper(src);

          auto const ncrse = crse.get_nodes().size(); assert(ncrse >= 1);
          auto const nfine = fine.get_nodes().size(); assert(nfine >= 1);

          auto crse_factory = crse.get_factory();
          auto fine_factory = fine.get_factory();

          EncapVecT crse_int(ncrse - 1), fine_int(nfine - 1), rstr_int(ncrse - 1);

          for (size_t m = 0; m < ncrse - 1; m++) { crse_int[m] = crse_factory->create(solution); }
          for (size_t m = 0; m < ncrse - 1; m++) { rstr_int[m] = crse_factory->create(solution); }
          for (size_t m = 0; m < nfine - 1; m++) { fine_int[m] = fine_factory->create(solution); }

          // compute '0 to node' integral on the coarse level
          crse.integrate(dt, crse_int);

          // compute '0 to node' integral on the fine level
          fine.integrate(dt, fine_int);

          // restrict '0 to node' fine integral
          int trat = (int(nfine) - 1) / (int(ncrse) - 1);
          for (size_t m = 1; m < ncrse; m++) {
            this->restrict(rstr_int[m - 1], fine_int[m * trat - 1]);
          }

          // compute 'node to node' tau correction
          EncapVecT tau(ncrse - 1), rstr_and_crse(2 * (ncrse - 1));
          for (size_t m = 0; m < ncrse - 1; m++) { tau[m] = crse.get_tau(m); }
          for (size_t m = 0; m < ncrse - 1; m++) { rstr_and_crse[m] = rstr_int[m]; }
          for (size_t m = 0; m < ncrse - 1; m++) { rstr_and_crse[ncrse - 1 + m] = crse_int[m]; }

          if (fmat.rows() == 0) {
            fmat.resize(ncrse - 1, 2 * (ncrse - 1));
            fmat.fill(0.0);

            for (size_t m = 0; m < ncrse - 1; m++) {
              fmat(m, m) = 1.0;
              fmat(m, ncrse - 1 + m) = -1.0;

              for (size_t n = 0; n < m; n++) {
                fmat(m, n) = -1.0;
                fmat(m, ncrse - 1 + n) = 1.0;
              }
            }
          }

          tau[0]->mat_apply(tau, 1.0, fmat, rstr_and_crse, true);
        }

    };

  }  // ::pfasst::encap
}  // ::pfasst

#endif
