/*
 * Polynomial time interpolation mixin.
 */

#ifndef _PFASST_ENCAP_POLYINTERP_HPP_
#define _PFASST_ENCAP_POLYINTERP_HPP_

#include <vector>
#include <cassert>
#include <memory>

#include "../interfaces.hpp"
#include "encap_sweeper.hpp"

namespace pfasst
{
  namespace encap
  {

    template<typename time = time_precision>
    class PolyInterpMixin
      : public pfasst::ITransfer<time>
    {
        using EncapVecT = vector<shared_ptr<Encapsulation<time>>>;
        matrix<time> tmat, fmat;

      public:
        virtual ~PolyInterpMixin() { }

        virtual void interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                         shared_ptr<const ISweeper<time>> src)
        {
          auto& fine = as_encap_sweeper(dst);
          auto& crse = as_encap_sweeper(src);

          auto crse_factory = crse.get_factory();
          auto fine_factory = fine.get_factory();

          auto crse_delta = crse_factory->create(solution);
          this->restrict(crse_delta, fine.get_state(0));
          crse_delta->saxpy(-1.0, crse.get_state(0));

          auto fine_delta = fine_factory->create(solution);
          this->interpolate(fine_delta, crse_delta);
          fine.get_state(0)->saxpy(-1.0, fine_delta);

          fine.evaluate(0);
        }

        virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                                 shared_ptr<const ISweeper<time>> src,
                                 bool interp_initial)
        {
          auto& fine = as_encap_sweeper(dst);
          auto& crse = as_encap_sweeper(src);

          if (tmat.size1() == 0) {
            tmat = pfasst::compute_interp<time>(fine.get_nodes(), crse.get_nodes());
          }

          size_t nfine = fine.get_nodes().size();
          size_t ncrse = crse.get_nodes().size();

          auto crse_factory = crse.get_factory();
          auto fine_factory = fine.get_factory();

          EncapVecT fine_state(nfine), fine_delta(ncrse);

          for (size_t m = 0; m < nfine; m++) { fine_state[m] = fine.get_state(m); }
          for (size_t m = 0; m < ncrse; m++) { fine_delta[m] = fine_factory->create(solution); }

          // if (interp_delta_from_initial) {
          //   for (size_t m = 1; m < nfine; m++) {
          //     fine_state[m]->copy(fine_state[0]);
          //   }
          // }

          auto crse_delta = crse_factory->create(solution);
          size_t m0 = interp_initial ? 0 : 1;
          for (size_t m = m0; m < ncrse; m++) {
            crse_delta->copy(crse.get_state(m));
            // if (interp_delta_from_initial) {
            //   crse_delta->saxpy(-1.0, crse->get_saved_state(0));
            //   // crse_delta->saxpy(-1.0, crse->get_state(0));
            // } else {
              crse_delta->saxpy(-1.0, crse.get_saved_state(m));
            // }
            interpolate(fine_delta[m], crse_delta);
          }

          if (!interp_initial) {
            fine_delta[0]->zero();
          }

          fine.get_state(0)->mat_apply(fine_state, 1.0, tmat, fine_delta, false);

          for (size_t m = m0; m < nfine; m++) { fine.evaluate(m); }
        }

        virtual void restrict(shared_ptr<ISweeper<time>> dst,
                              shared_ptr<const ISweeper<time>> src,
                              bool restrict_initial,
                              bool restrict_initial_only)
        {
          shared_ptr<EncapSweeper<time>> crse = dynamic_pointer_cast<EncapSweeper<time>>(dst);
          assert(crse);
          shared_ptr<const EncapSweeper<time>> fine = dynamic_pointer_cast<const EncapSweeper<time>>(src);
          assert(fine);

          this->restrict(crse, fine, restrict_initial, restrict_initial_only);
        }

        virtual void restrict(shared_ptr<EncapSweeper<time>> crse,
                              shared_ptr<const EncapSweeper<time>> fine,
                              bool restrict_initial,
                              bool restrict_initial_only)
        {
          if (restrict_initial_only) {
            this->restrict(crse->get_state(0), fine->get_state(0));
            return;
          }

          auto dnodes = crse->get_nodes();
          auto snodes = fine->get_nodes();

          size_t ncrse = crse->get_nodes().size();
          assert(ncrse > 1);
          size_t nfine = fine->get_nodes().size();

          int trat = (int(nfine) - 1) / (int(ncrse) - 1);

          int m0 = restrict_initial ? 0 : 1;
          for (size_t m = m0; m < ncrse; m++) {
            if (dnodes[m] != snodes[m * trat]) {
              throw NotImplementedYet("coarse nodes must be nested");
            }
            this->restrict(crse->get_state(m), fine->get_state(m * trat));
          }

          for (size_t m = m0; m < ncrse; m++) { crse->evaluate(m); }
        }

        virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                         shared_ptr<const ISweeper<time>> src)
        {
          shared_ptr<EncapSweeper<time>> crse = dynamic_pointer_cast<EncapSweeper<time>>(dst);
          assert(crse);
          shared_ptr<const EncapSweeper<time>> fine = dynamic_pointer_cast<const EncapSweeper<time>>(src);
          assert(fine);

          this->fas(dt, crse, fine);
        }

        virtual void fas(time dt, shared_ptr<EncapSweeper<time>> crse,
                         shared_ptr<const EncapSweeper<time>> fine)
        {
          size_t ncrse = crse->get_nodes().size(); assert(ncrse >= 1);
          size_t nfine = fine->get_nodes().size(); assert(nfine >= 1);

          auto crse_factory = crse->get_factory();
          auto fine_factory = fine->get_factory();

          EncapVecT crse_z2n(ncrse - 1), fine_z2n(nfine - 1), rstr_z2n(ncrse - 1);

          for (size_t m = 0; m < ncrse - 1; m++) { crse_z2n[m] = crse_factory->create(solution); }
          for (size_t m = 0; m < ncrse - 1; m++) { rstr_z2n[m] = crse_factory->create(solution); }
          for (size_t m = 0; m < nfine - 1; m++) { fine_z2n[m] = fine_factory->create(solution); }

          // compute '0 to node' integral on the coarse level
          crse->integrate(dt, crse_z2n);
          for (size_t m = 1; m < ncrse - 1; m++) {
            crse_z2n[m]->saxpy(1.0, crse_z2n[m - 1]);
          }

          // compute '0 to node' integral on the fine level
          fine->integrate(dt, fine_z2n);
          for (size_t m = 1; m < nfine - 1; m++) {
            fine_z2n[m]->saxpy(1.0, fine_z2n[m - 1]);
          }

          // restrict '0 to node' fine integral
          int trat = (int(nfine) - 1) / (int(ncrse) - 1);
          for (size_t m = 1; m < ncrse; m++) {
            this->restrict(rstr_z2n[m - 1], fine_z2n[m * trat - 1]);
          }

          // compute 'node to node' tau correction
          EncapVecT tau(ncrse - 1), rstr_and_crse(2 * (ncrse - 1));
          for (size_t m = 0; m < ncrse - 1; m++) { tau[m] = crse->get_tau(m); }
          for (size_t m = 0; m < ncrse - 1; m++) { rstr_and_crse[m] = rstr_z2n[m]; }
          for (size_t m = 0; m < ncrse - 1; m++) { rstr_and_crse[ncrse - 1 + m] = crse_z2n[m]; }

          if (fmat.size1() == 0) {
            fmat.resize(ncrse - 1, 2 * (ncrse - 1));
            fmat.clear();

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

        // required for interp/restrict helpers
        virtual void interpolate(shared_ptr<Encapsulation<time>> /*dst*/,
                                 shared_ptr<const Encapsulation<time>> /*src*/)
        {
          throw NotImplementedYet("mlsdc/pfasst");
        }

        virtual void restrict(shared_ptr<Encapsulation<time>> /*dst*/,
                              shared_ptr<const Encapsulation<time>> /*src*/)
        {
          throw NotImplementedYet("mlsdc/pfasst");
        }

    };

  }
}

#endif
