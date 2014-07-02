/*
 * Polynomial time interpolation mixin.
 */

#ifndef _PFASST_ENCAP_POLYINTERP_HPP_
#define _PFASST_ENCAP_POLYINTERP_HPP_

#include <vector>
#include <cassert>

#include "../interfaces.hpp"

namespace pfasst
{
  namespace encap
  {

    template<typename time = time_precision>
    class PolyInterpMixin : public pfasst::ITransfer<time>
    {
        matrix<time> tmat, fmat;

      public:
        virtual ~PolyInterpMixin() { }

        virtual void interpolate(ISweeper<time>* dst, const ISweeper<time>* src,
                                 bool interp_delta_from_initial,
                                 bool interp_initial)
        {
          auto* fine = dynamic_cast<EncapSweeper<time>*>(dst);
          auto* crse = dynamic_cast<const EncapSweeper<time>*>(src);

          if (tmat.size1() == 0) {
            tmat = pfasst::compute_interp<time>(fine->get_nodes(), crse->get_nodes());
          }

          size_t nfine = fine->get_nodes().size();
          size_t ncrse = crse->get_nodes().size();

          auto* crse_factory = crse->get_factory();
          auto* fine_factory = fine->get_factory();

          vector<Encapsulation<time>*> fine_state(nfine), fine_delta(ncrse);

          for (size_t m = 0; m < nfine; m++) { fine_state[m] = fine->get_state(m); }
          for (size_t m = 0; m < ncrse; m++) { fine_delta[m] = fine_factory->create(solution); }

          if (interp_delta_from_initial) {
            for (size_t m = 1; m < nfine; m++) {
              fine_state[m]->copy(fine_state[0]);
            }
          }

          auto* crse_delta = crse_factory->create(solution);
          size_t m0 = interp_initial ? 0 : 1;
          for (size_t m = m0; m < ncrse; m++) {
            crse_delta->copy(crse->get_state(m));
            if (interp_delta_from_initial) {
              crse_delta->saxpy(-1.0, crse->get_state(0));
            } else {
              crse_delta->saxpy(-1.0, crse->get_saved_state(m));
            }
            interpolate(fine_delta[m], crse_delta);
          }
          delete crse_delta;

          if (!interp_initial) {
            fine_delta[0]->zero();
          }

          fine->get_state(0)->mat_apply(fine_state, 1.0, tmat, fine_delta, false);

          for (size_t m = 0; m < ncrse; m++) { delete fine_delta[m]; }
          for (size_t m = m0; m < nfine; m++) { fine->evaluate(m); }
        }

        virtual void restrict(ISweeper<time>* dst, const ISweeper<time>* src, bool restrict_initial)
        {
          auto* crse = dynamic_cast<EncapSweeper<time>*>(dst);
          auto* fine = dynamic_cast<const EncapSweeper<time>*>(src);

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

        virtual void fas(time dt, ISweeper<time>* dst, const ISweeper<time>* src)
        {
          auto* crse = dynamic_cast<EncapSweeper<time>*>(dst);
          auto* fine = dynamic_cast<const EncapSweeper<time>*>(src);

          size_t ncrse = crse->get_nodes().size();
          assert(ncrse > 1);
          size_t nfine = fine->get_nodes().size();
          assert(nfine >= 1);

          auto* crse_factory = crse->get_factory();
          auto* fine_factory = fine->get_factory();

          vector<Encapsulation<time>*> crse_z2n(ncrse - 1), fine_z2n(nfine - 1), rstr_z2n(ncrse - 1);
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
          vector<Encapsulation<time>*> tau(ncrse - 1);
          for (size_t m = 0; m < ncrse - 1; m++) {
            tau[m] = crse->get_tau(m);
          }

          tau[0]->copy(rstr_z2n[0]);
          tau[0]->saxpy(-1.0, crse_z2n[0]);

          for (size_t m = 1; m < ncrse - 1; m++) {
            tau[m]->copy(rstr_z2n[m]);
            tau[m]->saxpy(-1.0, rstr_z2n[m - 1]);

            tau[m]->saxpy(-1.0, crse_z2n[m]);
            tau[m]->saxpy(1.0, crse_z2n[m - 1]);
          }

          for (size_t m = 0; m < ncrse - 1; m++) { delete crse_z2n[m]; }
          for (size_t m = 0; m < ncrse - 1; m++) { delete rstr_z2n[m]; }
          for (size_t m = 0; m < nfine - 1; m++) { delete fine_z2n[m]; }
        }

        // required for interp/restrict helpers
        virtual void interpolate(Encapsulation<time>* dst, const Encapsulation<time>* src)
        {
          throw NotImplementedYet("mlsdc/pfasst");
        }

        virtual void restrict(Encapsulation<time>* dst, const Encapsulation<time>* src)
        {
          throw NotImplementedYet("mlsdc/pfasst");
        }

    };

  }
}

#endif
