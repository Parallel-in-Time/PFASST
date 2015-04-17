/**
 * @file examples/boris/injective_transfer.hpp
 * @ingroup BorisFiles
 */
#ifndef _EXAMPLES__BORIS__INJECTIVE_TRANSFER__HPP_
#define _EXAMPLES__BORIS__INJECTIVE_TRANSFER__HPP_

#include <pfasst/logging.hpp>
#include <pfasst/interfaces.hpp>
#include <pfasst/encap/encapsulation.hpp>

#include "boris_sweeper.hpp"


namespace pfasst
{
  namespace examples
  {
    namespace boris
    {
      /**
       * @ingroup Boris
       */
      template<
        typename scalar,
        typename time = pfasst::time_precision
      >
      class InjectiveTransfer
        : public ITransfer<time>
      {
        public:
          typedef typename BorisSweeper<scalar, time>::encap_type encap_type;
          typedef typename BorisSweeper<scalar, time>::acceleration_type force_type;

        public:
          virtual ~InjectiveTransfer()
          {}


          virtual void interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                           shared_ptr<const ISweeper<time>> src) override
          {
            CVLOG(2, "BorisTransfer") << "interpolating initial particle only";
            auto fine = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(fine);
            auto crse = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(crse);
            CVLOG(5, "BorisTransfer") << "coarse:       " << crse->get_start_state();

            auto crse_factory = crse->get_factory();
            auto fine_factory = fine->get_factory();

            shared_ptr<encap_type> crse_delta = dynamic_pointer_cast<encap_type>(crse_factory->create(solution));
            this->restrict(crse_delta, dynamic_pointer_cast<const encap_type>(fine->get_start_state()));
            crse_delta->positions() -= dynamic_pointer_cast<const encap_type>(crse->get_start_state())->positions();
            crse_delta->velocities() -= dynamic_pointer_cast<const encap_type>(crse->get_start_state())->velocities();

            shared_ptr<encap_type> fine_delta = dynamic_pointer_cast<encap_type>(fine_factory->create(solution));
            this->interpolate(fine_delta, crse_delta);
            shared_ptr<encap_type> fine_start = dynamic_pointer_cast<encap_type>(fine->get_start_state());
            fine_start->positions() -= fine_delta->positions();
            fine_start->velocities() -= fine_delta->velocities();

            fine->set_start_state(fine_start);
            CVLOG(5, "BorisTransfer") << "interpolated: " << fine->get_start_state();
          }

          virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                                   shared_ptr<const ISweeper<time>> src,
                                   bool interp_initial = false) override
          {
            CVLOG(2, "BorisTransfer") << "interpolating";
            if (interp_initial) {
              this->interpolate_initial(dst, src);
            }

            auto fine = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(fine);
            auto crse = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(crse);

            auto const ncrse = crse->get_nodes().size(); assert(ncrse >= 1);
            auto const nfine = fine->get_nodes().size(); assert(nfine >= 1);

            auto crse_factory = crse->get_factory();
            auto fine_factory = fine->get_factory();

            vector<shared_ptr<encap_type>> fine_delta(nfine), crse_delta(ncrse);

            for (size_t m = 0; m < fine->get_nodes().size(); ++m) {
              shared_ptr<encap_type> crse_state = dynamic_pointer_cast<encap_type>(crse->get_state(m));
              shared_ptr<encap_type> fine_state = dynamic_pointer_cast<encap_type>(fine->get_state(m));
              CVLOG(5, "BorisTransfer") << "coarse[" << m << "]:       " << crse_state;
              CVLOG(5, "BorisTransfer") << "fine[" << m << "]:         " << fine_state;

              crse_delta[m] = dynamic_pointer_cast<encap_type>(crse_factory->create(solution));
              this->restrict(crse_delta[m], fine_state);
              crse_delta[m]->positions() -= crse_state->positions();
              crse_delta[m]->velocities() -= crse_state->velocities();
              CVLOG(5, "BorisTransfer") << "  crse_delta[" << m << "]: " << crse_delta[m];

              fine_delta[m] = dynamic_pointer_cast<encap_type>(fine_factory->create(solution));
              this->interpolate(fine_delta[m], crse_delta[m]);
              CVLOG(5, "BorisTransfer") << "  fine_delta[" << m << "]: " << fine_delta[m];
              fine_state->positions() -= fine_delta[m]->positions();
              fine_state->velocities() -= fine_delta[m]->velocities();

              fine->set_state(fine_state, m);
              fine->evaluate(m);
              CVLOG(5, "BorisTransfer") << "interpolated[" << m << "]: "
                                << dynamic_pointer_cast<encap_type>(fine->get_state(m));
            }
            fine->save();
          }

          virtual void interpolate(shared_ptr<ParticleCloud<scalar>> dst,
                                   shared_ptr<const ParticleCloud<scalar>> src)
          {
            CVLOG(5, "BorisTransfer") << "interpolate cloud: " << src;
            dst->copy(src);
            CVLOG(5, "BorisTransfer") << "               --> " << dst;
          }

          virtual void interpolate(shared_ptr<ParticleCloudComponent<scalar>> dst,
                                   shared_ptr<const ParticleCloudComponent<scalar>> src)
          {
            CVLOG(5, "BorisTransfer") << "interpolate cmpnt: <" << src << ">" << *(src.get());
            *(dst.get()) = *(src.get());
            CVLOG(5, "BorisTransfer") << "               --> <" << dst << ">" << *(dst.get());
          }


          virtual void restrict_initial(shared_ptr<ISweeper<time>> dst,
                                        shared_ptr<const ISweeper<time>> src) override
          {
            CVLOG(2, "BorisTransfer") << "restricting initial particle only";
            auto coarse = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(coarse);
            auto fine = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(fine);
            CVLOG(5, "BorisTransfer") << "fine:       " << fine->get_start_state();
            coarse->set_start_state(dynamic_pointer_cast<typename BorisSweeper<scalar, time>::encap_type>(fine->get_start_state()));
            CVLOG(5, "BorisTransfer") << "restricted: " << coarse->get_start_state();
          }

          virtual void restrict(shared_ptr<ISweeper<time>> dst,
                                shared_ptr<const ISweeper<time>> src,
                                bool restrict_initial = false) override
          {
            CVLOG(2, "BorisTransfer") << "restricting";
            if (restrict_initial) {
              this->restrict_initial(dst, src);
            }

            auto coarse = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(coarse);
            auto fine = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(fine);

            for (size_t m = 0; m < coarse->get_nodes().size(); ++m) {
              CVLOG(5, "BorisTransfer") << "fine[" << m << "]:       "
                      << dynamic_pointer_cast<ParticleCloud<scalar>>(fine->get_state(m));
              coarse->set_state(dynamic_pointer_cast<const encap_type>(fine->get_state(m)), m);
              CVLOG(5, "BorisTransfer") << "restricted[" << m << "]: "
                      << dynamic_pointer_cast<ParticleCloud<scalar>>(coarse->get_state(m));
              coarse->evaluate(m);
            }
          }

          virtual void restrict(shared_ptr<ParticleCloud<scalar>> dst,
                                shared_ptr<const ParticleCloud<scalar>> src)
          {
            CVLOG(5, "BorisTransfer") << "restrict cloud: " << src;
            dst->copy(src);
            CVLOG(5, "BorisTransfer") << "            --> " << dst;
          }

          virtual void restrict(shared_ptr<ParticleCloudComponent<scalar>> dst,
                                shared_ptr<const ParticleCloudComponent<scalar>> src)
          {
            CVLOG(5, "BorisTransfer") << "restrict cmpnt: <" << src << ">" << *(src.get());
            *(dst.get()) = *(src.get());
            CVLOG(5, "BorisTransfer") << "            --> <" << dst << ">" << *(dst.get());
          }


          virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                           shared_ptr<const ISweeper<time>> src) override
          {
            CVLOG(2, "BorisTransfer") << "computing FAS correction";
            auto crse = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(crse);
            auto fine = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(fine);

            typedef typename BorisSweeper<scalar, time>::acceleration_type force_type;

            auto const ncrse = crse->get_nodes().size(); assert(ncrse >= 1);
            auto const nfine = fine->get_nodes().size(); assert(nfine >= 1);

            auto crse_factory = dynamic_pointer_cast<ParticleCloudFactory<scalar>>(crse->get_factory());
            const size_t crse_nparticle = crse_factory->num_particles();
            const size_t crse_dim = crse_factory->dim();
  //           cloud_component_factory<scalar>(crse_nparticle, crse_dim);
            auto fine_factory = dynamic_pointer_cast<ParticleCloudFactory<scalar>>(fine->get_factory());
            const size_t fine_nparticle = fine_factory->num_particles();
            const size_t fine_dim = fine_factory->dim();
  //           cloud_component_factory<scalar>(fine_nparticle, fine_dim);

            vector<shared_ptr<force_type>> crse_int_q(ncrse), fine_int_q(nfine), rstr_int_q(ncrse);
            vector<shared_ptr<force_type>> crse_int_qq(ncrse), fine_int_qq(nfine), rstr_int_qq(ncrse);
            for (size_t m = 0; m < ncrse; m++) {
              crse_int_q[m] = make_shared<force_type>(cloud_component_factory<scalar>(crse_nparticle, crse_dim));
              crse_int_qq[m] = make_shared<force_type>(cloud_component_factory<scalar>(crse_nparticle, crse_dim));
              rstr_int_q[m] = make_shared<force_type>(cloud_component_factory<scalar>(crse_nparticle, crse_dim));
              rstr_int_qq[m] = make_shared<force_type>(cloud_component_factory<scalar>(crse_nparticle, crse_dim));
            }
            for (size_t m = 0; m < nfine; m++) {
              fine_int_q[m] = make_shared<force_type>(cloud_component_factory<scalar>(fine_nparticle, fine_dim));
              fine_int_qq[m] = make_shared<force_type>(cloud_component_factory<scalar>(fine_nparticle, fine_dim));
            }

            // compute '0 to node' integral on the coarse level
            CVLOG(5, "BorisTransfer") << "computing coarse integral";
            crse->integrate(dt, crse_int_q, crse_int_qq);

            // compute '0 to node' integral on the fine level
            CVLOG(5, "BorisTransfer") << "computing fine integral";
            fine->integrate(dt, fine_int_q, fine_int_qq);

            // restrict '0 to node' fine integral
            CVLOG(5, "BorisTransfer") << "restricting fine integral";
            int trat = (int(nfine) - 1) / (int(ncrse) - 1);
            for (size_t m = 0; m < ncrse; m++) {
              this->restrict(rstr_int_q[m], fine_int_q[m * trat]);
              this->restrict(rstr_int_qq[m], fine_int_qq[m * trat]);
            }

            CVLOG(5, "BorisTransfer") << "get previous node-to-node tau correction";
            // compute 'node to node' tau correction
            vector<shared_ptr<force_type>> tau_q(ncrse), tau_qq(ncrse);
            for (size_t m = 0; m < ncrse; m++) {
              tau_q[m] = crse->get_tau_q_as_force(m);
              tau_qq[m] = crse->get_tau_qq_as_force(m);
              CVLOG(5, "BorisTransfer") << "previous tau_q[" << m << "]:  <" << tau_q[m]  << ">" << *(tau_q[m].get());
              CVLOG(5, "BorisTransfer") << "previous tau_qq[" << m << "]: <" << tau_qq[m] << ">" << *(tau_qq[m].get());
              zero(*(tau_q[m].get()));
              zero(*(tau_qq[m].get()));
            }

            CVLOG(5, "BorisTransfer") << "convert node-to-node tau correction from 0-to-node";
            for (size_t m = 0; m < ncrse; ++m) {
              // compute 0-to-m FAS correction
//               CVLOG(5, "BorisTransfer") << "0-to-m FAS_q[" << m << "]:  "
//                       << *(tau_q[m].get()) << " += " << *(rstr_int_q[m].get()) << " - " << *(fine_int_q[m].get());
              *(tau_q[m].get()) = *(rstr_int_q[m].get()) - *(crse_int_q[m].get());
//               CVLOG(5, "BorisTransfer") << "   --> " << *(tau_q[m].get());
//               CVLOG(5, "BorisTransfer") << "0-to-m FAS_qq[" << m << "]: "
//                       << *(tau_qq[m].get()) << " += " << *(rstr_int_qq[m].get()) << " - " << *(fine_int_qq[m].get());
              *(tau_qq[m].get()) = *(rstr_int_qq[m].get()) - *(crse_int_qq[m].get());
//               CVLOG(5, "BorisTransfer") << "   --> " << *(tau_qq[m].get());

//               CVLOG(5, "BorisTransfer") << "  make it node-to-node";
              // make it a (m-1)-to-m FAS correction
//               for (size_t n = 0; n < m; ++n) {
//                 CVLOG(5, "BorisTransfer") << "  tau_q[" << m << "]  += " << *(rstr_int_q[n].get()) << " - " << *(crse_int_q[n].get());
//                 *(tau_q[m].get()) += *(rstr_int_q[n].get()) - *(crse_int_q[n].get());
//                 CVLOG(5, "BorisTransfer") << "     --> " << *(tau_q[m].get());
//                 CVLOG(5, "BorisTransfer") << "  tau_qq[" << m << "] += " << *(rstr_int_qq[n].get()) << " - " << *(crse_int_qq[n].get());
//                 *(tau_qq[m].get()) += *(rstr_int_qq[n].get()) - *(crse_int_qq[n].get());
//                 CVLOG(5, "BorisTransfer") << "     --> " << *(tau_qq[m].get());
//               }
            }

            for (size_t m = 0; m < ncrse; ++m) {
              CVLOG(5, "BorisTransfer") << "new tau_q[" << m << "]:  <" << tau_q[m]  << ">" << *(tau_q[m].get());
              CVLOG(5, "BorisTransfer") << "new tau_qq[" << m << "]: <" << tau_qq[m] << ">" << *(tau_qq[m].get());
            }
          }
          //! @}
      };
    }
  }
}

#endif  // _EXAMPLES__BORIS__INJECTIVE_TRANSFER__HPP_
