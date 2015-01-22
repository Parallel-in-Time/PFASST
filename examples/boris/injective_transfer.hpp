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
      template<
        typename scalar,
        typename time = pfasst::time_precision
      >
      class InjectiveTransfer
        : public ITransfer<time>
      {
        public:
          virtual ~InjectiveTransfer()
          {}


          virtual void interpolate_initial(shared_ptr<ISweeper<time>> dst,
                                           shared_ptr<const ISweeper<time>> src) override
          {
            VLOG_FUNC_START("InjectiveTransfer");
            auto fine = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(fine);
            auto crse = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(crse);
            VLOG(6) << LOG_INDENT << "coarse:       " << crse->get_start_state();

//             typedef typename BorisSweeper<scalar, time>::encap_type encap_type;
//             auto crse_factory = crse->get_factory();
//             auto fine_factory = fine->get_factory();

//             shared_ptr<encap_type> crse_delta = dynamic_pointer_cast<encap_type>(crse_factory->create(solution));
//             this->restrict(crse_delta, dynamic_pointer_cast<const encap_type>(fine->get_start_state()));
//             crse_delta->positions() -= dynamic_pointer_cast<const encap_type>(crse->get_start_state())->positions();
//             crse_delta->velocities() -= dynamic_pointer_cast<const encap_type>(crse->get_start_state())->velocities();

//             shared_ptr<encap_type> fine_delta = dynamic_pointer_cast<encap_type>(fine_factory->create(solution));
//             this->interpolate(fine_delta, crse_delta);
//             VLOG(6) << LOG_INDENT << "coarse correction: " << *(fine_delta.get());

//             dynamic_pointer_cast<encap_type>(fine->get_start_state())->positions() -= fine_delta->positions();
//             dynamic_pointer_cast<encap_type>(fine->get_start_state())->velocities() -= fine_delta->velocities();
            fine->set_start_state(crse->get_start_state());

            VLOG(6) << LOG_INDENT << "interpolated: " << fine->get_start_state();
            VLOG_FUNC_END("InjectiveTransfer");
          }

          virtual void interpolate(shared_ptr<ISweeper<time>> dst,
                                   shared_ptr<const ISweeper<time>> src,
                                   bool interp_initial = false) override
          {
            VLOG_FUNC_START("InjectiveTransfer");
            if (interp_initial) {
              this->interpolate_initial(dst, src);
            }

            auto fine = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(fine);
            auto crse = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(crse);

            typedef typename BorisSweeper<scalar, time>::encap_type encap_type;
            auto crse_factory = crse->get_factory();
            auto fine_factory = fine->get_factory();

            for (size_t m = 0; m < fine->get_nodes().size(); ++m) {
              VLOG(6) << LOG_INDENT << "coarse[" << m << "]:       "
                      << dynamic_pointer_cast<encap_type>(crse->get_state(m));

//               shared_ptr<encap_type> crse_delta = dynamic_pointer_cast<encap_type>(crse_factory->create(solution));
//               LOG(DEBUG) << LOG_INDENT << "crse_delta = restr(fine) - crse";
//               this->restrict(crse_delta, dynamic_pointer_cast<const encap_type>(fine->get_state(m)));
//               LOG(DEBUG) << LOG_INDENT << "  " << crse_delta->positions() << " - " << dynamic_pointer_cast<const encap_type>(crse->get_state(m))->positions();
//               crse_delta->positions() -= dynamic_pointer_cast<const encap_type>(crse->get_state(m))->positions();
//               LOG(DEBUG) << LOG_INDENT << "  " << crse_delta->velocities() << " - " << dynamic_pointer_cast<const encap_type>(crse->get_state(m))->velocities();
//               crse_delta->velocities() -= dynamic_pointer_cast<const encap_type>(crse->get_state(m))->velocities();

//               shared_ptr<encap_type> fine_delta = dynamic_pointer_cast<encap_type>(fine_factory->create(solution));
//               this->interpolate(fine_delta, crse_delta);
//               LOG(DEBUG) << LOG_INDENT << "coarse correction: " << fine_delta;

//               auto fine_state = fine->get_state(m);
//               LOG(DEBUG) << LOG_INDENT << "  " << dynamic_pointer_cast<encap_type>(fine_state)->positions() << " - " << fine_delta->positions();
//               dynamic_pointer_cast<encap_type>(fine_state)->positions() -= fine_delta->positions();
//               LOG(DEBUG) << LOG_INDENT << "  " << dynamic_pointer_cast<encap_type>(fine_state)->velocities() << " - " << fine_delta->velocities();
//               dynamic_pointer_cast<encap_type>(fine_state)->velocities() -= fine_delta->velocities();

              fine->set_state(crse->get_state(m), m);

              fine->evaluate(m);
              VLOG(6) << LOG_INDENT << "interpolated[" << m << "]: "
                      << dynamic_pointer_cast<encap_type>(fine->get_state(m));
            }
            fine->save();
            VLOG_FUNC_END("InjectiveTransfer");
          }

          virtual void interpolate(shared_ptr<ParticleCloud<scalar>> dst,
                                   shared_ptr<const ParticleCloud<scalar>> src)
          {
            VLOG_FUNC_START("InjectiveTransfer");
            *(dst.get()) = *(src.get());
            VLOG_FUNC_END("InjectiveTransfer");
          }

          virtual void interpolate(shared_ptr<ParticleCloudComponent<scalar>> dst,
                                   shared_ptr<const ParticleCloudComponent<scalar>> src)
          {
            VLOG_FUNC_START("InjectiveTransfer");
            *(dst.get()) = *(src.get());
            VLOG_FUNC_END("InjectiveTransfer");
          }


          virtual void restrict_initial(shared_ptr<ISweeper<time>> dst,
                                        shared_ptr<const ISweeper<time>> src) override
          {
            VLOG_FUNC_START("InjectiveTransfer");
            auto coarse = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(coarse);
            auto fine = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(fine);
            VLOG(6) << LOG_INDENT << "fine:       " << fine->get_start_state();
            coarse->set_start_state(fine->get_start_state());
            VLOG(6) << LOG_INDENT << "restricted: " << coarse->get_start_state();
            VLOG_FUNC_END("InjectiveTransfer");
          }

          virtual void restrict(shared_ptr<ISweeper<time>> dst,
                                shared_ptr<const ISweeper<time>> src,
                                bool restrict_initial = false) override
          {
            VLOG_FUNC_START("InjectiveTransfer");
            if (restrict_initial) {
              this->restrict_initial(dst, src);
            }

            auto coarse = dynamic_pointer_cast<BorisSweeper<scalar, time>>(dst);
            assert(coarse);
            auto fine = dynamic_pointer_cast<const BorisSweeper<scalar, time>>(src);
            assert(fine);

            for (size_t m = 0; m < coarse->get_nodes().size(); ++m) {
              VLOG(6) << LOG_INDENT << "fine[" << m << "]:       "
                      << dynamic_pointer_cast<ParticleCloud<scalar>>(fine->get_state(m));
              coarse->set_state(fine->get_state(m), m);
              VLOG(6) << LOG_INDENT << "restricted[" << m << "]: "
                      << dynamic_pointer_cast<ParticleCloud<scalar>>(coarse->get_state(m));
              coarse->evaluate(m);
            }
            coarse->save();
            VLOG_FUNC_END("InjectiveTransfer");
          }

          virtual void restrict(shared_ptr<ParticleCloudComponent<scalar>> dst,
                                shared_ptr<const ParticleCloudComponent<scalar>> src)
          {
            VLOG_FUNC_START("InjectiveTransfer");
            *(dst.get()) = *(src.get());
            VLOG_FUNC_END("InjectiveTransfer");
          }

          virtual void restrict(shared_ptr<ParticleCloud<scalar>> dst,
                                shared_ptr<const ParticleCloud<scalar>> src)
          {
            VLOG_FUNC_START("InjectiveTransfer");
            *(dst.get()) = *(src.get());
            VLOG_FUNC_END("InjectiveTransfer");
          }


          virtual void fas(time dt, shared_ptr<ISweeper<time>> dst,
                           shared_ptr<const ISweeper<time>> src) override
          {
            VLOG_FUNC_START("InjectiveTransfer");
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
            VLOG(6) << LOG_INDENT << "computing coarse integral";
            crse->integrate(dt, crse_int_q, crse_int_qq);

            // compute '0 to node' integral on the fine level
            VLOG(6) << LOG_INDENT << "computing fine integral";
            fine->integrate(dt, fine_int_q, fine_int_qq);

            // restrict '0 to node' fine integral
            int trat = (int(nfine) - 1) / (int(ncrse) - 1);
            for (size_t m = 0; m < ncrse; m++) {
              this->restrict(rstr_int_q[m], fine_int_q[m * trat]);
              this->restrict(rstr_int_qq[m], fine_int_qq[m * trat]);
            }

            // compute 'node to node' tau correction
            vector<shared_ptr<force_type>> tau_q(ncrse), tau_qq(ncrse);
            for (size_t m = 0; m < ncrse; m++) {
              tau_q[m] = crse->get_tau_q_as_force(m);
              tau_qq[m] = crse->get_tau_qq_as_force(m);
              zero(*(tau_q[m].get()));
              zero(*(tau_qq[m].get()));
              VLOG(6) << LOG_INDENT << "previous tau_q[" << m << "]: " << *(tau_q[m].get());
              VLOG(6) << LOG_INDENT << "previous tau_qq[" << m << "]: " << *(tau_qq[m].get());
            }

            for (size_t m = 0; m < ncrse; ++m) {
              // compute 0-to-m FAS correction
              VLOG(7) << LOG_INDENT << "0-to-m FAS_q[" << m << "]: "
                      << *(tau_q[m].get()) << " += " << *(rstr_int_q[m].get()) << " - " << *(fine_int_q[m].get());
              *(tau_q[m].get()) += *(rstr_int_q[m].get()) - *(crse_int_q[m].get());
              VLOG(7) << LOG_INDENT << "0-to-m FAS_qq[" << m << "]: "
                      << *(tau_qq[m].get()) << " += " << *(rstr_int_qq[m].get()) << " - " << *(fine_int_qq[m].get());
              *(tau_qq[m].get()) += *(rstr_int_qq[m].get()) - *(crse_int_qq[m].get());

              // make it a (m-1)-to-m FAS correction
              for (size_t n = 0; n < m; ++n) {
                *(tau_q[m].get()) += *(rstr_int_q[n].get()) - *(crse_int_q[n].get());
                *(tau_qq[m].get()) += *(rstr_int_qq[n].get()) - *(crse_int_qq[n].get());
              }
            }

            for (size_t m = 0; m < ncrse; ++m) {
              VLOG(5) << LOG_INDENT << "new tau_q[" << m << "]: " << *(tau_q[m].get());
              VLOG(5) << LOG_INDENT << "new tau_qq[" << m << "]: " << *(tau_qq[m].get());
            }
            VLOG_FUNC_END("InjectiveTransfer");
          }
          //! @}
      };
    }
  }
}

#endif  // _EXAMPLES__BORIS__INJECTIVE_TRANSFER__HPP_
