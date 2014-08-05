/*
 * PFASST controller.
 */

#ifndef _PFASST_PFASST_HPP_
#define _PFASST_PFASST_HPP_

#include "mlsdc.hpp"

using namespace std;

namespace pfasst
{

  template<typename time = pfasst::time_precision>
  class PFASST
    : public MLSDC<time>
  {
      ICommunicator* comm;

      typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      bool predict; //<! whether to use a 'predict' sweep
      bool initial; //<! whether we're sweeping from a new initial condition

      void perform_sweeps(size_t level)
      {
        auto sweeper = this->get_level(level);
        for (size_t s = 0; s < this->nsweeps[level]; s++) {
          if (predict) {
            sweeper->predict(initial & predict);
            predict = false;
          } else {
            sweeper->sweep();
          }
        }
      }

    public:
      void set_comm(ICommunicator* comm)
      {
        this->comm = comm;
      }

      /**
       * Predictor: restrict initial down, preform coarse sweeps, return to finest.
       */
      void predictor()
      {
        this->get_finest()->spread();

        // restrict fine initial condition
        for (auto l = this->finest() - 1; l >= this->coarsest(); --l) {
          auto crse = l.current();
          auto fine = l.fine();
          auto trns = l.transfer();
          trns->restrict(crse, fine, true, true);
          crse->spread();
          crse->save();
        }

        // perform sweeps on the coarse level based on rank
        predict = true;
        auto crse = this->coarsest().current();
        for (int nstep = 0; nstep < comm->rank() + 1; nstep++) {
          //          this->set_step(comm->rank());
          // XXX: set iteration?

          perform_sweeps(0);
          if (nstep < comm->rank()) {
            crse->advance();
          }
        }

        // return to finest level, sweeping as we go
        for (auto l = this->coarsest() + 1; l <= this->finest(); ++l) {
          auto crse = l.coarse();
          auto fine = l.current();
          auto trns = l.transfer();

          trns->interpolate(fine, crse, true);
          if (l < this->finest()) {
            perform_sweeps(l.level);
          }
        }
      }

      void broadcast()
      {
        this->get_finest()->broadcast(comm);
        initial = true;
      }

      /**
       * Evolve ODE using PFASST.
       *
       * This assumes that the user has set initial conditions on the
       * finest level.
       *
       * Currently uses "block mode" PFASST with the standard predictor.
       */
      void run()
      {
        int nblocks = int(this->get_end_time() / this->get_time_step()) / comm->size();

        if (nblocks == 0) {
          throw ValueError("invalid duration: there are more time processors than time steps");
        }

        initial = true;
        for (int nblock = 0; nblock < nblocks; nblock++) {
          this->set_step(nblock * comm->size() + comm->rank());

          predictor();

          for (this->set_iteration(0); this->get_iteration() < this->get_max_iterations();
               this->advance_iteration()) {
            for (auto l = this->coarsest() + 1; l <= this->finest(); ++l) {
              int tag = l.level * 10000 + this->get_iteration() + 10;
              l.current()->post(comm, tag);
            }

            perform_sweeps(this->nlevels() - 1);
            // XXX check convergence
            auto fine = this->get_level(this->nlevels() - 1);
            auto crse = this->get_level(this->nlevels() - 2);
            auto trns = this->get_transfer(this->nlevels() - 1);

            int tag = (this->nlevels() - 1) * 10000 + this->get_iteration() + 10;
            fine->send(comm, tag, false);
            trns->restrict(crse, fine, true, false);
            trns->fas(this->get_time_step(), crse, fine);
            crse->save();

            cycle_v(this->finest() - 1);

            trns->interpolate(fine, crse, true);
            fine->recv(comm, tag, false);
            trns->interpolate_initial(fine, crse);
            // XXX: call interpolate_q0(pf,F, G)
          }

          if (nblock < nblocks - 1) {
            broadcast();
          }
        }
      }

      /**
       * Cycle down: sweep on current (fine), restrict to coarse.
       */
      LevelIter cycle_down(LevelIter l)
      {
        auto fine = l.current();
        auto crse = l.coarse();
        auto trns = l.transfer();

        perform_sweeps(l.level);

        int tag = l.level * 10000 + this->get_iteration() + 10;
        fine->send(comm, tag, false);

        auto dt = this->get_time_step();
        trns->restrict(crse, fine, true, false);
        trns->fas(dt, crse, fine);
        crse->save();

        return l - 1;
      }

      /**
       * Cycle up: interpolate coarse correction to fine, sweep on
       * current (fine).
       *
       * Note that if the fine level corresponds to the finest MLSDC
       * level, we don't perform a sweep.  In this case the only
       * operation that is performed here is interpolation.
       */
      LevelIter cycle_up(LevelIter l)
      {
        auto fine = l.current();
        auto crse = l.coarse();
        auto trns = l.transfer();

        trns->interpolate(fine, crse, true);

        int tag = l.level * 10000 + this->get_iteration() + 10;
        fine->recv(comm, tag, false);
        // XXX          call interpolate_q0(pf,F, G)
        trns->interpolate_initial(fine, crse);

        if (l < this->finest()) {
          perform_sweeps(l.level);
        }

        return l + 1;
      }

      /**
       * Cycle bottom: sweep on the current (coarsest) level.
       */
      LevelIter cycle_bottom(LevelIter l)
      {
        auto crse = l.current();

        int tag = l.level * 10000 + this->get_iteration() + 10;
        crse->recv(comm, tag, true);
        this->perform_sweeps(l.level);
        crse->send(comm, tag, true);
        return l + 1;
      }

      /**
       * Perform an MLSDC V-cycle.
       */
      LevelIter cycle_v(LevelIter l)
      {
        // this v-cycle is a bit different than in mlsdc

        if (l.level == 0) {
          l = cycle_bottom(l);
        } else {
          l = cycle_down(l);
          l = cycle_v(l);
          l = cycle_up(l);
        }
        return l;
      }
  };

}  // ::pfasst

#endif
