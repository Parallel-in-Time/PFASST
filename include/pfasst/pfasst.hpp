/*
 * PFASST controller.
 */

#ifndef _PFASST_PFASST_HPP_
#define _PFASST_PFASST_HPP_

#include "mlsdc.hpp"

using namespace std;

namespace pfasst
{

  /**
   * implementation of the PFASST algorithm as described in \cite emmett_pfasst_2012
   */
  template<typename time = pfasst::time_precision>
  class PFASST
    : public MLSDC<time>
  {
      ICommunicator* comm;

      typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      bool predict; //<! whether to use a 'predict' sweep

      void perform_sweeps(size_t level)
      {
        auto sweeper = this->get_level(level);
        for (size_t s = 0; s < this->nsweeps[level]; s++) {
          if (predict) {
            sweeper->predict(predict);
            sweeper->post_predict();
            predict = false;
          } else {
            sweeper->sweep();
            sweeper->post_sweep();
          }
        }
      }

    public:
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
        if (this->comm->size() == 1) {
          pfasst::MLSDC<time>::run();
          return;
        }

        int nblocks = int(this->get_end_time() / this->get_time_step()) / comm->size();

        if (nblocks == 0) {
          throw ValueError("invalid duration: there are more time processors than time steps");
        }

        for (int nblock = 0; nblock < nblocks; nblock++) {
          this->set_step(nblock * comm->size() + comm->rank());

          if (this->comm->size() == 1) {
            predict = true;
          } else {
            predictor();
          }

          for (this->set_iteration(0);
               this->get_iteration() < this->get_max_iterations();
               this->advance_iteration()) {
            post();
            cycle_v(this->finest());
          }

          for (auto l = this->finest(); l >= this->coarsest(); --l) {
            l.current()->post_step();
          }

          if (nblock < nblocks - 1) {
            broadcast();
          }
        }
      }

    private:
      /**
       * Cycle down: sweep on current (fine), restrict to coarse.
       */
      LevelIter cycle_down(LevelIter l)
      {
        auto fine = l.current();
        auto crse = l.coarse();
        auto trns = l.transfer();

        perform_sweeps(l.level);

        if (l == this->finest()) {
          // note: convergence tests belong here
        }

        fine->send(comm, tag(l), false);

        trns->restrict(crse, fine, true);

        trns->fas(this->get_time_step(), crse, fine);
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

        fine->recv(comm, tag(l), false);
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

        crse->recv(comm, tag(l), true);
        this->perform_sweeps(l.level);
        crse->send(comm, tag(l), true);
        return l + 1;
      }

      /**
       * Perform an MLSDC V-cycle.
       */
      LevelIter cycle_v(LevelIter l)
      {
        if (l.level == 0) {
          l = cycle_bottom(l);
        } else {
          l = cycle_down(l);
          l = cycle_v(l);
          l = cycle_up(l);
        }
        return l;
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
          trns->restrict_initial(crse, fine);
          crse->spread();
          crse->save();
        }

        // perform sweeps on the coarse level based on rank
        predict = true;
        auto crse = this->coarsest().current();
        for (int nstep = 0; nstep < comm->rank() + 1; nstep++) {
          // XXX: set iteration and step?
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
      }

      int tag(LevelIter l)
      {
        return l.level * 10000 + this->get_iteration() + 10;
      }

      void post()
      {
        for (auto l = this->coarsest() + 1; l <= this->finest(); ++l) {
          l.current()->post(comm, tag(l));
        }
      }

    public:
      void set_comm(ICommunicator* comm)
      {
        this->comm = comm;
      }

  };

}  // ::pfasst

#endif
