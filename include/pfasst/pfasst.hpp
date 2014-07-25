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
  class PFASST : public MLSDC<time> {
    ICommunicator* comm;

    typedef typename pfasst::Controller<time>::LevelIter LevelIter;

      bool predict, initial;

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
      // restrict fine initial condition
      for (auto leviter = this->finest() - 1; leviter >= this->coarsest(); --leviter) {
	auto crse = leviter.current();
	auto fine = leviter.fine();
	auto trns = leviter.transfer();
	trns->restrict(crse, fine, false, true);
	crse->save(true);
      }

      // perform sweeps on the coarse level based on rank
      predict = true;
      auto crse = this->coarsest().current();
      for (int nstep=0; nstep<comm->rank()+1; nstep++) {
        this->set_step(comm->rank());
        // XXX: set iteration?

	initial = nstep == 0;
	perform_sweeps(0);
	crse->advance();
      }

      // return to finest level, sweeping as we go
      for (auto leviter=this->coarsest()+1; leviter<=this->finest(); ++leviter) {
	auto crse = leviter.coarse();
	auto fine = leviter.current();
	auto trns = leviter.transfer();

	trns->interpolate(fine, crse, true, true);
	if (leviter < this->finest())
	  perform_sweeps(leviter.level);
      }
    }

    void broadcast()
    {
      throw NotImplementedYet("broadcast");
      this->get_level(this->nlevels()-1)->advance();
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
      for (int nblock=0; nblock<nblocks; nblock++) {
        this->set_step(nblock * comm->size() + comm->rank());

	predictor();

        for (this->set_iteration(0); this->get_iteration() < this->get_max_iterations(); this->advance_iteration()) {
	  for (auto leviter=this->coarsest()+1; leviter<=this->finest(); ++leviter) {
	    int tag = leviter.level*10000 + this->get_iteration() + 10;
	    leviter.current()->post(comm, tag);
	  }

	  perform_sweeps(this->nlevels()-1);
	  // XXX check convergence
	  auto fine = this->get_level(this->nlevels()-1);
	  auto crse = this->get_level(this->nlevels()-2);
	  auto trns = this->get_transfer(this->nlevels()-1);

	  int tag = (this->nlevels()-1)*10000 + this->get_iteration() + 10;
	  fine->send(comm, tag, false);
	  trns->restrict(crse, fine, true, false);
	  trns->fas(this->get_time_step(), crse, fine);
	  crse->save();

	  cycle_v(this->finest()-1);

	  trns->interpolate(fine, crse, false, true, false);
	  fine->recv(comm, tag, false);
          trns->interpolate(fine, crse, false, false, true);
          // XXX: call interpolate_q0(pf,F, G)
	}

	if (nblock < nblocks-1)
	  broadcast();
      }
    }

    /**
     * Cycle down: sweep on current (fine), restrict to coarse.
     */
    LevelIter cycle_down(LevelIter leviter)
    {
      auto fine = leviter.current();
      auto crse = leviter.coarse();
      auto trns = leviter.transfer();

      perform_sweeps(leviter.level);

      int tag = leviter.level*10000 + this->get_iteration() + 10;
      fine->send(comm, tag, false);

      auto dt = this->get_time_step();
      trns->restrict(crse, fine, true, false);
      trns->fas(dt, crse, fine);
      crse->save();

      return leviter - 1;
    }

    /**
     * Cycle up: interpolate coarse correction to fine, sweep on
     * current (fine).
     *
     * Note that if the fine level corresponds to the finest MLSDC
     * level, we don't perform a sweep.  In this case the only
     * operation that is performed here is interpolation.
     */
    LevelIter cycle_up(LevelIter leviter)
    {
      auto fine = leviter.current();
      auto crse = leviter.coarse();
      auto trns = leviter.transfer();

      trns->interpolate(fine, crse, false, true, false);

      int tag = leviter.level*10000 + this->get_iteration() + 10;
      fine->recv(comm, tag, false);
      // XXX          call interpolate_q0(pf,F, G)
      trns->interpolate(fine, crse, false, false, true);

      if (leviter < this->finest()) {
	perform_sweeps(leviter.level);
      }

      return leviter + 1;
    }

    /**
     * Cycle bottom: sweep on the current (coarsest) level.
     */
    LevelIter cycle_bottom(LevelIter leviter)
    {
      auto crse = leviter.current();

      int tag = leviter.level*10000 + this->get_iteration() + 10;
      crse->recv(comm, tag, true);
      this->perform_sweeps(leviter.level);
      crse->send(comm, tag, true);
      return leviter + 1;
    }

    /**
     * Perform an MLSDC V-cycle.
     */
    LevelIter cycle_v(LevelIter leviter)
    {
      // this v-cycle is a bit different than in mlsdc

      if (leviter.level == 0) {
      	leviter = cycle_bottom(leviter);
      } else {
      	leviter = cycle_down(leviter);
      	leviter = cycle_v(leviter);
      	leviter = cycle_up(leviter);
      }
      return leviter;
    }

  };

}

#endif
