/*
 * PFASST controller.
 */

#ifndef _PFASST_PFASST_HPP_
#define _PFASST_PFASST_HPP_

#include <mpi.h>

#include "mlsdc.hpp"

using namespace std;

namespace pfasst {

  template<typename time>
  class PFASST : public MLSDC<time> {
    ICommunicator* comm;

    bool predict, initial;

    using LevelIter = typename pfasst::Controller<time>::LevelIter;

    void perform_sweeps(LevelIter leviter, time t, time dt)
    {
      auto* sweeper = leviter.current();
      for (int s=0; s<this->nsweeps[leviter.level]; s++)
	if (predict) {
	  sweeper->predict(t, dt, initial & predict);
	  predict = false;
	} else {
	  sweeper->sweep(t, dt);
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
    void predictor(time t, time dt)
    {
      // restrict fine initial condition
      for (auto leviter=this->finest()-1; leviter>=this->coarsest(); --leviter) {
	auto* crse = leviter.current();
	auto* fine = leviter.fine();
	auto* trns = leviter.transfer();
	trns->restrict(crse, fine, false, true);
	crse->save(true);
      }

      // perform sweeps on the coarse level based on rank
      predict = true;
      auto crse_leviter = this->coarsest();
      for (int nstep=0; nstep<comm->rank()+1; nstep++) {
	initial = nstep == 0;
	perform_sweeps(crse_leviter, t, dt);
	crse_leviter.current()->advance();
      }

      // return to finest level, sweeping as we go
      for (auto leviter=this->coarsest()+1; leviter<=this->finest(); ++leviter) {
	auto* crse = leviter.coarse();
	auto* fine = leviter.current();
	auto* trns = leviter.transfer();

	trns->interpolate(fine, crse, true, true);
	if (leviter < this->finest())
	  perform_sweeps(leviter, t, dt);
      }
    }

    void broadcast()
    {
      cout << "BROADCAST" << endl;
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
      int nblocks = this->nsteps / comm->size();

      initial = true;

      for (int nblock=0; nblock<nblocks; nblock++) {
	int  nstep = nblock * comm->size() + comm->rank();
	time t     = nstep * this->dt;

	predictor(t, this->dt);

	for (int niter=0; niter<this->niters; niter++) {
	  for (auto leviter=this->coarsest(); leviter<=this->finest(); ++leviter)
	    leviter.current()->post(comm);

	  cycle_v(this->finest(), t, this->dt);
	}

	if (nblock < nblocks-1)
	  broadcast();
      }
    }

    /**
     * Cycle down: sweep on current (fine), restrict to coarse.
     */
    LevelIter cycle_down(LevelIter leviter, double t, double dt)
    {
      auto* fine = leviter.current();
      auto* crse = leviter.coarse();
      auto* trns = leviter.transfer();

      // for (int s=0; s<nsweeps[leviter.level]; s++)
      // 	fine->sweep(t, dt);

      trns->restrict(crse, fine);
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
    LevelIter cycle_up(LevelIter leviter, double t, double dt)
    {
      auto* fine = leviter.current();
      auto* crse = leviter.coarse();
      auto* trns = leviter.transfer();

      trns->interpolate(fine, crse, false);

      // if (leviter < this->finest())
      // 	for (int s=0; s<nsweeps[leviter.level]; s++)
      // 	  fine->sweep(t, dt);

      return leviter + 1;
    }

    /**
     * Cycle bottom: sweep on the current (coarsest) level.
     */
    LevelIter cycle_bottom(LevelIter leviter, double t, double dt)
    {
      auto* crse = leviter.current();
      crse->recv(comm);
      perform_sweeps(leviter, t, dt);
      crse->send(comm);
      return leviter + 1;
    }

    /**
     * Perform an MLSDC V-cycle.
     */
    LevelIter cycle_v(LevelIter leviter, double t, double dt)
    {
      if (leviter.level == 0) {
      	leviter = cycle_bottom(leviter, t, dt);
      } else {
      	leviter = cycle_down(leviter, t, dt);
      	leviter = cycle_v(leviter, t, dt);
      	leviter = cycle_up(leviter, t, dt);
      }
      return leviter;
    }

  };

}

#endif
