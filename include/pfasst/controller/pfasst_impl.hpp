#include "pfasst/controller/pfasst.hpp"

#include "pfasst/logging.hpp"


namespace pfasst
{
  template<typename time>
  void PFASST<time>::perform_sweeps(size_t level)
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

  /**
   * @internals
   * @note No consistency checks on validity of given communicator are done.
   * @endinternals
   */
  template<typename time>
  void PFASST<time>::set_comm(ICommunicator* comm)
  {
    this->comm = comm;
  }

  /**
   * @note Currently uses _block mode_ PFASST with the standard predictor
   *   (see PFASST::predictor()).
   */
  template<typename time>
  void PFASST<time>::run()
  {
    if (this->comm->size() == 1) {
      pfasst::MLSDC<time>::run();
      return;
    }

    int nblocks = int(this->get_end_time() / this->get_step_size()) / comm->size();

    if (nblocks == 0) {
      ML_CLOG(INFO, "Controller", "invalid duration: there are more time processors than time steps");
      throw ValueError("invalid duration: there are more time processors than time steps");
    }

    if (nblocks * comm->size() * this->get_step_size() < this->get_end_time()) {
      ML_CLOG(INFO, "Controller", "invalid duration: mismatch between number of time processors and time steps");
      throw ValueError("invalid duration: mismatch between number of time processors and time steps");
    }

    for (int nblock = 0; nblock < nblocks; nblock++) {
      this->set_step(nblock * comm->size() + comm->rank());

      predictor();

      ML_CLOG(DEBUG, "Controller", "iterating on step " << this->get_step()
                                   << " (0/" << this->get_max_iterations() << ")");
      for (this->set_iteration(0);
           this->get_iteration() < this->get_max_iterations() && this->comm->status->keep_iterating();
           this->advance_iteration()) {

        if (this->comm->status->previous_is_iterating()) {
          post();
        }
        cycle_v(this->finest());
      }
      ML_CLOG(DEBUG, "Controller", "done iterating on step " << this->get_step()
                                   << " (" << this->get_iteration() << "/"
                                   << this->get_max_iterations() << ")");

      for (auto l = this->finest(); l >= this->coarsest(); --l) {
        l.current()->post_step();
      }

      if (nblock < nblocks - 1) {
        broadcast();
      }

      this->comm->status->clear();
    }
  }

  template<typename time>
  typename PFASST<time>::LevelIter PFASST<time>::cycle_down(typename PFASST<time>::LevelIter l)
  {
    auto fine = l.current();
    auto crse = l.coarse();
    auto trns = l.transfer();

    perform_sweeps(l.level);

    if (l == this->finest() && fine->converged()) {
      this->comm->status->set_converged(true);
    }

    fine->send(comm, tag(l), false);
    trns->restrict(crse, fine, true);
    trns->fas(this->get_step_size(), crse, fine);
    crse->save();

    return l - 1;
  }

  template<typename time>
  typename PFASST<time>::LevelIter PFASST<time>::cycle_up(typename PFASST<time>::LevelIter level_iter)
  {
    auto fine = level_iter.current();
    auto crse = level_iter.coarse();
    auto trns = level_iter.transfer();

    trns->interpolate(fine, crse, true);
    fine->recv(comm, tag(level_iter), false);
    trns->interpolate_initial(fine, crse);

    if (level_iter < this->finest()) {
      perform_sweeps(level_iter.level);
    }

    return level_iter + 1;
  }

  template<typename time>
  typename PFASST<time>::LevelIter PFASST<time>::cycle_bottom(typename PFASST<time>::LevelIter level_iter)
  {
    auto crse = level_iter.current();

    if (this->comm->status->previous_is_iterating()) {
      crse->recv(comm, tag(level_iter), true);
    }
    this->comm->status->recv(stag(level_iter));
    this->perform_sweeps(level_iter.level);
    crse->send(comm, tag(level_iter), true);
    this->comm->status->set_converged(!this->comm->status->keep_iterating());
    this->comm->status->send(stag(level_iter));
    return level_iter + 1;
  }

  template<typename time>
  typename PFASST<time>::LevelIter PFASST<time>::cycle_v(typename PFASST<time>::LevelIter level_iter)
  {
    if (level_iter.level == 0) {
      level_iter = cycle_bottom(level_iter);
    } else {
      level_iter = cycle_down(level_iter);
      level_iter = cycle_v(level_iter);
      level_iter = cycle_up(level_iter);
    }
    return level_iter;
  }

  template<typename time>
  void PFASST<time>::predictor()
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

  template<typename time>
  void PFASST<time>::broadcast()
  {
    this->get_finest()->broadcast(this->comm);
  }

  /**
   * @internals
   * A simple formula is used with current level index \\( L \\) (provided by @p level_iter) and
   * current iteration number \\( I \\):
   * \\[ (L+1) * 10000 + I \\]
   * @endinternals
   */
  template<typename time>
  int PFASST<time>::tag(LevelIter level_iter)
  {
    return (level_iter.level+1) * 10000 + this->get_iteration();
  }

  template<typename time>
  int PFASST<time>::stag(LevelIter level_iter)
  {
    return level_iter.level * 1000 + this->get_iteration();
  }

  /**
   * @see IStatus::post() for details and implementations of posting current status
   * @see ISweeper::post() for details and implementations of posting current level
   */
  template<typename time>
  void PFASST<time>::post()
  {
    if (this->comm->status->previous_is_iterating()) {
      this->comm->status->post(0);
      for (auto l = this->coarsest() + 1; l <= this->finest(); ++l) {
        l.current()->post(comm, tag(l));
      }
    }
  }
}  // ::pfasst
