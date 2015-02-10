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

  template<typename time>
  void PFASST<time>::run()
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
           this->get_iteration() < this->get_max_iterations() && this->comm->status->keep_iterating();
           this->advance_iteration()) {

        if (this->comm->status->previous_is_iterating()) {
          post();
        }
        cycle_v(this->finest());
      }

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
  void PFASST<time>::set_comm(ICommunicator* comm)
  {
    this->comm = comm;
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

    trns->fas(this->get_time_step(), crse, fine);
    crse->save();

    return l - 1;
  }

  template<typename time>
  typename PFASST<time>::LevelIter PFASST<time>::cycle_up(typename PFASST<time>::LevelIter l)
  {
    auto fine = l.current();
    auto crse = l.coarse();
    auto trns = l.transfer();

    trns->interpolate(fine, crse, true);

    if (this->comm->status->previous_is_iterating()) {
      fine->recv(comm, tag(l), false);
      trns->interpolate_initial(fine, crse);
    }

    if (l < this->finest()) {
      perform_sweeps(l.level);
    }

    return l + 1;
  }

  template<typename time>
  typename PFASST<time>::LevelIter PFASST<time>::cycle_bottom(typename PFASST<time>::LevelIter l)
  {
    auto crse = l.current();

    if (this->comm->status->previous_is_iterating()) {
      crse->recv(comm, tag(l), true);
    }
    this->comm->status->recv();
    this->perform_sweeps(l.level);
    crse->send(comm, tag(l), true);
    this->comm->status->send();
    return l + 1;
  }

  template<typename time>
  typename PFASST<time>::LevelIter PFASST<time>::cycle_v(typename PFASST<time>::LevelIter l)
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
    this->get_finest()->broadcast(comm);
  }

  template<typename time>
  int PFASST<time>::tag(LevelIter l)
  {
    return l.level * 10000 + this->get_iteration() + 10;
  }

  template<typename time>
  void PFASST<time>::post()
  {
    this->comm->status->post();
    for (auto l = this->coarsest() + 1; l <= this->finest(); ++l) {
      l.current()->post(comm, tag(l));
    }
  }
}  // ::pfasst
