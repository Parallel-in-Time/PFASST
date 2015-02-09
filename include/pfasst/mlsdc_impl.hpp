#include "mlsdc.hpp"

#include <algorithm>
using namespace std;

#include "logging.hpp"


namespace pfasst
{
  template<typename time>
  void MLSDC<time>::perform_sweeps(size_t level)
  {
    auto sweeper = this->get_level(level);
    CLOG(INFO, "Controller") << "on level " << level + 1 << "/" << this->nlevels();
    for (size_t s = 0; s < this->nsweeps[level]; s++) {
      if (predict) {
        sweeper->predict(initial & predict);
        sweeper->post_predict();
        predict = false;
      } else {
        sweeper->sweep();
        sweeper->post_sweep();
      }
    }
  }

  template<typename time>
  void MLSDC<time>::setup()
  {
    nsweeps.resize(this->nlevels());
    fill(nsweeps.begin(), nsweeps.end(), 1);
    for (auto leviter = this->coarsest(); leviter <= this->finest(); ++leviter) {
      leviter.current()->set_controller(this);
      leviter.current()->setup(leviter != this->finest());
    }
  }

  template<typename time>
  void MLSDC<time>::set_nsweeps(vector<size_t> nsweeps)
  {
    this->nsweeps = nsweeps;
  }

  template<typename time>
  void MLSDC<time>::run()
  {
    for (; this->get_time() < this->get_end_time(); this->advance_time()) {
      predict = true;
      initial = true;
      converged = false;

      for (this->set_iteration(0);
           this->get_iteration() < this->get_max_iterations() && !converged;
           this->advance_iteration()) {
        cycle_v(this->finest());
        initial = false;
      }

      perform_sweeps(this->finest().level);

      for (auto l = this->finest(); l >= this->coarsest(); --l) {
        l.current()->post_step();
      }

      if (this->get_time() + this->get_time_step() < this->get_end_time()) {
        this->get_finest()->advance();
      }
    }
  }

  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_down(typename MLSDC<time>::LevelIter l)
  {
    auto fine = l.current();
    auto crse = l.coarse();
    auto trns = l.transfer();

    perform_sweeps(l.level);

    if (l == this->finest() && fine->converged()) {
      converged = true;
      return l;
    }

    CVLOG(1, "Controller") << "Cycle down onto level " << l.level << "/" << this->nlevels();
    trns->restrict(crse, fine, initial);
    trns->fas(this->get_time_step(), crse, fine);
    crse->save();

    return l - 1;
  }

  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_up(typename MLSDC<time>::LevelIter l)
  {
    auto fine = l.current();
    auto crse = l.coarse();
    auto trns = l.transfer();

    CVLOG(1, "Controller") << "Cycle up onto level " << l.level + 1 << "/" << this->nlevels();
    trns->interpolate(fine, crse);

    if (l < this->finest()) {
      perform_sweeps(l.level);
    }

    return l + 1;
  }

  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_bottom(typename MLSDC<time>::LevelIter l)
  {
    perform_sweeps(l.level);
    return l + 1;
  }

  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_v(typename MLSDC<time>::LevelIter l)
  {
    if (l.level == 0) {
      l = cycle_bottom(l);
    } else {
      l = cycle_down(l);
      if (converged) {
        return l;
      }
      l = cycle_v(l);
      l = cycle_up(l);
    }
    return l;
  }
}  // ::pfasst
