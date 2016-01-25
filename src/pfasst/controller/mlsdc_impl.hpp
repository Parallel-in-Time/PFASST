#include "pfasst/controller/mlsdc.hpp"

#include <algorithm>
using namespace std;

#include "pfasst/logging.hpp"


namespace pfasst
{
  /**
   * One sweep is either a call to ISweeper::predict() if MLSDC::predict is `true` followed by
   * triggering the ISweeper::post_predict() hook.
   * In case MLSDC::predict is `false`, ISweeper::sweep() and subsequent ISweeper::post_sweep() are
   * called.
   *
   * @note There is no check on convergence here.
   */
  template<typename time>
  void MLSDC<time>::perform_sweeps(size_t level)
  {
    auto sweeper = this->get_level(level);
    ML_CVLOG(1, "Controller", "on level " << level + 1 << "/" << this->nlevels());
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

  /**
   * @note To overwrite the default behaviour of a single sweep per level, the user must call
   *   MLSDC::set_nsweeps() after MLSDC::setup().
   *
   * @internals
   * On the finest level, ISweeper::setup() is called with `false`, on all other levels with `true`.
   * @endinternals
   */
  template<typename time>
  void MLSDC<time>::setup()
  {
    this->nsweeps.resize(this->nlevels());
    fill(this->nsweeps.begin(), this->nsweeps.end(), 1);
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

      if (this->get_time() + this->get_step_size() < this->get_end_time()) {
        this->get_finest()->advance();
      }
    }
  }

  /**
   * First, MLSDC::perform_sweeps() with current level is called followed by a check on convergence
   * via ISweeper::converged() on the same level.
   *
   * If it has not converged the transfer operator of the current level is used to restrict from
   * current to coarse via ITransfer::restrict() (under consideration of MLSDC::initial).
   * The FAS correction is computed via the same transfer operators ITransfer::fas() method.
   *
   * A call to ISweeper::save() finalizes the restriction on the coarser level.
   */
  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_down(typename MLSDC<time>::LevelIter level_iter)
  {
    auto fine = level_iter.current();
    auto crse = level_iter.coarse();
    auto trns = level_iter.transfer();

    perform_sweeps(level_iter.level);

    if (level_iter == this->finest() && fine->converged()) {
      converged = true;
      return level_iter;
    }

    ML_CVLOG(1, "Controller", "Cycle down onto level " << level_iter.level << "/" << this->nlevels());
    trns->restrict(crse, fine, initial);
    trns->fas(this->get_step_size(), crse, fine);
    crse->save();

    return level_iter - 1;
  }

  /**
   * First the transfer operator of the current level is used (via ITransfer::interpolate()) to
   * interpolate from the coarser level to the current (finer) level.
   *
   * @note If the fine level corresponds to the finest MLSDC level, we don't perform a sweep.
   *   In this case the only operation that is performed here is interpolation.
   */
  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_up(typename MLSDC<time>::LevelIter level_iter)
  {
    auto fine = level_iter.current();
    auto crse = level_iter.coarse();
    auto trns = level_iter.transfer();

    ML_CVLOG(1, "Controller", "Cycle up onto level " << level_iter.level + 1 << "/" << this->nlevels());
    trns->interpolate(fine, crse);

    if (level_iter < this->finest()) {
      perform_sweeps(level_iter.level);
    }

    return level_iter + 1;
  }

  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_bottom(typename MLSDC<time>::LevelIter level_iter)
  {
    perform_sweeps(level_iter.level);
    return level_iter + 1;
  }

  /**
   * @note This method is recursive with two exit points.
   *   The two exit points:
   *     1. if @p level_iter points to the coarsest level
   *     2. if MLSDC::cycle_down() results in a converged state
   */
  template<typename time>
  typename MLSDC<time>::LevelIter MLSDC<time>::cycle_v(typename MLSDC<time>::LevelIter level_iter)
  {
    if (level_iter.level == 0) {
      level_iter = cycle_bottom(level_iter);
    } else {
      level_iter = cycle_down(level_iter);
      if (converged) {
        return level_iter;
      }
      level_iter = cycle_v(level_iter);
      level_iter = cycle_up(level_iter);
    }
    return level_iter;
  }
}  // ::pfasst
