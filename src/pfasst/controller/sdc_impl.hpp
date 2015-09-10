#include "pfasst/controller/sdc.hpp"


namespace pfasst
{
  /**
   * For each time step at most a total of Controller::get_max_iterations() iterations will be done
   * on the configured sweeper (which happens to be the only in Controller::levels).
   *
   * On the first iteration on a time step, ISweeper::predict() will get called and
   * ISweeper::sweep() on all subsequent iterations on the same time step.
   * After each ISweeper::post_predict() or ISweeper::post_sweep() call ISweeper::converged() will
   * get checked to shortcut the iterations.
   * The next iteration is introduced by Controller::advance_iteration().
   *
   * Once the maximum number of iterations is reached (or the sweeper has converged before) that,
   * the ISweeper::post_step() hook is triggered followed by ISweeper::advance() to complete the
   * time step and advance to the next via Controller::advance_time() until Controller::get_time()
   * is equal or greater Controller::get_end_time().
   */
  template<typename time>
  void SDC<time>::run()
  {
    auto sweeper = this->get_level(0);

    for (; this->get_time() < this->get_end_time(); this->advance_time()) {
      bool initial = this->get_step() == 0;
      for (this->set_iteration(0);
           this->get_iteration() < this->get_max_iterations();
           this->advance_iteration()) {
        bool predict = this->get_iteration() == 0;
        if (predict) {
          sweeper->predict(initial);
          sweeper->post_predict();
        } else {
          sweeper->sweep();
          sweeper->post_sweep();
        }
        if (sweeper->converged()) {
          break;
        }
      }
      sweeper->post_step();
      if (this->get_time() + this->get_step_size() < this->get_end_time()) {
        sweeper->advance();
      }
    }
  }
}  // ::pfasst
