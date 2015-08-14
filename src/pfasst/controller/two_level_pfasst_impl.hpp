#include "pfasst/controller/two_level_pfasst.hpp"

#include <cassert>
using namespace std;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT, class CommT>
  TwoLevelPfasst<TransferT, CommT>::TwoLevelPfasst()
    : TwoLevelMLSDC<TransferT, CommT>()
  {
    TwoLevelPfasst<TransferT, CommT>::init_loggers();
    this->set_logger_id("PFASST");
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::init_loggers()
  {
    log::add_custom_logger("PFASST");
    log::add_custom_logger("LVL_COARSE");
    log::add_custom_logger("LVL_FINE");
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::set_options()
  {
    TwoLevelMLSDC<TransferT, CommT>::set_options();
  }


  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::setup()
  {
    assert(this->get_communicator() != nullptr);

    TwoLevelMLSDC<TransferT, CommT>::setup();

    assert(this->get_transfer() != nullptr);

    if (this->get_num_levels() != 2) {
      CLOG(ERROR, "CONTROL") << "Two levels (Sweeper) must have been added for Two-Level-PFASST.";
      throw logic_error("Two-Level-PFASST requires two levels");
    }

    if (this->get_communicator()->get_size() < 2) {
      CLOG(ERROR, "CONTROL") << "Two-Level-PFASST requires at least two processes.";
      throw logic_error("two processes required for Two-Level-PFASST");
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::run()
  {
    Controller<TransferT, CommT>::run();

    assert(this->get_communicator() != nullptr);

    const size_t num_blocks = this->get_num_steps() / this->get_communicator()->get_size();

    if (num_blocks == 0) {
      CLOG(ERROR, "CONTROL") << "Invalid Duration: There are more time processes than time steps.";
      throw logic_error("invalid duration: too many time processes for given time steps");
    }

    for (size_t curr_block = 0; curr_block < num_blocks; ++curr_block) {
      this->status()->step() = curr_block * this->get_communicator()->get_size() + this->get_communicator()->get_rank();
      CLOG(DEBUG, "CONTROL") << "processing time step " << this->get_status()->get_step()
        << " in time block " << curr_block;

      this->predictor();

      // iterate on each time step
      do {
        shared_ptr<Status<typename TransferT::traits::fine_time_type>> prev_rank_status;
        if (!this->get_communicator()->is_first()) {
          prev_rank_status = make_shared<Status<typename TransferT::traits::fine_time_type>>();
        }

        this->status()->state() = State::ITERATING;

        this->sweep_fine();

        this->cycle_down();

        if (!this->get_communicator()->is_first()) {
          if (this->get_status()->get_state() > State::FAILED) {
            assert(this->get_coarse()->get_initial_state() != nullptr);
            this->get_coarse()->initial_state()->recv(this->get_communicator(),
                                                      this->get_communicator()->get_rank() - 1,
                                                      this->compute_tag(), true);
          }

          prev_rank_status->recv(this->get_communicator(),
                                 this->get_communicator()->get_rank() - 1,
                                 this->compute_tag(), true);
        }

        this->sweep_coarse();

        // TODO: set status here !?

        if (!this->get_communicator()->is_last()) {
          assert(this->get_coarse()->get_end_state() != nullptr);
          this->get_coarse()->get_end_state()->send(this->get_communicator(),
                                                    this->get_communicator()->get_rank() + 1,
                                                    this->compute_tag(), true);

          this->get_status()->send(this->get_communicator(),
                                   this->get_communicator()->get_rank() - 1,
                                   this->compute_tag(), true);
        }

        this->cycle_up();

        // convergence check
        const bool fine_converged = this->get_fine()->converged();
        const bool previous_done = (this->get_communicator()->is_first())
                                     ? true : prev_rank_status->get_state() <= State::FAILED;
        if (previous_done && fine_converged) {
          this->status()->state() = State::CONVERGED;
        }
      } while(this->advance_iteration());

      this->broadcast();
    }
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelPfasst<TransferT, CommT>::advance_time(const size_t& num_steps)
  {
    return Controller<TransferT, CommT>::advance_time(num_steps);
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelPfasst<TransferT, CommT>::advance_iteration()
  {
    return Controller<TransferT, CommT>::advance_iteration();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::cycle_down()
  {
    // TODO: check convergence state here !?

    if (!this->get_communicator()->is_last()) {
      assert(this->get_fine()->get_end_state() != nullptr);
      this->get_fine()->get_end_state()->send(this->get_communicator(),
                                              this->get_communicator()->get_rank() + 1,
                                              this->compute_tag(), false);
    }

    this->get_transfer()->restrict(this->get_fine(), this->get_coarse(), true);
    this->get_transfer()->fas(this->get_status()->get_dt(), this->get_fine(), this->get_coarse());
    this->get_coarse()->save();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::cycle_up()
  {
    this->get_transfer()->interpolate(this->get_coarse(), this->get_fine(), true);
    if (!this->get_communicator()->is_first()) {
      assert(this->get_fine()->get_initial_state() != nullptr);
      this->get_fine()->initial_state()->recv(this->get_communicator(),
                                              this->get_communicator()->get_rank() - 1,
                                              this->compute_tag(), false);
    }
    this->get_transfer()->interpolate_initial(this->get_coarse(), this->get_fine());
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::predictor()
  {
    assert(this->get_status()->get_iteration() == 0);
    this->status()->state() = State::PREDICTING;

    this->get_fine()->spread();

    // restrict fine initial condition
    this->get_transfer()->restrict_initial(this->get_fine(), this->get_coarse());
    this->get_coarse()->spread();
    this->get_coarse()->save();

    // perform sweeps on coarse level
    for (size_t predict_step = 0;
         predict_step <= this->get_communicator()->get_rank();
         ++predict_step) {
      if (predict_step == 0) {
        this->predict_coarse();
      } else {
        this->sweep_coarse();
      }

      this->get_coarse()->save();
    }

    // return to fine level
    this->get_transfer()->interpolate(this->get_coarse(), this->get_fine(), true);
    this->sweep_fine();
    this->get_fine()->save();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::broadcast()
  {
    this->get_fine()->get_end_state()->bcast(this->get_communicator(),
                                             this->get_communicator()->get_root());
  }

  template<class TransferT, class CommT>
  int
  TwoLevelPfasst<TransferT, CommT>::compute_tag() const
  {
    // TODO: come up with a good way of computing unique but meaningful tags
    return 0;
  }
}  // ::pfasst
