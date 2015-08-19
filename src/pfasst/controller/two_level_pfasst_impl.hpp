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
    this->_prev_status = make_shared<Status<time_type>>();
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
      CLOG(ERROR, this->get_logger_id()) << "Two levels (Sweeper) must have been added for Two-Level-PFASST.";
      throw logic_error("Two-Level-PFASST requires two levels");
    }

    if (this->get_communicator()->get_size() < 2) {
      CLOG(ERROR, this->get_logger_id()) << "Two-Level-PFASST requires at least two processes.";
      throw logic_error("two processes required for Two-Level-PFASST");
    }
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::run()
  {
    Controller<TransferT, CommT>::run();

    assert(this->get_communicator() != nullptr);

    if (this->get_num_steps() % this->get_communicator()->get_size() != 0) {
      CLOG(ERROR, this->get_logger_id()) << "Number of time steps (" << this->get_num_steps()
                                         << ") must be a multiple of the number of processors ("
                                         << this->get_communicator()->get_size() << ").";
      throw logic_error("number time steps must be multiple of number processors");
    }
    
    const size_t num_blocks = this->get_num_steps() / this->get_communicator()->get_size();

    if (num_blocks == 0) {
      CLOG(ERROR, this->get_logger_id()) << "Invalid Duration: There are more time processes ("
                                         << this->get_communicator()->get_size() << ") than time steps ("
                                         << this->get_num_steps() << ").";
      throw logic_error("invalid duration: too many time processes for given time steps");
    }

    do {
      this->status()->step() = this->_time_block * this->get_communicator()->get_size() + this->get_communicator()->get_rank();

      CLOG(INFO, this->get_logger_id()) << "";
      CLOG(INFO, this->get_logger_id()) << "Time Step " << (this->get_status()->get_step() + 1)
      << " of " << this->get_num_steps();

      this->status()->state() = State::PREDICTING;

      this->predictor();

      this->status()->iteration()++;

      this->status()->state() = State::ITERATING;

      if (!this->get_communicator()->is_last()) {
        CVLOG(1, this->get_logger_id()) << "sending status: " << to_string(this->get_status());
        this->get_status()->send(this->get_communicator(),
                                 this->get_communicator()->get_rank() + 1,
                                 this->compute_tag(1, true), true);
      }

      // iterate on each time step
      do {
        CLOG(INFO, this->get_logger_id()) << "";
        CLOG(INFO, this->get_logger_id()) << "Iteration " << this->get_status()->get_iteration();

        this->status()->state() = State::ITERATING;

        if (!this->get_communicator()->is_first()) {
          CVLOG(1, this->get_logger_id()) << "looking for state of previous process";
          this->_prev_status->recv(this->get_communicator(),
                                   this->get_communicator()->get_rank() - 1,
                                   this->compute_tag(1, true), true);
          CVLOG(1, this->get_logger_id()) << "Status received: " << to_string(this->_prev_status);
        }

        if (this->_prev_status->get_state() == State::FAILED) {
          CLOG(WARNING, this->get_logger_id()) << "previous process failed";
          CLOG(ERROR, this->get_logger_id()) << "We are aborting here. (read: TODO)";
          this->get_communicator()->abort(-1);
        } else if (this->_prev_status->get_state() == State::CONVERGED) {
          // TODO: do we really need to continue here ?!?
//          CLOG(WARNING, this->get_logger_id()) << "previous has converged and this needs to be propagated";
          continue;
        }

        this->cycle_down();

        if (!this->get_communicator()->is_first() && this->_prev_status->get_state() > State::FAILED) {
          this->get_coarse()->initial_state()->recv(this->get_communicator(),
                                                    this->get_communicator()->get_rank() - 1,
                                                    this->compute_tag(0, false), true);
        }
        this->sweep_coarse();

        if (!this->get_communicator()->is_last()) {
          this->get_coarse()->get_end_state()->send(this->get_communicator(),
                                                    this->get_communicator()->get_rank() + 1,
                                                    this->compute_tag(0, false), true);
        }

        this->cycle_up();

        this->sweep_fine();
        
        if (!this->get_communicator()->is_last()) {
          this->get_fine()->get_end_state()->send(this->get_communicator(),
                                                  this->get_communicator()->get_rank() + 1,
                                                  this->compute_tag(1, false), false);
        }

        // convergence check
      } while(this->advance_iteration());
      
      CLOG(INFO, this->get_logger_id()) << "Time Step done.";

      this->broadcast();
    } while(this->advance_time(this->get_communicator()->get_size()));
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelPfasst<TransferT, CommT>::advance_time(const size_t& num_steps)
  {
    CLOG(INFO, this->get_logger_id()) << "";

    if (TwoLevelMLSDC<TransferT, CommT>::advance_time(num_steps)) {
      this->_time_block++;
      return true;
    } else {
      return false;
    }
  }

  template<class TransferT, class CommT>
  bool
  TwoLevelPfasst<TransferT, CommT>::advance_iteration()
  {
    const bool fine_converged = this->get_fine()->converged();
    const bool previous_done = (this->get_communicator()->is_first())
                               ? true : this->_prev_status->get_state() <= State::FAILED;

    if (previous_done && fine_converged) {
      CLOG(INFO, this->get_logger_id()) << "FINE sweeper has converged as well as previous process.";
      this->status()->state() = State::CONVERGED;

    } else {
      if (Controller<TransferT, CommT>::advance_iteration()) {
        CLOG(INFO, this->get_logger_id()) << "FINE sweeper has not yet converged and additional iterations to do.";
        this->get_fine()->save();
        this->get_coarse()->save();
        this->status()->state() = State::ITERATING;

      } else {
        CLOG(WARNING, this->get_logger_id()) << "FINE sweeper has not yet converged and iterations threshold reached.";
        this->status()->state() = State::FAILED;

      }
    }

    if (!this->get_communicator()->is_last()) {
      CVLOG(1, this->get_logger_id()) << "sending status: " << to_string(this->get_status());
      this->get_status()->send(this->get_communicator(),
                               this->get_communicator()->get_rank() + 1,
                               this->compute_tag(1, true), true);
    }

    return (this->get_status()->get_state() > State::FAILED);
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::predict_coarse()
  {
    CLOG(INFO, this->get_logger_id()) << "Predicting on COARSE level";

    this->status()->state() = State::PRE_ITER_COARSE;
    this->get_coarse()->pre_predict();

    this->status()->state() = State::ITER_COARSE;
    this->get_coarse()->predict();

    this->status()->state() = State::POST_ITER_COARSE;
    this->get_coarse()->post_predict();

    this->status()->state() = State::PREDICTING;
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::predict_fine()
  {
    CLOG(INFO, this->get_logger_id()) << "Predicting on FINE level";

    this->status()->state() = State::PRE_ITER_FINE;
    this->get_fine()->pre_predict();

    this->status()->state() = State::ITER_FINE;
    this->get_fine()->predict();

    this->status()->state() = State::POST_ITER_FINE;
    this->get_fine()->post_predict();

    this->status()->state() = State::PREDICTING;
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::sweep_coarse()
  {
    CLOG(INFO, this->get_logger_id()) << "Sweeping on COARSE level";

    this->status()->state() = State::PRE_ITER_COARSE;
    this->get_coarse()->pre_sweep();

    this->status()->state() = State::ITER_COARSE;
    this->get_coarse()->sweep();

    this->status()->state() = State::POST_ITER_COARSE;
    this->get_coarse()->post_sweep();

    this->status()->state() = State::ITERATING;
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::sweep_fine()
  {
    CLOG(INFO, this->get_logger_id()) << "Sweeping on FINE level";

    this->status()->state() = State::PRE_ITER_FINE;
    this->get_fine()->pre_sweep();

    this->status()->state() = State::ITER_FINE;
    this->get_fine()->sweep();

    this->status()->state() = State::POST_ITER_FINE;
    this->get_fine()->post_sweep();

    this->status()->state() = State::ITERATING;
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::cycle_down()
  {
    CVLOG(1, this->get_logger_id()) << "cycle down to coarse level";

    this->get_transfer()->restrict(this->get_fine(), this->get_coarse(), true);
    this->get_transfer()->fas(this->get_status()->get_dt(), this->get_fine(), this->get_coarse());
    this->get_coarse()->save();
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::cycle_up()
  {
    CVLOG(1, this->get_logger_id()) << "cycle up to fine level";

    this->get_transfer()->interpolate(this->get_coarse(), this->get_fine(), true);

    if (!this->get_communicator()->is_first() && this->_prev_status->get_state() > State::FAILED) {
      assert(this->get_fine()->get_initial_state() != nullptr);
      CVLOG(1, this->get_logger_id()) << "looking for new initial value of fine level";
      this->get_fine()->initial_state()->recv(this->get_communicator(),
                                              this->get_communicator()->get_rank() - 1,
                                              this->compute_tag(1, false), false);
    }

    this->get_transfer()->interpolate_initial(this->get_coarse(), this->get_fine());
  }

  template<class TransferT, class CommT>
  void
  TwoLevelPfasst<TransferT, CommT>::predictor()
  {
    assert(this->get_status()->get_iteration() == 0);

    CLOG(INFO, this->get_logger_id()) << "";
    CLOG(INFO, this->get_logger_id()) << "PFASST Prediction step";

    // TODO: spread initial values to subsequent processes
    // restrict fine initial condition ...
    this->get_transfer()->restrict_initial(this->get_fine(), this->get_coarse());
    // ... and spread it to all nodes on the coarse level
    this->get_coarse()->spread();
    this->get_coarse()->save();

    // perform PFASST prediction sweeps on coarse level
    for (size_t predict_step = 0;
         predict_step <= this->get_communicator()->get_rank();
         ++predict_step) {
      // do the sweeper's prediction once ...
      if (predict_step == 0) {
        this->predict_coarse();
      } else {
        // and default sweeps for subsequent processes
        
        if (!this->get_communicator()->is_first()) {
          CVLOG(1, this->get_logger_id()) << "receiving coarse initial value";

          this->get_coarse()->initial_state()->recv(this->communicator(), this->get_communicator()->get_rank() - 1,
                                                    this->compute_tag(0, false), true);
        }

        this->sweep_coarse();
      }

    if (!this->get_communicator()->is_last()) {
      CVLOG(1, this->get_logger_id()) << "sending coarse end";
      this->get_coarse()->get_end_state()->send(this->communicator(), this->get_communicator()->get_rank() + 1,
                                                this->compute_tag(0, false), true);
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
    this->get_fine()->get_end_state()->bcast(this->get_communicator(), this->get_communicator()->get_size() - 1);
    
    // TODO: cleanup dangling request handlers of communicator here
  }

  template<class TransferT, class CommT>
  int
  TwoLevelPfasst<TransferT, CommT>::compute_tag(const size_t& level, const bool for_status) const
  {
    int tag = (level + 1) * 1000000;
    if (!for_status) {
      tag += this->get_status()->get_iteration() * 100;
    }
    CVLOG(2, this->get_logger_id()) << "tag for level " << level
                                    << " in iteration " << this->get_status()->get_iteration()
                                    << " for " << ((for_status) ? "status" : "data") << "communication"
                                    << " --> " << tag;
    return tag;
  }
}  // ::pfasst
