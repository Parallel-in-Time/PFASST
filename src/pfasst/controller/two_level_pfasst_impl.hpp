#include "pfasst/controller/two_level_pfasst.hpp"

#include <cassert>
using namespace std;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT>
  TwoLevelPfasst<TransferT>::TwoLevelPfasst()
    : Controller<TransferT>()
  {}

  template<class TransferT>
  shared_ptr<comm::Communicator>&
  TwoLevelPfasst<TransferT>::communicator()
  {
    return this->_comm;
  }

  template<class TransferT>
  const shared_ptr<comm::Communicator>
  TwoLevelPfasst<TransferT>::get_communicator() const
  {
    return this->_comm;
  }

  template<class TransferT>
  size_t
  TwoLevelPfasst<TransferT>::get_num_levels() const
  {
    size_t num = 0;
    if (this->_coarse_level != nullptr) {
      num++;
    }
    if (this->_fine_level != nullptr) {
      num++;
    }
    return num;
  }

  template<class TransferT>
  template<class SweeperT>
  void
  TwoLevelPfasst<TransferT>::add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse)
  {
    static_assert(is_same<SweeperT, typename TransferT::traits::fine_sweeper_type>::value
                  || is_same<SweeperT, typename TransferT::traits::coarse_sweeper_type>::value,
                  "Sweeper must be either a Coarse or Fine Sweeper Type.");

    if (as_coarse) {
      if (is_same<SweeperT, typename transfer_type::traits::coarse_sweeper_type>::value) {
        this->_coarse_level = sweeper;
      } else {
        CLOG(ERROR, "CONTROL") << "Type of given Sweeper ("
          << typeid(SweeperT).name() << ") is not applicable as Coarse Sweeper ("
          << typeid(typename transfer_type::traits::coarse_sweeper_type).name() << ").";
        throw logic_error("given sweeper can not be used as coarse sweeper");
      }
    } else {
      if (is_same<SweeperT, typename transfer_type::traits::fine_sweeper_type>::value) {
        this->_fine_level = sweeper;
      } else {
        CLOG(ERROR, "CONTROL") << "Type of given Sweeper ("
          << typeid(SweeperT).name() << ") is not applicable as Fine Sweeper ("
          << typeid(typename transfer_type::traits::fine_sweeper_type).name() << ").";
        throw logic_error("given sweeper can not be used as fine sweeper");
      }
    }
  }

  template<class TransferT>
  const shared_ptr<typename TransferT::traits::coarse_sweeper_type>
  TwoLevelPfasst<TransferT>::get_coarse() const
  {
    return this->_coarse_level;
  }

  template<class TransferT>
  shared_ptr<typename TransferT::traits::coarse_sweeper_type>
  TwoLevelPfasst<TransferT>::get_coarse()
  {
    return this->_coarse_level;
  }

  template<class TransferT>
  const shared_ptr<typename TransferT::traits::fine_sweeper_type>
  TwoLevelPfasst<TransferT>::get_fine() const
  {
    return this->_fine_level;
  }

  template<class TransferT>
  shared_ptr<typename TransferT::traits::fine_sweeper_type>
  TwoLevelPfasst<TransferT>::get_fine()
  {
    return this->_fine_level;
  }

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::set_options()
  {
    Controller<TransferT>::set_options();
  }


  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::setup()
  {
    assert(this->get_communicator() != nullptr);
    assert(this->get_transfer() != nullptr);

    Controller<TransferT>::setup();

    if (this->get_num_levels() != 2) {
      CLOG(ERROR, "CONTROL") << "Two levels (Sweeper) must have been added for Two-Level-PFASST.";
      throw logic_error("Two-Level-PFASST requires two levels");
    }

    if (this->get_communicator()->get_size() < 2) {
      CLOG(ERROR, "CONTROL") << "Two-Level-PFASST requires at least two processes.";
      throw logic_error("two processes required for Two-Level-PFASST");
    }
  }

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::run()
  {
    Controller<TransferT>::run();

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

  template<class TransferT>
  bool
  TwoLevelPfasst<TransferT>::advance_time(const size_t& num_steps)
  {
    return Controller<TransferT>::advance_time(num_steps);
  }

  template<class TransferT>
  bool
  TwoLevelPfasst<TransferT>::advance_iteration()
  {
    return Controller<TransferT>::advance_iteration();
  }


  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::predict_coarse()
  {
    this->status()->state() = State::PRE_ITER_COARSE;
    this->get_coarse()->pre_predict();

    this->status()->state() = State::ITER_COARSE;
    this->get_coarse()->predict();

    this->status()->state() = State::POST_ITER_COARSE;
    this->get_coarse()->post_predict();

    this->status()->state() = State::PREDICTING;
  }

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::predict_fine()
  {
    this->status()->state() = State::PRE_ITER_FINE;
    this->get_fine()->pre_predict();

    this->status()->state() = State::ITER_FINE;
    this->get_fine()->predict();

    this->status()->state() = State::POST_ITER_FINE;
    this->get_fine()->post_predict();

    this->status()->state() = State::PREDICTING;
  }

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::sweep_coarse()
  {
    this->status()->state() = State::PRE_ITER_COARSE;
    this->get_coarse()->pre_sweep();

    this->status()->state() = State::ITER_COARSE;
    this->get_coarse()->sweep();

    this->status()->state() = State::POST_ITER_COARSE;
    this->get_coarse()->post_sweep();

    this->status()->state() = State::ITERATING;
  }

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::sweep_fine()
  {
    this->status()->state() = State::PRE_ITER_FINE;
    this->get_fine()->pre_sweep();

    this->status()->state() = State::ITER_FINE;
    this->get_fine()->sweep();

    this->status()->state() = State::POST_ITER_FINE;
    this->get_fine()->post_sweep();

    this->status()->state() = State::ITERATING;
  }

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::cycle_down()
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

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::cycle_up()
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

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::predictor()
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

  template<class TransferT>
  void
  TwoLevelPfasst<TransferT>::broadcast()
  {
    this->get_fine()->get_end_state()->bcast(this->get_communicator(),
                                             this->get_communicator()->get_root());
  }

  template<class TransferT>
  int
  TwoLevelPfasst<TransferT>::compute_tag() const
  {
    // TODO: come up with a good way of computing unique but meaningful tags
    return 0;
  }
}  // ::pfasst
