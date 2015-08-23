#include "pfasst/controller/sdc.hpp"

#include <stdexcept>
using namespace std;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT>
  SDC<TransferT>::SDC()
    : Controller<TransferT>()
  {
    SDC<TransferT>::init_loggers();
    this->set_logger_id("SDC");
  }

  template<class TransferT>
  void
  SDC<TransferT>::init_loggers()
  {
    log::add_custom_logger("SDC");
  }

  template<class TransferT>
  size_t
  SDC<TransferT>::get_num_levels() const
  {
    return ((this->_sweeper != nullptr) ? 1 : 0);
  }

  template<class TransferT>
  template<class SweeperT>
  void
  SDC<TransferT>::add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse)
  {
    this->add_sweeper(sweeper);
  }

  template<class TransferT>
  template<class SweeperT>
  void
  SDC<TransferT>::add_sweeper(shared_ptr<SweeperT> sweeper)
  {
    static_assert(is_same<SweeperT, typename TransferT::traits::fine_sweeper_type>::value,
                  "Sweeper must be a Fine Sweeper Type.");

    this->_sweeper = sweeper;
  }

  template<class TransferT>
  void
  SDC<TransferT>::add_transfer(shared_ptr<TransferT> transfer)
  {
    CLOG(WARNING, this->get_logger_id()) << "SDC Controller does not require a transfer operator.";
  }

  template<class TransferT>
  const shared_ptr<typename TransferT::traits::fine_sweeper_type>
  SDC<TransferT>::get_sweeper() const
  {
    return this->_sweeper;
  }

  template<class TransferT>
  shared_ptr<typename TransferT::traits::fine_sweeper_type>
  SDC<TransferT>::get_sweeper()
  {
    return this->_sweeper;
  }

  template<class TransferT>
  void
  SDC<TransferT>::set_options()
  {
    Controller<TransferT>::set_options();

    this->get_sweeper()->set_options();
  }

  template<class TransferT>
  void
  SDC<TransferT>::setup()
  {
    Controller<TransferT>::setup();

    if (this->get_num_levels() != 1) {
      CLOG(ERROR, this->get_logger_id()) << "One level (Sweeper) must have been added for SDC.";
      throw logic_error("SDC requires one level");
    }

    this->get_sweeper()->status() = this->get_status();
    this->get_sweeper()->setup();
  }

  template<class TransferT>
  void
  SDC<TransferT>::run()
  {
    Controller<TransferT>::run();

    CLOG(INFO, this->get_logger_id()) << "";
    CLOG(INFO, this->get_logger_id()) << "Sequential SDC";
    CLOG(INFO, this->get_logger_id()) << "  t0:        " << LOG_FIXED << this->get_status()->get_time();
    CLOG(INFO, this->get_logger_id()) << "  dt:        " << LOG_FIXED << this->get_status()->get_dt();
    CLOG(INFO, this->get_logger_id()) << "  T:         " << LOG_FIXED << this->get_status()->get_t_end();
    CLOG(INFO, this->get_logger_id()) << "  num steps: " << LOG_FIXED << this->get_num_steps();
    CLOG(INFO, this->get_logger_id()) << "  max iter:  " << LOG_FIXED << this->get_status()->get_max_iterations();
    CLOG(INFO, this->get_logger_id()) << "  Initial Value: " << to_string(this->get_sweeper()->get_initial_state());

    // iterate over time steps
    do {
      const bool do_initial = this->get_status()->get_step() == 0;
      CLOG(INFO, this->get_logger_id()) << "";
      CLOG(INFO, this->get_logger_id()) << "Time Step " << (this->get_status()->get_step() + 1)
                                        << " of " << this->get_num_steps();

      // iterate on current time step
      do {
        const bool do_prediction = this->get_status()->get_iteration() == 0;

        if (do_prediction) {
          CLOG(INFO, this->get_logger_id()) << "";
          CLOG(INFO, this->get_logger_id()) << "SDC Prediction step";
          this->get_sweeper()->pre_predict();
          this->get_sweeper()->predict();
          this->get_sweeper()->post_predict();
        } else {
          CLOG(INFO, this->get_logger_id()) << "";
          CLOG(INFO, this->get_logger_id()) << "Iteration " << this->get_status()->get_iteration();
          this->get_sweeper()->pre_sweep();
          this->get_sweeper()->sweep();
          this->get_sweeper()->post_sweep();
        }
      } while(this->advance_iteration());
    } while(this->advance_time());
  }

  template<class TransferT>
  bool
  SDC<TransferT>::advance_time(const size_t& num_steps)
  {
    this->get_sweeper()->post_step();

    if (Controller<TransferT>::advance_time(num_steps)) {
      this->get_sweeper()->advance(num_steps);
      return true;
    } else {
      return false;
    }
  }

  template<class TransferT>
  bool
  SDC<TransferT>::advance_iteration()
  {
    if (this->get_sweeper()->converged()) {
      CLOG(INFO, this->get_logger_id()) << "Sweeper has converged.";
      return false;
    } else if (Controller<TransferT>::advance_iteration()) {
      CLOG(INFO, this->get_logger_id()) << "Sweeper has not yet converged and additional iterations to do.";
      this->get_sweeper()->save();
      return true;
    } else {
      CLOG(INFO, this->get_logger_id()) << "Sweeper has not yet converged and no more iterations to do.";
      return false;
    }
  }
}  // ::pfasst
