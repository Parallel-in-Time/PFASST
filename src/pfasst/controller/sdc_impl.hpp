#include "pfasst/controller/sdc.hpp"

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT>
  SDC<TransferT>::SDC()
    : Controller<TransferT>()
  {}

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
    CLOG(WARNING, "CONTROL") << "SDC Controller does not require a transfer operator.";
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
  SDC<TransferT>::setup()
  {
    Controller<TransferT>::setup();

    if (this->get_num_levels() != 1) {
      CLOG(ERROR, "CONTROL") << "One level (Sweeper) must have been added for SDC.";
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

    // iterate over time steps
    do {
      const bool do_initial = this->get_status()->get_step() == 0;

      // iterate on current time step
      do {
        const bool do_prediction = this->get_status()->get_iteration() == 0;

        if (do_prediction) {
          this->get_sweeper()->pre_predict();
          this->get_sweeper()->predict();
          this->get_sweeper()->post_predict();
        } else {
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
      this->get_sweeper()->advance();
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
      return false;
    } else if (Controller<TransferT>::advance_iteration()) {
      this->get_sweeper()->save();
      return true;
    } else {
      return false;
    }
  }
}  // ::pfasst
