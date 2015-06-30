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
  const shared_ptr<typename TransferT::traits::coarse_sweeper_type>
  SDC<TransferT>::get_coarse() const
  {
    CLOG(WARNING, "CONTROL") << "SDC Controller has no Coarse Level.";
    return this->get_sweeper();
  }

  template<class TransferT>
  shared_ptr<typename TransferT::traits::coarse_sweeper_type>
  SDC<TransferT>::get_coarse()
  {
    CLOG(WARNING, "CONTROL") << "SDC Controller has no Coarse Level.";
    return this->get_sweeper();
  }

  template<class TransferT>
  const shared_ptr<typename TransferT::traits::fine_sweeper_type>
  SDC<TransferT>::get_fine() const
  {
    CLOG(WARNING, "CONTROL") << "SDC Controller has no Fine Level.";
    return this->get_sweeper();
  }

  template<class TransferT>
  shared_ptr<typename TransferT::traits::fine_sweeper_type>
  SDC<TransferT>::get_fine()
  {
    CLOG(WARNING, "CONTROL") << "SDC Controller has no Fine Level.";
    return this->get_sweeper();
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
  }

  template<class TransferT>
  void
  SDC<TransferT>::run()
  {
    Controller<TransferT>::run();

    // iterate over time steps
    for(; this->get_status()->get_time() < this->get_status()->get_t_end(); this->advance_time()) {
      const bool do_initial = this->get_status()->get_step() == 0;

      // iterate on current time step
      for (
        this->status()->iteration() = 0;
        this->get_status()->get_iteration() < this->get_status()->get_max_iterations();
      ) {
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

        if (this->get_sweeper()->converged()) {
          break;
        }

        if (this->advance_iteration()) {
          this->get_sweeper()->save();
        }
      }

      this->get_sweeper()->post_step();

      if (this->advance_time()) {
        this->get_sweeper()->advance();
      }
    }
  }

}  // ::pfasst
