#include "pfasst/controller/controller.hpp"

#include <cmath>
using namespace std;

#include "pfasst/util.hpp"
#include "pfasst/config.hpp"
#include "pfasst/logging.hpp"


namespace pfasst
{
  template<class TransferT, class CommT>
  Controller<TransferT, CommT>::Controller()
    :   _status(make_shared<Status<typename TransferT::traits::fine_time_type>>())
      , _ready(false)
      , _logger_id("CONTROL")
  {}

  template<class TransferT, class CommT>
  shared_ptr<CommT>&
  Controller<TransferT, CommT>::communicator()
  {
    return this->_comm;
  }

  template<class TransferT, class CommT>
  const shared_ptr<CommT>
  Controller<TransferT, CommT>::get_communicator() const
  {
    return this->_comm;
  }

  template<class TransferT, class CommT>
  shared_ptr<Status<typename TransferT::traits::fine_time_type>>&
  Controller<TransferT, CommT>::status()
  {
    return this->_status;
  }

  template<class TransferT, class CommT>
  const shared_ptr<Status<typename TransferT::traits::fine_time_type>>
  Controller<TransferT, CommT>::get_status() const
  {
    return this->_status;
  }

  template<class TransferT, class CommT>
  size_t
  Controller<TransferT, CommT>::get_num_levels() const
  {
    return 0;
  }

  template<class TransferT, class CommT>
  size_t
  Controller<TransferT, CommT>::get_num_steps() const
  {
    if (this->get_status()->get_t_end() <= 0) {
      CLOG(ERROR, this->get_logger_id()) << "Time end point must be non-zero positive."
        << " NOT " << this->get_status()->get_t_end();
      throw logic_error("time end point must be non-zero positive");
    }

    if (this->get_status()->get_dt() <= 0.0) {
      CLOG(ERROR, this->get_logger_id()) << "Time delta must be non-zero positive."
        << " NOT " << this->get_status()->get_dt();
      throw logic_error("time delta must be non-zero positive");
    }

    const auto div = this->get_status()->get_t_end() / this->get_status()->get_dt();
    CLOG_IF(!almost_equal(div * this->get_status()->get_dt(),
                          (size_t)div * this->get_status()->get_dt()), WARNING, this->get_logger_id())
      << "End time point not an integral multiple of time delta: "
      << this->get_status()->get_t_end() << " / " << this->get_status()->get_dt()
      << " = " << div << " != " << (size_t)div;

    return (size_t)(this->get_status()->get_t_end() / this->get_status()->get_dt());
  }

  template<class TransferT, class CommT>
  bool&
  Controller<TransferT, CommT>::ready()
  {
    return this->_ready;
  }

  template<class TransferT, class CommT>
  bool
  Controller<TransferT, CommT>::is_ready() const
  {
    return this->_ready;
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::set_logger_id(const string& logger_id)
  {
    this->_logger_id = logger_id;
  }

  template<class TransferT, class CommT>
  const char*
  Controller<TransferT, CommT>::get_logger_id() const
  {
    return this->_logger_id.c_str();
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::set_options()
  {
    this->status()->max_iterations() = config::get_value<size_t>("num_iters", this->get_status()->get_max_iterations());
    this->status()->t_end() = config::get_value<typename TransferT::traits::fine_time_type>("t_end", this->get_status()->get_t_end());
  }

  template<class TransferT, class CommT>
  template<class SweeperT>
  void
  Controller<TransferT, CommT>::add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse)
  {}

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::add_transfer(shared_ptr<TransferT> transfer)
  {
    this->_transfer = transfer;
  }

  template<class TransferT, class CommT>
  const shared_ptr<TransferT>
  Controller<TransferT, CommT>::get_transfer() const
  {
    return this->_transfer;
  }

  template<class TransferT, class CommT>
  shared_ptr<TransferT>
  Controller<TransferT, CommT>::get_transfer()
  {
    return this->_transfer;
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::setup()
  {
    CLOG_IF(this->is_ready(), WARNING, this->get_logger_id())
      << "Controller has already been setup.";

    CVLOG(1, this->get_logger_id()) << "setting up controller";

    if (this->get_status()->get_t_end() <= 0.0) {
      CLOG(ERROR, this->get_logger_id()) << "End time point must be larger than zero."
        << " (" << this->get_status()->get_t_end() << ")";
      throw logic_error("end time point must be larger zero");
    }

    const auto num_steps = this->get_num_steps();
    if (num_steps * this->get_status()->get_dt() != this->get_status()->get_t_end()) {
      CLOG(ERROR, this->get_logger_id()) << "End time point not an integral multiple of time delta. "
        << " (" << num_steps << " * " << this->get_status()->get_dt()
        << " = " << num_steps * this->get_status()->get_dt() << " != " << this->get_status()->get_t_end() << ")";
      throw logic_error("time end point is not an integral multiple of time delta");
    }

    CLOG_IF(this->get_status()->get_max_iterations() == 0, WARNING, this->get_logger_id())
      << "You sould define a maximum number of iterations to avoid endless runs."
      << " (" << this->get_status()->get_max_iterations() << ")";

    this->ready() = true;
  }

  template<class TransferT, class CommT>
  void
  Controller<TransferT, CommT>::run()
  {
    if (!this->is_ready()) {
      CLOG(ERROR, this->get_logger_id()) << "Controller is not ready to run. setup() not called yet.";
      throw logic_error("controller not ready to run");
    }
  }

  template<class TransferT, class CommT>
  bool
  Controller<TransferT, CommT>::advance_time(const size_t& num_steps)
  {
    const time_type delta_time = num_steps * this->get_status()->get_dt();
    const time_type new_time = this->get_status()->get_time() + delta_time;


    if (new_time > this->get_status()->get_t_end()) {
      CLOG(WARNING, this->get_logger_id()) << "Not advancing " << num_steps
                                           << ((num_steps > 1) ? " time steps " : " time step ")
                                           << "with dt=" << this->get_status()->get_dt() << " to t=" << new_time
                                           << " as it will exceed T_end=" << this->get_status()->get_t_end() << " by "
                                           << (new_time - this->get_status()->get_t_end());

      return false;

    } else if(almost_equal(new_time, this->get_status()->get_t_end())) {
      CLOG(INFO, this->get_logger_id()) << "End time point reached: " << this->get_status()->get_t_end();

      return false;

    } else {
      CLOG(INFO, this->get_logger_id()) << "Advancing " << num_steps
                                        << ((num_steps > 1) ? " time steps " : " time step ")
                                        << "with dt=" << this->get_status()->get_dt() << " to t=" << new_time;

      this->status()->time() += delta_time;
      this->status()->step() += num_steps;
      this->status()->iteration() = 0;

      return true;
    }
  }

  template<class TransferT, class CommT>
  bool
  Controller<TransferT, CommT>::advance_iteration()
  {
    if (this->get_status()->get_iteration() + 1 > this->get_status()->get_max_iterations()) {
      CLOG(WARNING, this->get_logger_id()) << "Not advancing to next iteration ("
                                           << (this->get_status()->get_iteration() + 1)
                                           << ") as it will exceed maximum number of allowed iterations ("
                                           << this->get_status()->get_max_iterations() << ")";

      return false;

    } else {
      CLOG(INFO, this->get_logger_id()) << "Advancing to next iteration -> "
                                        << (this->get_status()->get_iteration() + 1);

      this->status()->iteration()++;

      return true;
    }
  }
}  // ::pfasst
