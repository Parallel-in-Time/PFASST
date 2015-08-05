#ifndef _PFASST__CONTROLLER__SDC_HPP_
#define _PFASST__CONTROLLER__SDC_HPP_

#include <memory>
using namespace std;

#include "pfasst/globals.hpp"
#include "pfasst/exceptions.hpp"
#include "pfasst/controller/status.hpp"
#include "pfasst/controller/controller.hpp"


namespace pfasst
{
  template<
    class TransferT
  >
  class SDC
    : public Controller<TransferT>
  {
    static_assert(is_same<
                    integral_constant<size_t, TransferT::traits::num_levels>,
                    integral_constant<size_t, 1>>::value,
                  "SDC only works for single sweeper setups.");

    public:
      typedef          TransferT                             transfer_type;
      typedef typename transfer_type::traits::fine_time_type time_type;

    protected:
      shared_ptr<typename transfer_type::traits::fine_sweeper_type> _sweeper;
      shared_ptr<void>                                              _coarse_level;
      shared_ptr<void>                                              _fine_level;

      shared_ptr<void>               _transfer;
      shared_ptr<Status<time_type>>  _status;
      shared_ptr<comm::Communicator> _comm;
      bool                           _ready;

    public:
      SDC();
      SDC(const SDC<TransferT>& other) = default;
      SDC(SDC<TransferT>&& other) = default;
      virtual ~SDC() = default;
      SDC<TransferT>& operator=(const SDC<TransferT>& other) = default;
      SDC<TransferT>& operator=(SDC<TransferT>&& other) = default;

      virtual size_t get_num_levels() const override;

      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse);
      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper);

      virtual void add_transfer(shared_ptr<TransferT> transfer) override;

      virtual const shared_ptr<typename TransferT::traits::fine_sweeper_type> get_sweeper() const;
      virtual       shared_ptr<typename TransferT::traits::fine_sweeper_type> get_sweeper();

      virtual void setup() override;
      virtual void run() override;

      virtual bool advance_time(const size_t& num_steps = 1) override;
      virtual bool advance_iteration() override;
  };
}  // ::pfasst

#include "pfasst/controller/sdc_impl.hpp"

#endif  // _PFASST__CONTROLLER__INTERFACE_HPP_
