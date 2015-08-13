#ifndef _PFASST__CONTROLLER__TWO_LEVEL_MLSDC_HPP_
#define _PFASST__CONTROLLER__TWO_LEVEL_MLSDC_HPP_

#include <memory>
using namespace std;

#include "pfasst/controller/controller.hpp"


namespace pfasst
{
  template<
    class TransferT
  >
  class TwoLevelMLSDC
    : public Controller<TransferT>
  {
    public:
      typedef          TransferT                             transfer_type;
      typedef typename transfer_type::traits::fine_time_type time_type;

      static void init_loggers();

    protected:
      shared_ptr<typename transfer_type::traits::coarse_sweeper_type> _coarse_level;
      shared_ptr<typename transfer_type::traits::fine_sweeper_type>   _fine_level;

      virtual void predict_coarse();
      virtual void predict_fine();
      virtual void sweep_coarse();
      virtual void sweep_fine();

      virtual void cycle_down();
      virtual void cycle_up();

    public:
      TwoLevelMLSDC();
      TwoLevelMLSDC(const TwoLevelMLSDC<TransferT>& other) = default;
      TwoLevelMLSDC(TwoLevelMLSDC<TransferT>&& other) = default;
      virtual ~TwoLevelMLSDC() = default;
      TwoLevelMLSDC<TransferT>& operator=(const TwoLevelMLSDC<TransferT>& other) = default;
      TwoLevelMLSDC<TransferT>& operator=(TwoLevelMLSDC<TransferT>&& other) = default;

      virtual size_t get_num_levels() const override;

      template<class SweeperT>
      void add_sweeper(shared_ptr<SweeperT> sweeper, const bool as_coarse);

      virtual const shared_ptr<typename TransferT::traits::coarse_sweeper_type> get_coarse() const;
      virtual       shared_ptr<typename TransferT::traits::coarse_sweeper_type> get_coarse();
      virtual const shared_ptr<typename TransferT::traits::fine_sweeper_type> get_fine() const;
      virtual       shared_ptr<typename TransferT::traits::fine_sweeper_type> get_fine();

      virtual void set_options() override;

      virtual void setup() override;
      virtual void run() override;

      virtual bool advance_time(const size_t& num_steps = 1) override;
      virtual bool advance_iteration() override;
  };
}  // ::pfasst

#include "pfasst/controller/two_level_mlsdc_impl.hpp"

#endif  // _PFASST__CONTROLLER__TWO_LEVEL_MLSDC_HPP_
