#ifndef _PFASST__CONTROLLER__TWO_LEVEL_PFASST_HPP_
#define _PFASST__CONTROLLER__TWO_LEVEL_PFASST_HPP_

#include <memory>
using namespace std;

#include "pfasst/controller/two_level_mlsdc.hpp"
#include "pfasst/comm/communicator.hpp"


namespace pfasst
{
  template<
    class TransferT
  >
  class TwoLevelPfasst
    : public TwoLevelMLSDC<TransferT>
  {
    public:
      typedef          TransferT                             transfer_type;
      typedef typename transfer_type::traits::fine_time_type time_type;

      static void init_loggers() override;

    protected:
      shared_ptr<comm::Communicator> _comm;

      virtual void predictor() override;
      virtual void cycle_down() override;
      virtual void cycle_up() override;

      virtual void broadcast();
      virtual int compute_tag() const;

    public:
      TwoLevelPfasst();
      TwoLevelPfasst(const TwoLevelPfasst<TransferT>& other) = default;
      TwoLevelPfasst(TwoLevelPfasst<TransferT>&& other) = default;
      virtual ~TwoLevelPfasst() = default;
      TwoLevelPfasst<TransferT>& operator=(const TwoLevelPfasst<TransferT>& other) = default;
      TwoLevelPfasst<TransferT>& operator=(TwoLevelPfasst<TransferT>&& other) = default;

      virtual       shared_ptr<comm::Communicator>& communicator();
      virtual const shared_ptr<comm::Communicator>  get_communicator() const;

      virtual void set_options() override;

      virtual void setup() override;
      virtual void run() override;

      virtual bool advance_time(const size_t& num_steps = 1) override;
      virtual bool advance_iteration() override;
  };
}  // ::pfasst

#include "pfasst/controller/two_level_pfasst_impl.hpp"

#endif  // _PFASST__CONTROLLER__TWO_LEVEL_PFASST_HPP_
