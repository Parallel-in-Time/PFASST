#include "fixtures/test_helpers.hpp"

#include <pfasst/controller/two_level_mlsdc.hpp>
using pfasst::TwoLevelMLSDC;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>

#include <pfasst/transfer/traits.hpp>

#include "comm/mocks.hpp"
#include "sweeper/mocks.hpp"
#include "transfer/mocks.hpp"

typedef pfasst::vector_encap_traits<double, double>                     VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait>                  VectorEncapsulation;
typedef NiceMock<SweeperMock<pfasst::sweeper_traits<VectorEncapTrait>>> SweeperType;
typedef pfasst::transfer_traits<SweeperType, SweeperType, 2>            TransferTraits;
typedef NiceMock<TransferMock<TransferTraits>>                          TransferType;
typedef NiceMock<CommMock>                                              CommunicatorType;


typedef ::testing::Types<TwoLevelMLSDC<TransferType>> ControllerTypes;
INSTANTIATE_TYPED_TEST_CASE_P(TwoLevelMLSDC, Concepts, ControllerTypes);


class Interface
  : public ::testing::Test
{
  protected:
    shared_ptr<TwoLevelMLSDC<TransferType>> controller;

    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<CommunicatorType> comm;

    virtual void SetUp()
    {
      this->controller = make_shared<TwoLevelMLSDC<TransferType>>();
      this->status = make_shared<pfasst::Status<double>>();
      this->comm = make_shared<CommunicatorType>();
    }
};

TEST_F(Interface, has_a_status)
{
  EXPECT_THAT(controller->get_status(), NotNull());
}

TEST_F(Interface, status_can_be_assigned)
{
  controller->status() = status;
  EXPECT_THAT(controller->get_status(), Eq(status));
}

TEST_F(Interface, status_can_be_modified)
{
  controller->status()->time() = 42.0;
  EXPECT_THAT(controller->get_status()->get_time(), Eq(42.0));
}

TEST_F(Interface, has_no_communicator_after_instantiation)
{
  EXPECT_THAT(controller->get_communicator(), IsNull());
}

TEST_F(Interface, communicator_can_be_assigned)
{
  ASSERT_THAT(controller->get_communicator(), Not(Eq(comm)));

  controller->communicator() = comm;
  EXPECT_THAT(controller->get_communicator(), Eq(comm));
}


class Setup
  : public ::testing::Test
{
  protected:
    shared_ptr<TwoLevelMLSDC<TransferType>> controller;

    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<SweeperType> sweeper1;
    shared_ptr<SweeperType> sweeper2;
    shared_ptr<TransferType> transfer;
    shared_ptr<VectorEncapsulation> sweeper1_initial;
    shared_ptr<VectorEncapsulation> sweeper1_end;
    shared_ptr<VectorEncapsulation> sweeper2_initial;
    shared_ptr<VectorEncapsulation> sweeper2_end;

    virtual void SetUp()
    {
      this->controller = make_shared<TwoLevelMLSDC<TransferType>>();
      this->transfer = make_shared<TransferType>();
      this->status = make_shared<pfasst::Status<double>>();

      this->sweeper1 = make_shared<SweeperType>();
      this->sweeper2 = make_shared<SweeperType>();

      this->sweeper1_initial = this->sweeper1->get_encap_factory()->create();
      this->sweeper1_end = this->sweeper1->get_encap_factory()->create();
      this->sweeper2_initial = this->sweeper2->get_encap_factory()->create();
      this->sweeper2_end = this->sweeper2->get_encap_factory()->create();

      ON_CALL(*(this->sweeper1.get()), get_initial_state())
        .WillByDefault(Return(this->sweeper1_initial));
      ON_CALL(*(this->sweeper1.get()), initial_state())
        .WillByDefault(ReturnRef(this->sweeper1_initial));
      ON_CALL(*(this->sweeper1.get()), get_end_state())
        .WillByDefault(Return(this->sweeper1_end));

      ON_CALL(*(this->sweeper2.get()), get_initial_state())
        .WillByDefault(Return(this->sweeper2_initial));
      ON_CALL(*(this->sweeper2.get()), initial_state())
        .WillByDefault(ReturnRef(this->sweeper2_initial));
      ON_CALL(*(this->sweeper2.get()), get_end_state())
        .WillByDefault(Return(this->sweeper2_end));
    }
};

TEST_F(Setup, adding_coarser_level)
{
  ASSERT_THAT(controller->get_num_levels(), Eq(0));

  controller->add_sweeper(sweeper1, false);
  EXPECT_THAT(controller->get_num_levels(), Eq(1));

  controller->add_sweeper(sweeper2, true);
  EXPECT_THAT(controller->get_num_levels(), Eq(2));
}

TEST_F(Setup, adding_finer_level)
{
  ASSERT_THAT(controller->get_num_levels(), Eq(0));

  controller->add_sweeper(sweeper1, true);
  EXPECT_THAT(controller->get_num_levels(), Eq(1));

  controller->add_sweeper(sweeper2, false);
  EXPECT_THAT(controller->get_num_levels(), Eq(2));
}

TEST_F(Setup, exactly_two_levels_must_be_added)
{
  controller->status() = status;
  controller->status()->t_end() = 0.1;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;
  controller->add_transfer(transfer);

  ASSERT_THAT(controller->get_num_levels(), Eq(0));
  EXPECT_THROW(controller->setup(), logic_error);

  controller = make_shared<TwoLevelMLSDC<TransferType>>();
  controller->status() = status;
  controller->status()->t_end() = 0.1;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;
  controller->add_transfer(transfer);
  controller->add_sweeper(sweeper1, true);
  ASSERT_THAT(controller->get_num_levels(), Eq(1));
  EXPECT_THROW(controller->setup(), logic_error);

  controller = make_shared<TwoLevelMLSDC<TransferType>>();
  controller->status() = status;
  controller->status()->t_end() = 0.1;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;
  controller->add_transfer(transfer);
  controller->add_sweeper(sweeper1, true);
  controller->add_sweeper(sweeper2, false);
  EXPECT_CALL(*(sweeper1.get()), status()).Times(AnyNumber()).WillRepeatedly(ReturnRef(status));
  EXPECT_CALL(*(sweeper2.get()), status()).Times(AnyNumber()).WillRepeatedly(ReturnRef(status));
  EXPECT_CALL(*(sweeper1.get()), setup()).Times(1);
  EXPECT_CALL(*(sweeper2.get()), setup()).Times(1);
  ASSERT_THAT(controller->get_num_levels(), Eq(2));
  controller->setup();
}

TEST_F(Setup, setup_required_for_running)
{
  controller->status()->t_end() = 0.1;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;
  controller->add_sweeper(sweeper1, true);
  controller->add_sweeper(sweeper1, false);
  controller->add_transfer(transfer);

  ASSERT_FALSE(controller->is_ready());
  EXPECT_THROW(controller->run(), logic_error);

  EXPECT_CALL(*(sweeper1.get()), status()).Times(AnyNumber()).WillRepeatedly(ReturnRef(status));
  EXPECT_CALL(*(sweeper2.get()), status()).Times(AnyNumber()).WillRepeatedly(ReturnRef(status));
  EXPECT_CALL(*(sweeper1.get()), setup()).Times(AnyNumber());
  EXPECT_CALL(*(sweeper2.get()), setup()).Times(AnyNumber());

  controller->setup();
  EXPECT_TRUE(controller->is_ready());
  controller->run();
}


class Logic
  : public ::testing::Test
{
  protected:
    shared_ptr<TwoLevelMLSDC<TransferType>> controller;

    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<SweeperType> sweeper1;
    shared_ptr<SweeperType> sweeper2;
    shared_ptr<TransferType> transfer;

    virtual void SetUp()
    {
      this->controller = make_shared<TwoLevelMLSDC<TransferType>>();
      this->transfer = make_shared<TransferType>();
      this->status = make_shared<pfasst::Status<double>>();
      this->sweeper1 = make_shared<SweeperType>();
      this->sweeper2 = make_shared<SweeperType>();
      this->controller->add_sweeper(this->sweeper1, true);
      this->controller->add_sweeper(this->sweeper2, false);
      this->controller->add_transfer(this->transfer);
    }
};

TEST_F(Logic, advance_in_time_with_sufficient_t_end)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 1.2;
  
  EXPECT_CALL(*(sweeper1.get()), advance(1)).Times(1);
  EXPECT_CALL(*(sweeper2.get()), advance(1)).Times(1);

  EXPECT_TRUE(controller->advance_time(1));
  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.1));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(2));
}

TEST_F(Logic, advance_in_time_with_insufficient_t_end)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 1.1;

  EXPECT_FALSE(controller->advance_time());
  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.0));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(1));
}

TEST_F(Logic, advance_in_time_multiple_steps_at_once)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 1.4;

  EXPECT_TRUE(controller->advance_time(3));
  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.3));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(4));
}


TEST_F(Logic, advance_iteration_with_exceeding_max_iteration_threshold)
{
  controller->status()->iteration() = 1;
  controller->status()->max_iterations() = 1;
  ASSERT_THAT(controller->get_status()->get_iteration(), Eq(1));
  ASSERT_THAT(controller->get_status()->get_max_iterations(), Eq(1));

  EXPECT_FALSE(controller->advance_iteration());
  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(1));
}

TEST_F(Logic, advance_iteration)
{
  controller->status()->iteration() = 1;
  controller->status()->max_iterations() = 5;
  ASSERT_THAT(controller->get_status()->get_iteration(), Eq(1));
  ASSERT_THAT(controller->get_status()->get_max_iterations(), Eq(5));

  EXPECT_TRUE(controller->advance_iteration());
  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(2));
}


TEST_MAIN()
