#include "fixtures/test_helpers.hpp"

#include <pfasst/controller/sdc.hpp>
using pfasst::SDC;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>

#include <pfasst/transfer/traits.hpp>
#include <pfasst/transfer/polynomial.hpp>

#include "sweeper/mocks.hpp"
#include "quadrature/mocks.hpp"
#include "transfer/mocks.hpp"

typedef pfasst::vector_encap_traits<double, double>                     VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait>                  VectorEncapsulation;
typedef NiceMock<SweeperMock<pfasst::sweeper_traits<VectorEncapTrait>>> SweeperType;
typedef pfasst::transfer_traits<SweeperType, SweeperType, 1>            TransferTraits;
typedef TransferMock<TransferTraits>                                    TransferType;
typedef NiceMock<QuadratureMock<double>>                                QuadType;


typedef ::testing::Types<SDC<TransferType>> SDCTypes;
INSTANTIATE_TYPED_TEST_CASE_P(SDC, Concepts, SDCTypes);


class Interface
  : public ::testing::Test
{
  protected:
    shared_ptr<SDC<TransferType>> controller;

    shared_ptr<pfasst::Status<double>> status;

    virtual void SetUp()
    {
      this->controller = make_shared<SDC<TransferType>>();
      this->status = make_shared<pfasst::Status<double>>();
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
  controller->status() = status;
  controller->status()->time() = 42.0;
  EXPECT_THAT(controller->get_status()->get_time(), Eq(42.0));
}

TEST_F(Interface, computes_number_steps_fails_if_tend_or_dt_not_set)
{
  controller->status() = status;
  EXPECT_THROW(controller->get_num_steps(), logic_error);

  controller->status() = status;
  controller->status()->t_end() = 4.2;
  EXPECT_THROW(controller->get_num_steps(), logic_error);
}

TEST_F(Interface, computes_number_steps)
{
  controller->status() = status;
  controller->status()->t_end() = 4.2;
  controller->status()->dt() = 0.1;
  EXPECT_THAT(controller->get_num_steps(), Eq(42));
}


class Setup
  : public ::testing::Test
{
  protected:
    shared_ptr<SDC<TransferType>> controller;

    vector<double> nodes{0.0, 0.5, 1.0};
    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<SweeperType> sweeper;
    shared_ptr<TransferType> transfer;
    shared_ptr<QuadType> quad;

    virtual void SetUp()
    {
      this->controller = make_shared<SDC<TransferType>>();
      this->status = make_shared<pfasst::Status<double>>();
      this->controller->status() = status;

      this->quad = make_shared<QuadType>();
      ON_CALL(*(this->quad.get()), right_is_node()).WillByDefault(Return(true));
      ON_CALL(*(this->quad.get()), get_nodes()).WillByDefault(ReturnRef(this->nodes));
      ON_CALL(*(this->quad.get()), get_num_nodes()).WillByDefault(Return(3));

      this->sweeper = make_shared<SweeperType>();
      ON_CALL(*(this->sweeper.get()), get_quadrature()).WillByDefault(Return(this->quad));
      ON_CALL(*(this->sweeper.get()), status()).WillByDefault(ReturnRef(this->status));
      ON_CALL(*(this->sweeper.get()), get_status()).WillByDefault(Return(this->status));
    }
};

TEST_F(Setup, adding_coarser_level)
{
  ASSERT_THAT(controller->get_num_levels(), Eq(0));

  controller->add_sweeper(sweeper);
  EXPECT_THAT(controller->get_num_levels(), Eq(1));
}

TEST_F(Setup, a_level_must_be_added)
{
  controller->status()->t_end() = 0.2;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;

  EXPECT_THROW(controller->setup(), logic_error);

  controller->add_sweeper(sweeper);

  EXPECT_CALL(*(sweeper.get()), status()).Times(1);
  EXPECT_CALL(*(sweeper.get()), setup()).Times(1);
  controller->setup();
}

TEST_F(Setup, setup_required_for_running)
{
  controller->status()->t_end() = 0.1;
  controller->status()->dt() = 0.1;
  controller->status()->max_iterations() = 1;

  controller->add_sweeper(sweeper);

  ASSERT_FALSE(controller->is_ready());
  EXPECT_THROW(controller->run(), logic_error);

  EXPECT_CALL(*(sweeper.get()), status()).Times(1);
  EXPECT_CALL(*(sweeper.get()), setup()).Times(1);
  controller->setup();

  EXPECT_TRUE(controller->is_ready());
  controller->run();
}


class Logic
  : public ::testing::Test
{
  protected:
    shared_ptr<SDC<TransferType>> controller;

    vector<double> nodes{0.0, 0.5, 1.0};
    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<SweeperType> sweeper;
    shared_ptr<QuadType> quad;

    virtual void SetUp()
    {
      this->controller = make_shared<SDC<TransferType>>();
      this->status = make_shared<pfasst::Status<double>>();
      this->controller->status() = status;
      this->quad = make_shared<QuadType>();
      ON_CALL(*(this->quad.get()), right_is_node()).WillByDefault(Return(true));
      ON_CALL(*(this->quad.get()), get_nodes()).WillByDefault(ReturnRef(this->nodes));
      ON_CALL(*(this->quad.get()), get_num_nodes()).WillByDefault(Return(3));

      this->sweeper = make_shared<SweeperType>();
      ON_CALL(*(this->sweeper.get()), get_quadrature()).WillByDefault(Return(this->quad));
      ON_CALL(*(this->sweeper.get()), status()).WillByDefault(ReturnRef(this->status));
      ON_CALL(*(this->sweeper.get()), get_status()).WillByDefault(Return(this->status));

      this->controller->add_sweeper(this->sweeper);
    }
};

TEST_F(Logic, advance_in_time_with_sufficient_t_end)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 2.0;

  EXPECT_CALL(*(sweeper.get()), advance()).Times(1);

  EXPECT_TRUE(controller->advance_time());

  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.1));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(2));
}

TEST_F(Logic, advance_in_time_with_insufficient_t_end)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 1.0;

  EXPECT_CALL(*(sweeper.get()), advance()).Times(0);

  EXPECT_FALSE(controller->advance_time());

  EXPECT_THAT(controller->get_status()->get_time(), Eq(1.0));
  EXPECT_THAT(controller->get_status()->get_step(), Eq(1));
}

TEST_F(Logic, advance_in_time_multiple_steps_at_once)
{
  controller->status()->dt() = 0.1;
  controller->status()->time() = 1.0;
  controller->status()->step() = 1;
  controller->status()->t_end() = 2.0;

  EXPECT_CALL(*(sweeper.get()), advance()).Times(1);

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

  EXPECT_CALL(*(sweeper.get()), converged()).Times(1).WillRepeatedly(Return(false));

  EXPECT_FALSE(controller->advance_iteration());

  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(1));
}

TEST_F(Logic, advance_iteration)
{
  controller->status()->max_iterations() = 5;
  ASSERT_THAT(controller->get_status()->get_iteration(), Eq(0));
  ASSERT_THAT(controller->get_status()->get_max_iterations(), Eq(5));

  EXPECT_CALL(*(sweeper.get()), converged()).Times(1);
  EXPECT_CALL(*(sweeper.get()), save()).Times(1);

  EXPECT_TRUE(controller->advance_iteration());

  EXPECT_THAT(controller->get_status()->get_iteration(), Eq(1));
}

TEST_F(Logic, single_time_step_sdc)
{
  controller->status()->max_iterations() = 3;
  controller->status()->dt() = 0.1;
  controller->status()->time() = 0.0;
  controller->status()->t_end() = 0.1;

  EXPECT_CALL(*(sweeper.get()), status()).Times(1);
  EXPECT_CALL(*(sweeper.get()), setup()).Times(1);
  controller->setup();
  Mock::VerifyAndClearExpectations(&(*(sweeper.get())));

  EXPECT_CALL(*(sweeper.get()), converged()).Times(4);
  EXPECT_CALL(*(sweeper.get()), save()).Times(3);

  EXPECT_CALL(*(sweeper.get()), pre_predict()).Times(1);
  EXPECT_CALL(*(sweeper.get()), predict()).Times(1);
  EXPECT_CALL(*(sweeper.get()), post_predict()).Times(1);

  EXPECT_CALL(*(sweeper.get()), pre_sweep()).Times(3);
  EXPECT_CALL(*(sweeper.get()), sweep()).Times(3);
  EXPECT_CALL(*(sweeper.get()), post_sweep()).Times(3);

  EXPECT_CALL(*(sweeper.get()), advance()).Times(0);

  controller->run();
  EXPECT_THAT(controller->status()->get_step(), Eq(0));
  EXPECT_THAT(controller->status()->get_iteration(), Eq(3));
}


TEST_MAIN()
