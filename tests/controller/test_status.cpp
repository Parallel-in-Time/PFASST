#include "fixtures/test_helpers.hpp"

#include <pfasst/controller/status.hpp>
#include <pfasst/comm/communicator.hpp>

typedef pfasst::comm::Communicator CommType;


typedef ::testing::Types<pfasst::Status<double>> StatusTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Status, Concepts, StatusTypes);


class Interface
  : public ::testing::Test
{
  protected:
    pfasst::Status<double> status;
};

TEST_F(Interface, has_a_step)
{
  EXPECT_THAT(status.get_step(), Eq(0));

  status.step() = 1;
  EXPECT_THAT(status.get_step(), Eq(1));
}

TEST_F(Interface, has_an_iteration)
{
  EXPECT_THAT(status.get_iteration(), Eq(0));

  status.iteration() = 1;
  EXPECT_THAT(status.get_iteration(), Eq(1));
}

TEST_F(Interface, has_a_time_point)
{
  EXPECT_THAT(status.get_time(), Eq(0.0));

  status.time() = 1.42;
  EXPECT_THAT(status.get_time(), Eq(1.42));
}

TEST_F(Interface, has_a_time_delta)
{
  EXPECT_THAT(status.get_dt(), Eq(0.0));

  status.dt() = 0.42;
  EXPECT_THAT(status.get_dt(), Eq(0.42));
}

TEST_F(Interface, has_a_state)
{
  EXPECT_THAT(status.get_state(), Eq(pfasst::State::UNKNOWN));

  status.state() = pfasst::State::CONVERGED;
  EXPECT_THAT(status.get_state(), Eq(pfasst::State::CONVERGED));
}

TEST_F(Interface, has_a_residual)
{
  EXPECT_THAT(status.get_residual(), Eq(0.0));

  status.residual() = 0.1;
  EXPECT_THAT(status.get_residual(), Eq(0.1));
}


class Communication
  : public ::testing::Test
{
  protected:
    shared_ptr<pfasst::Status<double>> status;
    shared_ptr<CommType> comm;

  public:
    virtual void SetUp()
    {
      this->status = make_shared<pfasst::Status<double>>();
      this->comm = make_shared<CommType>();
    }
};

TEST_F(Communication, can_be_send) {
  status->send(comm, 1, 0, true);

  status->send(comm, 1, 0, false);
}

TEST_F(Communication, can_be_received) {
  status->recv(comm, 1, 0, true);

  status->recv(comm, 1, 0, false);
}

TEST_F(Communication, can_be_broadcasted) {
  status->bcast(comm, 0);
}


TEST_MAIN()
