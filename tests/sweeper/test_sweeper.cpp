#include "fixtures/test_helpers.hpp"

#include <vector>
using namespace std;

#include <pfasst/sweeper/interface.hpp>
using pfasst::Sweeper;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
typedef pfasst::vector_encap_traits<double, double> VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait> VectorEncapsulation;

#include "quadrature/mocks.hpp"
#include "controller/mocks.hpp"


typedef ::testing::Types<Sweeper<pfasst::sweeper_traits<VectorEncapTrait>>> SweeperTypes;
INSTANTIATE_TYPED_TEST_CASE_P(Sweeper, Concepts, SweeperTypes);


class Setup
  : public ::testing::Test
{
  protected:
    typedef          Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> sweeper_type;
    typedef typename sweeper_type::encap_type                          encap_type;

    sweeper_type sweeper;

    vector<double> nodes{0.0, 0.5, 1.0};
    shared_ptr<NiceMock<QuadratureMock<double>>> quadrature = make_shared<NiceMock<QuadratureMock<double>>>();

    virtual void SetUp()
    {
      ON_CALL(*(quadrature.get()), get_num_nodes()).WillByDefault(Return(nodes.size()));
      ON_CALL(*(quadrature.get()), get_nodes()).WillByDefault(ReturnRef(nodes));
    }
};

TEST_F(Setup, quadrature_is_required_for_setup)
{
  EXPECT_THAT(sweeper.quadrature(), IsNull());
  EXPECT_THAT(sweeper.get_quadrature(), IsNull());

  EXPECT_THROW(sweeper.setup(), runtime_error);

  EXPECT_CALL(*(quadrature.get()), get_num_nodes()).Times(1);
  EXPECT_CALL(*(quadrature.get()), get_nodes()).Times(1);

  sweeper.quadrature() = quadrature;
  EXPECT_THAT(sweeper.quadrature(), NotNull());
  EXPECT_THAT(sweeper.get_quadrature(), NotNull());

  sweeper.setup();
}

TEST_F(Setup, state_data_initialized_after_setup)
{
  EXPECT_THROW(sweeper.get_initial_state(), runtime_error);
  EXPECT_THROW(sweeper.initial_state(), runtime_error);
  EXPECT_THAT(sweeper.get_end_state(), IsNull());
  EXPECT_THAT(sweeper.get_states(), IsEmpty());
  EXPECT_THAT(sweeper.get_previous_states(), IsEmpty());
  EXPECT_THAT(sweeper.tau(), IsEmpty());
  EXPECT_THAT(sweeper.get_residuals(), IsEmpty());

  sweeper.quadrature() = quadrature;
  auto num_nodes = quadrature->get_num_nodes();
  sweeper.setup();

  EXPECT_THAT(sweeper.get_status(), IsNull());

  EXPECT_THAT(sweeper.get_initial_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state(), NotNull());

  EXPECT_THAT(sweeper.get_states(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_states(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_states(), Not(MutuallyEqual()));

  EXPECT_THAT(sweeper.get_previous_states(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_previous_states(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_previous_states(), Not(MutuallyEqual()));

  EXPECT_THAT(sweeper.get_tau(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_tau(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_tau(), Not(MutuallyEqual()));

  EXPECT_THAT(sweeper.get_residuals(), SizeIs(num_nodes + 1));
  EXPECT_THAT(sweeper.get_residuals(), Each(NotNull()));
  EXPECT_THAT(sweeper.get_residuals(), Not(MutuallyEqual()));
}


class DataAccess
  : public ::testing::Test
{
  protected:
    typedef          Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> sweeper_type;
    typedef typename sweeper_type::encap_type                          encap_type;

    sweeper_type sweeper;

    vector<double> nodes{0.0, 0.5, 1.0};
    shared_ptr<NiceMock<QuadratureMock<double>>> quadrature = make_shared<NiceMock<QuadratureMock<double>>>();

    shared_ptr<encap_type> encap = make_shared<encap_type>(vector<double>{1.0, 2.0, 3.0});
    shared_ptr<StatusMock<double>> status = make_shared<StatusMock<double>>();

    virtual void SetUp()
    {
      sweeper.encap_factory()->set_size(3);
      ON_CALL(*(quadrature.get()), get_num_nodes()).WillByDefault(Return(nodes.size()));
      ON_CALL(*(quadrature.get()), get_nodes()).WillByDefault(ReturnRef(nodes));
      sweeper.quadrature() = quadrature;
      sweeper.setup();
    }
};

TEST_F(DataAccess, initial_state_for_modification)
{
  sweeper.initial_state() = encap;
  EXPECT_THAT(sweeper.initial_state()->data(), Pointwise(Eq(), encap->data()));

  *(sweeper.get_initial_state().get()) = vector<double>{1.0, 1.0, 1.0};
  EXPECT_THAT(sweeper.initial_state()->data(), Each(1.0));
}

TEST_F(DataAccess, tau_for_modification)
{
  sweeper.tau() = vector<shared_ptr<encap_type>>{encap, encap, encap};
  EXPECT_THAT(sweeper.get_tau(), Each(Eq(encap)));
}

TEST_F(DataAccess, states_after_spreading_initial_state)
{
  sweeper.initial_state() = encap;
  sweeper.spread();
  EXPECT_THAT(sweeper.get_initial_state(), Eq(encap));
  EXPECT_THAT(sweeper.get_initial_state()->data(), Pointwise(Eq(), encap->data()));
  EXPECT_THAT(sweeper.get_states(), Not(MutuallyEqual()));
  for (auto state : sweeper.get_states()) {
    EXPECT_THAT(state->data(), Pointwise(Eq(), encap->data()));
  }
}

TEST_F(DataAccess, previous_states_after_spreading_initial_state_and_saving)
{
  sweeper.initial_state() = encap;
  sweeper.spread();
  sweeper.save();
  EXPECT_THAT(sweeper.get_initial_state(), Eq(encap));
  EXPECT_THAT(sweeper.get_initial_state()->data(), Pointwise(Eq(), encap->data()));
  EXPECT_THAT(sweeper.get_previous_states(), Not(MutuallyEqual()));
  for (auto state : sweeper.get_previous_states()) {
    EXPECT_THAT(state->data(), Pointwise(Eq(), encap->data()));
  }
}

TEST_F(DataAccess, status_for_modification)
{
  ASSERT_THAT(sweeper.get_status(), IsNull());

  sweeper.status() = status;
  EXPECT_THAT(sweeper.get_status(), Eq(status));
}


class Interface
  : public ::testing::Test
{
  protected:
    typedef          Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> sweeper_type;
    typedef typename sweeper_type::encap_type                          encap_type;

    sweeper_type sweeper;
    shared_ptr<NiceMock<StatusMock<double>>> status = make_shared<NiceMock<StatusMock<double>>>();

    virtual void SetUp()
    {
      sweeper.status() = status;
    }
};

TEST_F(Interface, no_implementation_of_advancing)
{
  EXPECT_THROW(sweeper.advance(), pfasst::NotImplementedYet);
}

TEST_F(Interface, no_implementation_of_reevaluation)
{
  EXPECT_THROW(sweeper.reevaluate(), pfasst::NotImplementedYet);
}

TEST_F(Interface, no_implementation_of_residual_computation)
{
  EXPECT_FALSE(sweeper.converged());

  sweeper.set_abs_residual_tol(1.0);
  EXPECT_THROW(sweeper.converged(), pfasst::NotImplementedYet);

  sweeper.set_abs_residual_tol(0.0);
  EXPECT_FALSE(sweeper.converged());

  sweeper.set_rel_residual_tol(1.0);
  EXPECT_THROW(sweeper.converged(), pfasst::NotImplementedYet);
}


class Logic
  : public ::testing::Test
{
  protected:
    typedef          Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> sweeper_type;
    typedef typename sweeper_type::encap_type                          encap_type;

    sweeper_type sweeper;

    vector<double> nodes{0.0, 0.5, 1.0};
    shared_ptr<NiceMock<QuadratureMock<double>>> quadrature = make_shared<NiceMock<QuadratureMock<double>>>();
    shared_ptr<NiceMock<StatusMock<double>>> status = make_shared<NiceMock<StatusMock<double>>>();

    shared_ptr<encap_type> encap = make_shared<encap_type>(vector<double>{1.0, 2.0, 3.0});

    virtual void SetUp()
    {
      sweeper.encap_factory()->set_size(3);
      ON_CALL(*(quadrature.get()), get_num_nodes()).WillByDefault(Return(nodes.size()));
      ON_CALL(*(quadrature.get()), get_nodes()).WillByDefault(ReturnRef(nodes));
      sweeper.quadrature() = quadrature;
      sweeper.setup();
      sweeper.status() = status;
    }
};

TEST_F(Logic, post_predict_finishes_end_state_if_right_is_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(true));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  sweeper.post_predict();
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Eq(), encap->data()));
}

TEST_F(Logic, post_predict_fails_if_right_is_not_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(false));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  EXPECT_THROW(sweeper.post_predict(), pfasst::NotImplementedYet);
}

TEST_F(Logic, post_sweep_finishes_end_state_if_right_is_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(true));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  sweeper.post_sweep();
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Eq(), encap->data()));
}

TEST_F(Logic, post_sweep_fails_if_right_is_not_node)
{
  ON_CALL(*(quadrature.get()), right_is_node()).WillByDefault(Return(false));

  sweeper.initial_state() = encap;
  sweeper.spread();

  EXPECT_THAT(sweeper.get_end_state(), NotNull());
  EXPECT_THAT(sweeper.get_end_state()->data(), Pointwise(Not(Eq()), encap->data()));

  EXPECT_THROW(sweeper.post_sweep(), pfasst::NotImplementedYet);
}


TEST_MAIN()
