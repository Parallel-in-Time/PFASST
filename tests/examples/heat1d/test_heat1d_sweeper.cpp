#include "fixtures/test_helpers.hpp"

#include <vector>
using namespace std;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
typedef pfasst::vector_encap_traits<double, double> VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait> VectorEncapsulation;

#include "examples/heat1d/heat1d_sweeper.hpp"
template<class... Ts>
using Sweeper = pfasst::examples::heat1d::Heat1D<Ts...>;


class ProblemSetup
  : public ::testing::Test
{
  protected:
    typedef          Sweeper<pfasst::sweeper_traits<VectorEncapTrait>> sweeper_type;
    typedef typename sweeper_type::encap_type                          encap_type;

    shared_ptr<sweeper_type> sweeper;

    shared_ptr<pfasst::Status<double>> status = make_shared<pfasst::Status<double>>();

    vector<double> exact_t0{   0.000000000000000000000000
                            ,  0.7071067811865474617150085
                            ,  1.000000000000000000000000
                            ,  0.7071067811865475727373109
                            ,  1.224646799147353207173764e-16
                            , -0.7071067811865474617150085
                            , -1.000000000000000000000000
                            , -0.7071067811865476837596134};
    vector<double> exact_t001{   0.000000000000000000000000e+00
                              ,  0.7015456730923740336081096
                              ,  0.9921354055113971170953846
                              ,  0.7015456730923741446304120
                              ,  1.215015448680293835571449e-16
                              , -0.7015456730923740336081096
                              , -0.9921354055113971170953846
                              , -0.7015456730923742556527145};

    virtual void SetUp()
    {
      sweeper = make_shared<sweeper_type>(8);
      sweeper->status() = status;
    }
};

TEST_F(ProblemSetup, computes_exact_solution_at_t0)
{
  auto exact = sweeper->exact(0.0);

  EXPECT_THAT(exact->get_data(), Pointwise(DoubleNear(), exact_t0));
}

TEST_F(ProblemSetup, computes_exact_solution_at_t01)
{
  auto exact = sweeper->exact(0.01);

  EXPECT_THAT(exact->get_data(), Pointwise(DoubleNear(), exact_t001));
}



TEST_MAIN()
