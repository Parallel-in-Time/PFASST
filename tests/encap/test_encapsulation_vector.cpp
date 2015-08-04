#include "fixtures/test_helpers.hpp"

#include <algorithm>
#include <memory>
#include <type_traits>
#include <vector>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <Eigen/Dense>
#include <leathers/pop>
template<typename precision>
using Matrix = Eigen::Matrix<precision, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

#include <pfasst/encap/traits.hpp>
#include <pfasst/encap/vector.hpp>
typedef pfasst::vector_encap_traits<double, double> VectorEncapTrait;
typedef pfasst::encap::Encapsulation<VectorEncapTrait> VectorEncapsulation;

#include "comm/mocks.hpp"


typedef ::testing::Types<VectorEncapsulation> EncapTypes;
INSTANTIATE_TYPED_TEST_CASE_P(VectorEncap, Concepts, EncapTypes);


TEST(Construction, empty_constructible)
{
  VectorEncapsulation vec;
  EXPECT_THAT(vec.get_data().size(), Eq(0));
  EXPECT_THAT(vec.get_data(), Eq(vector<double>{}));
}

TEST(Construction, data_constructible)
{
  VectorEncapsulation vec(vector<double>{1.0, 2.0, 3.0});

  EXPECT_THAT(vec.get_data().size(), Eq(3));
  EXPECT_THAT(vec.get_data(), Eq(vector<double>{1.0, 2.0, 3.0}));
}

TEST(DataAccession, assignable)
{
  VectorEncapsulation vec;
  vector<double> data{1.0, 2.0, 3.0};

  vec = data;
  EXPECT_THAT(vec.get_data().size(), Eq(data.size()));
  EXPECT_THAT(vec.get_data(), Eq(data));
}


TEST(Operation, zeroing_out)
{
  vector<double> vec{1.0, 2.0, 3.0};
  VectorEncapsulation x(vec);
  EXPECT_THAT(x.get_data(), Pointwise(Eq(), vec));

  x.zero();
  EXPECT_THAT(x.get_data(), Each(Eq(0.0)));
}

TEST(Operation, in_place_axpy)
{
  VectorEncapsulation vec_x(vector<double>{1.0, 2.0, 3.0});
  shared_ptr<VectorEncapsulation> vec_y = \
    make_shared<VectorEncapsulation>(vector<double>{1.0, 1.0, 1.0});

  vec_x.scaled_add(0.5, vec_y);
  EXPECT_THAT(vec_x.get_data(), Pointwise(Eq(), vector<double>{1.5, 2.5, 3.5}));
}

TEST(Operation, global_axpy)
{
  shared_ptr<VectorEncapsulation> vec_x = \
    make_shared<VectorEncapsulation>(vector<double>{1.0, 2.0, 3.0});
  shared_ptr<VectorEncapsulation> vec_y = \
    make_shared<VectorEncapsulation>(vector<double>{1.0, 1.0, 1.0});

  auto result = pfasst::encap::axpy(0.5, vec_x, vec_y);
  EXPECT_THAT(result->get_data(), Pointwise(Eq(), vector<double>{1.5, 2.0, 2.5}));
}

TEST(Operation, norm0_as_member)
{
  VectorEncapsulation vec_x(vector<double>{1.0, -4.0, 3.0});
  EXPECT_THAT(vec_x.norm0(), Eq(4.0));
}

TEST(Operation, global_norm0)
{
  shared_ptr<VectorEncapsulation> vec_x = \
    make_shared<VectorEncapsulation>(vector<double>{1.0, -4.0, 3.0});
  EXPECT_THAT(norm0(vec_x), Eq(4.0));
}


TEST(MatrixApplication, identity)
{
  vector<double> data{1.0, 2.0, 3.0};
  vector<shared_ptr<VectorEncapsulation>> vec(3);
  generate(vec.begin(), vec.end(),
           [data]() { return make_shared<VectorEncapsulation>(data); });
  for_each(vec.cbegin(), vec.cend(), [data](const shared_ptr<VectorEncapsulation>& xi) {
    EXPECT_THAT(xi->get_data(), Pointwise(Eq(), data));
  });
  Matrix<double> mat = Matrix<double>::Identity(3, 3);

  auto result_mat_mul_vec = pfasst::encap::mat_mul_vec(1.0, mat, vec);
  for_each(result_mat_mul_vec.cbegin(), result_mat_mul_vec.cend(),
           [data](const shared_ptr<VectorEncapsulation>& xi) {
             EXPECT_THAT(xi->get_data(), Pointwise(Eq(), data));
           });

  vector<shared_ptr<VectorEncapsulation>> result_mat_apply(result_mat_mul_vec);
  pfasst::encap::mat_apply(result_mat_apply, 1.0, mat, vec, true);
  for_each(result_mat_apply.cbegin(), result_mat_apply.cend(),
           [data](const shared_ptr<VectorEncapsulation>& xi) {
             EXPECT_THAT(xi->get_data(), Pointwise(Eq(), data));
           });
}

TEST(MatrixApplication, zero_matrix)
{
  vector<double> data{1.0, 2.0, 3.0};
  vector<shared_ptr<VectorEncapsulation>> vec(3);
  generate(vec.begin(), vec.end(),
           [data]() { return make_shared<VectorEncapsulation>(data); });
  for_each(vec.cbegin(), vec.cend(), [data](const shared_ptr<VectorEncapsulation>& xi) {
    EXPECT_THAT(xi->get_data(), Pointwise(Eq(), data));
  });
  Matrix<double> mat = Matrix<double>::Zero(3, 3);

  auto result_mat_mul_vec = pfasst::encap::mat_mul_vec(1.0, mat, vec);
  for_each(result_mat_mul_vec.cbegin(), result_mat_mul_vec.cend(),
           [](const shared_ptr<VectorEncapsulation>& xi) {
             EXPECT_THAT(xi->get_data(), Pointwise(Eq(), vector<double>{0.0, 0.0, 0.0}));
           });

  vector<shared_ptr<VectorEncapsulation>> result_mat_apply(result_mat_mul_vec);
  pfasst::encap::mat_apply(result_mat_apply, 1.0, mat, vec, true);
  for_each(result_mat_apply.cbegin(), result_mat_apply.cend(),
           [data](const shared_ptr<VectorEncapsulation>& xi) {
             EXPECT_THAT(xi->get_data(), Pointwise(Eq(), vector<double>{0.0, 0.0, 0.0}));
           });
}

TEST(MatrixApplication, all_ones)
{
  vector<double> data{1.0, 2.0, 3.0};
  vector<shared_ptr<VectorEncapsulation>> vec(3);
  generate(vec.begin(), vec.end(),
           [data]() { return make_shared<VectorEncapsulation>(data); });
  for_each(vec.cbegin(), vec.cend(), [data](const shared_ptr<VectorEncapsulation>& xi) {
    EXPECT_THAT(xi->get_data(), Pointwise(Eq(), data));
  });
  Matrix<double> mat = Matrix<double>::Ones(3, 3);

  auto result_mat_mul_vec = pfasst::encap::mat_mul_vec(1.0, mat, vec);
  for_each(result_mat_mul_vec.cbegin(), result_mat_mul_vec.cend(),
           [](const shared_ptr<VectorEncapsulation>& xi) {
             EXPECT_THAT(xi->get_data(), Pointwise(Eq(), vector<double>{3.0, 6.0, 9.0}));
           });

  vector<shared_ptr<VectorEncapsulation>> result_mat_apply(result_mat_mul_vec);
  pfasst::encap::mat_apply(result_mat_apply, 1.0, mat, vec, true);
  for_each(result_mat_apply.cbegin(), result_mat_apply.cend(),
           [data](const shared_ptr<VectorEncapsulation>& xi) {
             EXPECT_THAT(xi->get_data(), Pointwise(Eq(), vector<double>{3.0, 6.0, 9.0}));
           });
}


TEST(Communication, sending)
{
  shared_ptr<NiceMock<CommMock>> comm = make_shared<NiceMock<CommMock>>();
  shared_ptr<VectorEncapsulation> vec = \
    make_shared<VectorEncapsulation>(vector<double>{1.0, 2.0, 3.0});
  EXPECT_CALL(*(comm.get()), send(Matcher<const double*>(vec->get_data().data()), 3, 1, 0)).Times(1);
  vec->send(comm, 1, 0, true);

  EXPECT_CALL(*(comm.get()), isend(Matcher<const double*>(vec->get_data().data()), 3, 1, 0)).Times(1);
  vec->send(comm, 1, 0, false);
}

TEST(Communication, receiving)
{
  shared_ptr<NiceMock<CommMock>> comm = make_shared<NiceMock<CommMock>>();
  shared_ptr<VectorEncapsulation> vec = \
    make_shared<VectorEncapsulation>(vector<double>{1.0, 2.0, 3.0});
  EXPECT_CALL(*(comm.get()), recv(Matcher<double*>(vec->data().data()), 3, 1, 0)).Times(1);
  vec->recv(comm, 1, 0, true);

  EXPECT_CALL(*(comm.get()), irecv(Matcher<double*>(vec->data().data()), 3, 1, 0)).Times(1);
  vec->recv(comm, 1, 0, false);
}

TEST(Communication, broadcasting)
{
  shared_ptr<NiceMock<CommMock>> comm = make_shared<NiceMock<CommMock>>();
  shared_ptr<VectorEncapsulation> vec = \
    make_shared<VectorEncapsulation>(vector<double>{1.0, 2.0, 3.0});
  EXPECT_CALL(*(comm.get()), bcast(Matcher<double*>(vec->data().data()), 3, 0)).Times(1);
  vec->bcast(comm, 0);
}


TEST(Factory, predefine_size)
{
  pfasst::encap::EncapsulationFactory<VectorEncapTrait> null_factory;
  EXPECT_THAT(null_factory.size(), Eq(0));

  pfasst::encap::EncapsulationFactory<VectorEncapTrait> sized_factory(3);
  EXPECT_THAT(sized_factory.size(), Eq(3));

  null_factory.set_size(3);
  EXPECT_THAT(null_factory.size(), Eq(3));
}

TEST(Factory, create_vector_encap)
{
  pfasst::encap::EncapsulationFactory<VectorEncapTrait> factory(3);
  auto encap = factory.create();
  EXPECT_THAT(encap->get_data().size(), Eq(3));

  factory.set_size(5);
  EXPECT_THAT(factory.size(), Eq(5));
  auto encap_5 = factory.create();
  EXPECT_THAT(encap_5->get_data().size(), Eq(5));
}

TEST_MAIN()
