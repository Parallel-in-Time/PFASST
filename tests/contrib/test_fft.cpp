#include "fixtures/test_helpers.hpp"

#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <vector>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <boost/math/constants/constants.hpp>
#include <leathers/pop>
using namespace boost::math::constants;


#include <pfasst/contrib/fft.hpp>
using pfasst::contrib::FFT;

#include <pfasst/encap/vector.hpp>
typedef pfasst::encap::VectorEncapsulation<double, double> VectorType;
#include <pfasst/logging.hpp>


typedef ::testing::Types<FFT<double>> FFTTypes;
INSTANTIATE_TYPED_TEST_CASE_P(FFT, Concepts, FFTTypes);

class Interface
  : public ::testing::Test
{
  protected:
    typedef FFT<double> fft_type;

    fft_type fft;
};

TEST_F(Interface, query_z_pointer_for_specific_num_dofs)
{
  complex<double>* z_ptr = fft.get_workspace(1)->z;
}


class DiscreteFastFourierTransform
  : public ::testing::TestWithParam<shared_ptr<VectorType>>
{
  protected:
    typedef FFT<double> fft_type;

    fft_type fft;
    shared_ptr<VectorType> values;

    virtual void SetUp()
    {
      this->values = GetParam();
    }

    vector<complex<double>> eval_base_function(const shared_ptr<VectorType> vec, const size_t k)
    {
      const size_t ndofs = vec->get_data().size();
      vector<complex<double>> result(ndofs);

      transform(vec->get_data().cbegin(), vec->get_data().cend(),
                result.begin(),
                [ndofs, k](const double& t) {
                  return exp(complex<double>(0.0, 1.0) * (two_pi<double>() / double(ndofs) * k * t));
                });

      return result;
    }

    vector<double> two_pi_k_t(const shared_ptr<VectorType> vec, const size_t& k)
    {
      const size_t ndofs = vec->get_data().size();
      vector<double> result(ndofs);

      transform(vec->get_data().cbegin(), vec->get_data().cend(),
                result.begin(),
                [k](const double& t) {
                  return cos(two_pi<double>() * k * t);
                });

      return result;
    }
};

TEST_P(DiscreteFastFourierTransform, forward_transform)
{
  size_t ndofs = values->get_data().size();
//   LOG(INFO) << "values: " << values->get_data();
  for (size_t k = 0; k < ndofs; ++k) {
    const double precision = k * ndofs * numeric_limits<double>::epsilon();

//     LOG(INFO) << "k=" << k;
//     LOG(INFO) << "  cos(2pi*k*t): " << two_pi_k_t(values, k);
    vector<complex<double>> forward(ndofs);
    auto* fft_forward = fft.forward(make_shared<VectorType>(two_pi_k_t(values, k)));
    for (size_t i = 0; i < ndofs; ++i) {
      forward[i] = fft_forward[i];
    }
//     LOG(INFO) << "       --> fft: " << forward;

    for (size_t i = 0; i < ndofs; ++i) {
      EXPECT_THAT(forward[i].imag(), DoubleNear(0.0, precision));

      if (i != k && i != (ndofs - k)) {
        EXPECT_THAT(forward[i].real(), DoubleNear(0.0, precision));
      } else if (ndofs % 2 == 0 && (i == 0 || i == ndofs / 2)) {
          EXPECT_THAT(forward[i].real(), DoubleNear(ndofs, precision));
      } else {
        // TODO: test the other cases where entries are of ndofs/2
      }
    }
  }
}


TEST_P(DiscreteFastFourierTransform, backward_transform)
{
  size_t ndofs = values->get_data().size();
//   LOG(INFO) << "values: " << values->get_data();
  for (size_t k = 0; k < ndofs; ++k) {
    const double precision = k * ndofs * numeric_limits<double>::epsilon();

//     LOG(INFO) << "k=" << k;
    shared_ptr<VectorType> test_values = make_shared<VectorType>(two_pi_k_t(values, k));
//     LOG(INFO) << "  cos(2pi*k*t): " << test_values->get_data();
    fft.forward(test_values);

    shared_ptr<VectorType> backward = make_shared<VectorType>();
    backward->data().resize(ndofs);
    backward->zero();
    shared_ptr<VectorType> expected = make_shared<VectorType>(backward->get_data());
    expected->scaled_add(ndofs, test_values);

    fft.backward(backward);

//     LOG(INFO) << "  --> backward: " << backward->get_data();
    for (size_t i = 0; i < ndofs; ++i) {
      EXPECT_THAT(backward->get_data()[i], DoubleNear(expected->get_data()[i], precision));
    }
  }
}

auto values_3 = make_shared<VectorType>(vector<double>{0.0, third<double>(), two_thirds<double>()});
auto values_4 = make_shared<VectorType>(vector<double>{0.0, 0.25, 0.5, 0.75});
auto values_5 = make_shared<VectorType>(vector<double>{0.0, 0.2, 0.4, 0.6, 0.8});
auto values_10 = make_shared<VectorType>(vector<double>{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9});

INSTANTIATE_TEST_CASE_P(DFT,
                        DiscreteFastFourierTransform,
                        ::testing::Values(values_4, values_5, values_10));


TEST_MAIN()
