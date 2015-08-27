#include "fixtures/test_helpers.hpp"

#include <pfasst/quadrature/quadrature.hpp>


template<typename precision>
class QuadratureMock
  : public pfasst::quadrature::IQuadrature<precision>
{
  public:
    MOCK_CONST_METHOD0_T(get_nodes, vector<precision>&());
    MOCK_CONST_METHOD0_T(get_num_nodes, size_t());
    MOCK_CONST_METHOD0_T(left_is_node, bool());
    MOCK_CONST_METHOD0_T(right_is_node, bool());
    MOCK_CONST_METHOD0_T(expected_error, precision());

    MOCK_CONST_METHOD0_T(get_b_mat, Matrix<double>&());
};
