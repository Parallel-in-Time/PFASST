#include "fixtures/test_helpers.hpp"

#include <pfasst/contrib/fft.hpp>
using pfasst::contrib::FFT;

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


class HowShouldINameThis
  : public ::testing::Test
{
  protected:
    typedef FFT<double> fft_type;

    fft_type fft;
};

TEST_F(HowShouldINameThis, query_non_existend_workspace)
{
  
}


TEST_MAIN()
