#ifndef _PFASST__CONTRIB__FFTW_HPP_
#define _PFASST__CONTRIB__FFTW_HPP_

#include <complex>
using std::real;
#include <map>
#include <memory>
#include <utility>
using namespace std;

#include <leathers/push>
#include <leathers/all>
#include <fftw3.h>
#include <leathers/pop>

#include "pfasst/encap/vector.hpp"


namespace pfasst
{
  namespace contrib
  {
    //! TODO: rewrite to get rid of side-effects and use real RAII
    template<
      class precision
    >
    class FFT
    {
      protected:
        struct workspace {
          fftw_plan           ffft;
          fftw_plan           ifft;
          fftw_complex*       wk;
          complex<precision>* z;
        };

        map<size_t, shared_ptr<workspace>> workspaces;

      public:
        virtual ~FFT();

        shared_ptr<workspace> get_workspace(size_t ndofs);

        template<
          class time_precision
        >
        complex<precision>* forward(const shared_ptr<encap::VectorEncapsulation<time_precision, precision>> x);
        void backward(shared_ptr<encap::VectorEncapsulation<time_precision, precision>> x);
    };
  }  // ::pfasst::contrib
}  // ::pfasst

#include "pfasst/contrib/fft_impl.hpp"

#endif  // _PFASST__CONTRIB__FFTW_HPP_
