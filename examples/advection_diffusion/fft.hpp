/*
 * FFT helper class.
 *
 * Please note: side affects galore!  This is not my best work...
 */

#ifndef _FFT_HPP_
#define _FFT_HPP_

#include <map>

#include <pfasst/encap/vector.hpp>

#include <fftw3.h>

template<typename time>
class FFT {

  using dvector = pfasst::encap::VectorEncapsulation<double,time>;

  struct workspace {
    fftw_plan        ffft;
    fftw_plan        ifft;
    fftw_complex*    wk;
    complex<double>* z;
  };

  map<int,workspace*> workspaces;

public:

  ~FFT()
  {
    for (auto& x: workspaces) {
      auto* wk = std::get<1>(x);
      fftw_free(wk->wk);
      fftw_destroy_plan(wk->ffft);
      fftw_destroy_plan(wk->ifft);
      delete wk;
    }
    workspaces.clear();
  }

  workspace* get_workspace(int ndofs)
  {
    if (workspaces.find(ndofs) == workspaces.end()) {
      workspace* wk = new workspace;
      wk->wk = fftw_alloc_complex(ndofs);
      wk->ffft = fftw_plan_dft_1d(ndofs, wk->wk, wk->wk, FFTW_FORWARD, FFTW_ESTIMATE);
      wk->ifft = fftw_plan_dft_1d(ndofs, wk->wk, wk->wk, FFTW_BACKWARD, FFTW_ESTIMATE);
      wk->z = reinterpret_cast<complex<double>*>(wk->wk);
      workspaces.insert(pair<int,workspace*>(ndofs, wk));
    }

    return workspaces[ndofs];
  }

  complex<double>* forward(const dvector& x)
  {
    workspace *wk = get_workspace(x.size());
    for (size_t i=0; i<x.size(); i++)
      wk->z[i] = x[i];
    fftw_execute_dft(wk->ffft, wk->wk, wk->wk);
    return wk->z;
  }

  void backward(dvector &x)
  {
    workspace *wk = get_workspace(x.size());
    fftw_execute_dft(wk->ifft, wk->wk, wk->wk);
    for (size_t i=0; i<x.size(); i++)
      x[i] = real(wk->z[i]);
  }

};

#endif
