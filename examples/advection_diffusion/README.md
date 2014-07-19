# Advection / Diffusion                                         {#page_examples_advection_diffusion}

This directory contains several implementations of an advection/diffusion solver using the PFASST 
framework.

All of the solvers use the SDC sweeper defined in `advection_diffusion_sweeper.hpp`, and the FFT 
routines in `fft.hpp`.

The implementations are, in order of complexity:

* `vanilla_sdc.cpp` - basic example that uses an encapsulated IMEX sweeper.

* `serial_mlsdc.cpp` - basic multi-level version that uses polynomial interpolation in time and 
  spectral interpolation in space, as defined in `specrtal_transfer_1d.hpp`.

* `serial_mlsdc_autobuild.cpp` - same as above, but uses the "auto build" feature to shorten `main`.
