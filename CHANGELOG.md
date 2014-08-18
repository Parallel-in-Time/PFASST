# Changelog                                                                        {#page_changelog}

## v0.2.0 -- MPI PFASST (2014/08/XX)

XXX: DOI: [10.5281/zenodo.11047](http://dx.doi.org/10.5281/zenodo.11047)

### Notable Features

* Addition of MPI based PFASST.

### Details

* Addition of MPI based PFASST algorithm using the standard predictor stage in a block mode with 
  fixed iterations.
  ([#46][], [#57][], [#59][])

* Addition of a simple scalar example and appropriate tests.
  ([#61][], [#63][], [#76][])

* Further tests on proper calculation of quadrature nodes and weights
  ([#74][])

* Better handling of 3rd-party dependencies
  ([#53][], [#55][])

* Various tidying
  ([#56][], [#60][], [#77][])

* Basic Profiling support and test coverage report
  ([#78][])

[#46]: https://github.com/Parallel-in-Time/PFASST/pull/46
[#57]: https://github.com/Parallel-in-Time/PFASST/pull/56
[#59]: https://github.com/Parallel-in-Time/PFASST/pull/59
[#53]: https://github.com/Parallel-in-Time/PFASST/pull/53
[#55]: https://github.com/Parallel-in-Time/PFASST/pull/55
[#56]: https://github.com/Parallel-in-Time/PFASST/pull/56
[#60]: https://github.com/Parallel-in-Time/PFASST/pull/60
[#61]: https://github.com/Parallel-in-Time/PFASST/pull/61
[#63]: https://github.com/Parallel-in-Time/PFASST/pull/63
[#74]: https://github.com/Parallel-in-Time/PFASST/pull/74
[#76]: https://github.com/Parallel-in-Time/PFASST/pull/76
[#77]: https://github.com/Parallel-in-Time/PFASST/pull/77
[#78]: https://github.com/Parallel-in-Time/PFASST/pull/78

### Contributors

* Matthew Emmett, Lawrence Berkeley National Laboratory ([memmett][])
* Torbjörn Klatt, Jülich Supercomputing Centre ([torbjoernk][])
* Daniel Ruprecht, Institute of Computational Science, University of Lugano ([danielru][])
* Robert Speck, Jülich Supercomputing Centre ([pancetta][])

[memmett]: https://github.com/memmett
[torbjoernk]: https://github.com/torbjoernk
[danielru]: https://github.com/danielru
[pancetta]: https://github.com/pancetta


## v0.1.0 -- First Release (2014/07/25)

DOI: [10.5281/zenodo.11047](http://dx.doi.org/10.5281/zenodo.11047)

### Notable Features

* Initial release with basic implementations of SDC and MLSDC

### Details

* agreed on code style guideline using Astyle tool
  ([#2][], [#15][], [#21][], [#22][])
* few test cases for quadrature and advection-diffusion example
  ([#23][], [#25][], [#33][], [#37][])

[#2]: https://github.com/Parallel-in-Time/PFASST/pull/2
[#15]: https://github.com/Parallel-in-Time/PFASST/pull/15
[#21]: https://github.com/Parallel-in-Time/PFASST/pull/21
[#22]: https://github.com/Parallel-in-Time/PFASST/pull/22
[#23]: https://github.com/Parallel-in-Time/PFASST/pull/23
[#25]: https://github.com/Parallel-in-Time/PFASST/pull/25
[#33]: https://github.com/Parallel-in-Time/PFASST/pull/33
[#37]: https://github.com/Parallel-in-Time/PFASST/pull/37

### Contributors

* Matthew Emmett, Lawrence Berkeley National Laboratory ([memmett][])
* Torbjörn Klatt, Jülich Supercomputing Centre ([torbjoernk][])

[memmett]: https://github.com/memmett
[torbjoernk]: https://github.com/torbjoernk
