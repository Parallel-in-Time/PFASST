# Changelog                                                                        {#page_changelog}

## v0.4.0 -- The Big Improvement (2015/04/17)

DOI: [10.6084/m9.figshare.1381721](http://dx.doi.org/10.6084/m9.figshare.1381721)

    $ git diff --stat v0.3.0..v0.4.0
    > 131 files changed, 15182 insertions(+), 6940 deletions(-)

### Notable Features

* Iteration controll for PFASST via relative and absolute residual tolerances.
  ([#139][], [#141][], [#142][])

* Adding an implicit sweeper with support for LU decomposition of quadrature matrix.
  ([#150][])

* General improvements to the logging and configuration functionality.
  ([#144][], [#145][], [#146][], [#148][], [#163][], [#168][], [#169][], [#177][])

* Successful tests on supercomputers (i.e., JUQUEEN)
  ([#168][], [#176][])

* Improved code organization.
  ([#156][], [#157][], [#158][], [#170][], [#181][])

* Rewamped documentation.
  ([#149][], [#160][], [#179][], [#183][])

* Examples: Boris supports multiple particles and multi-level coarsening
  ([#152][], [#172][], [#180][])

### Details

* logging framework [Easylogging][] has been updated from v9.75 to v9.80.
  ([#163][])

* root finding algorithm for polynomials uses fixed number of iterations or residual tolerance.
  ([#162][])

* reworked how `pfasst::init(...)` is used.
  ([#148][])

* minimum _CMake_ version bumped to 2.8.*6*
  ([#153][])

* `git describe`-based version number as `pfasst::VERSION`
  ([#185][])

* script to generate and compile test coverage report was rewritten in Python 3.x
  ([#154][])

[#139]: https://github.com/Parallel-in-Time/PFASST/pull/139
[#141]: https://github.com/Parallel-in-Time/PFASST/pull/141
[#142]: https://github.com/Parallel-in-Time/PFASST/pull/142
[#144]: https://github.com/Parallel-in-Time/PFASST/pull/144
[#145]: https://github.com/Parallel-in-Time/PFASST/pull/145
[#146]: https://github.com/Parallel-in-Time/PFASST/pull/146
[#148]: https://github.com/Parallel-in-Time/PFASST/pull/148
[#149]: https://github.com/Parallel-in-Time/PFASST/pull/149
[#150]: https://github.com/Parallel-in-Time/PFASST/pull/150
[#152]: https://github.com/Parallel-in-Time/PFASST/pull/152
[#153]: https://github.com/Parallel-in-Time/PFASST/pull/153
[#154]: https://github.com/Parallel-in-Time/PFASST/pull/154
[#156]: https://github.com/Parallel-in-Time/PFASST/pull/156
[#157]: https://github.com/Parallel-in-Time/PFASST/pull/157
[#158]: https://github.com/Parallel-in-Time/PFASST/pull/158
[#160]: https://github.com/Parallel-in-Time/PFASST/pull/160
[#162]: https://github.com/Parallel-in-Time/PFASST/pull/162
[#163]: https://github.com/Parallel-in-Time/PFASST/pull/163
[#168]: https://github.com/Parallel-in-Time/PFASST/pull/168
[#169]: https://github.com/Parallel-in-Time/PFASST/pull/169
[#170]: https://github.com/Parallel-in-Time/PFASST/pull/170
[#172]: https://github.com/Parallel-in-Time/PFASST/pull/172
[#176]: https://github.com/Parallel-in-Time/PFASST/pull/176
[#177]: https://github.com/Parallel-in-Time/PFASST/pull/177
[#180]: https://github.com/Parallel-in-Time/PFASST/pull/180
[#181]: https://github.com/Parallel-in-Time/PFASST/pull/181
[#183]: https://github.com/Parallel-in-Time/PFASST/pull/183
[#185]: https://github.com/Parallel-in-Time/PFASST/pull/185

### Contributors

* Matthew Emmett, Lawrence Berkeley National Laboratory ([memmett][])
* Torbjörn Klatt, Jülich Supercomputing Centre ([torbjoernk][])
* Daniel Ruprecht, Institute of Computational Science, University of Lugano ([danielru][])
* Robert Speck, Jülich Supercomputing Centre ([pancetta][])
* Selman Terzi, Jülich Supercomputing Centre ([selmanTerzi][])

[memmett]: https://github.com/memmett
[torbjoernk]: https://github.com/torbjoernk
[danielru]: https://github.com/danielru
[pancetta]: https://github.com/pancetta
[selmanTerzi]: https://github.com/selmanTerzi

---

## v0.3.0 -- The Big Cleanup (2014/12/12)

DOI: [10.5281/zenodo.13221](http://dx.doi.org/10.5281/zenodo.13221)

    $ git diff --stat v0.2.0..v0.3.0
    > 84 files changed, 7687 insertions(+), 1968 deletions(-)

### Notable Features

* Complete rewrite of quadrature functions.
  ([#103][])

* Addition of iteration control mechanism for SDC and MLSDC. _(extention for PFASST is scheduled)_
  ([#108][])

* Introduced framework for passing command line parameters to programs using PFASST++.
  ([#88][])

* Versatile and colourful logging framework based on [Easylogging++][].
  ([#105][], [#123][])

* New Example: Boris-SDC.
  ([#113][])

### Details

* Using [Eigen3](http://eigen.tuxfamily.org/) matrix data types and functionality instead of
  [Boost's uBLAS](http://www.boost.org/doc/libs/1_57_0/libs/numeric/ublas/doc/index.html).
  ([#82][], [#83][], [#85][])

* Introduced callback interception points in sweepers: `post_step`, `post_sweep` and `post_predict`.
  ([#107][])

* In sweepers: `state(0)` is not assumed to be same as `start_state`.
  ([#132][])

* Examples are now in their own namespace `pfasst::examples::`
  ([#138][])

* Bunch of CMake-related fixes and extentions.
  ([#83][], [#96][], [#101][], [#109][], [#112][], [#119][], [#121][], [#131][], [#133][], [#135][])

* Example of using PFASST++ with a Makefile project (see advection-diffusion example).
  ([#129][])

* Support for [HashDist](https://github.com/hashdist/hashdist) to provide dependencies.
  ([#129][])

* A few corrections and extentions to the documentation.
  ([#117][], [#120][])

[#82]: https://github.com/Parallel-in-Time/PFASST/pull/82
[#83]: https://github.com/Parallel-in-Time/PFASST/pull/83
[#85]: https://github.com/Parallel-in-Time/PFASST/pull/85
[#88]: https://github.com/Parallel-in-Time/PFASST/pull/88
[#96]: https://github.com/Parallel-in-Time/PFASST/pull/96
[#101]: https://github.com/Parallel-in-Time/PFASST/pull/101
[#103]: https://github.com/Parallel-in-Time/PFASST/pull/103
[#105]: https://github.com/Parallel-in-Time/PFASST/pull/105
[#107]: https://github.com/Parallel-in-Time/PFASST/pull/107
[#108]: https://github.com/Parallel-in-Time/PFASST/pull/108
[#109]: https://github.com/Parallel-in-Time/PFASST/pull/109
[#112]: https://github.com/Parallel-in-Time/PFASST/pull/112
[#113]: https://github.com/Parallel-in-Time/PFASST/pull/113
[#117]: https://github.com/Parallel-in-Time/PFASST/pull/117
[#119]: https://github.com/Parallel-in-Time/PFASST/pull/119
[#120]: https://github.com/Parallel-in-Time/PFASST/pull/120
[#121]: https://github.com/Parallel-in-Time/PFASST/pull/121
[#123]: https://github.com/Parallel-in-Time/PFASST/pull/123
[#129]: https://github.com/Parallel-in-Time/PFASST/pull/129
[#131]: https://github.com/Parallel-in-Time/PFASST/pull/131
[#132]: https://github.com/Parallel-in-Time/PFASST/pull/132
[#133]: https://github.com/Parallel-in-Time/PFASST/pull/133
[#135]: https://github.com/Parallel-in-Time/PFASST/pull/135
[#138]: https://github.com/Parallel-in-Time/PFASST/pull/138

### Contributors

* Matthew Emmett, Lawrence Berkeley National Laboratory ([memmett][])
* Torbjörn Klatt, Jülich Supercomputing Centre ([torbjoernk][])
* Daniel Ruprecht, Institute of Computational Science, University of Lugano ([danielru][])
* Robert Speck, Jülich Supercomputing Centre ([pancetta][])

[memmett]: https://github.com/memmett
[torbjoernk]: https://github.com/torbjoernk
[danielru]: https://github.com/danielru
[pancetta]: https://github.com/pancetta

---

## v0.2.0 -- MPI PFASST (2014/08/29)

DOI: [10.5281/zenodo.11517](http://dx.doi.org/10.5281/zenodo.11517)

    $ git diff --stat v0.1.0..v0.2.0
    > 52 files changed, 3306 insertions(+), 875 deletions(-)

### Notable Features

* Addition of MPI based PFASST.

### Details

* Addition of MPI based PFASST algorithm using the standard predictor stage in a block mode with 
  fixed iterations.
  ([#46][], [#57][], [#59][])

* Addition of a simple scalar example and appropriate tests.
  ([#61][], [#63][], [#76][], [#81][])

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
[#81]: https://github.com/Parallel-in-Time/PFASST/pull/81

### Contributors

* Matthew Emmett, Lawrence Berkeley National Laboratory ([memmett][])
* Torbjörn Klatt, Jülich Supercomputing Centre ([torbjoernk][])
* Daniel Ruprecht, Institute of Computational Science, University of Lugano ([danielru][])
* Robert Speck, Jülich Supercomputing Centre ([pancetta][])

[memmett]: https://github.com/memmett
[torbjoernk]: https://github.com/torbjoernk
[danielru]: https://github.com/danielru
[pancetta]: https://github.com/pancetta

---

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


[Easylogging++]: https://github.com/easylogging/easyloggingpp
