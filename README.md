PFASST                                                                                  {#mainpage}
======

The PFASST algorithm is a time-parallel algorithm for solving ODEs and
PDEs.

The PFASST project is a C++ implementation of the 'parallel full
approximation scheme in space and time' (PFASST) algorithm.  It also
contains basic implementations of the 'spectral deferred correction'
(SDC) and 'multi-level spectral deferred correction' (MLSDC)
algorithms.


References
----------

XXX


Documentation
-------------

Doxygen generated documentation can be found [on the PinT server][documentation].  
Currently, it features the following content:

* \subpage #page_building_installing
* \subpage #page_examples
* \subpage #page_contributing
* \subpage #page_style_guide


Build status
------------

| Branch      | Status                              |
|-------------|-------------------------------------|
| Master      | ![master-status][]                  |
| Development | ![dev-status][]                     |

For details see [Travis][pfasst-travis].


[documentation]: https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST%20%28Docu%29/doxygen/
[pfasst-travis]: https://travis-ci.org/Parallel-in-Time/PFASST
[master-status]: https://travis-ci.org/Parallel-in-Time/PFASST.svg?branch=master
[dev-status]:    https://travis-ci.org/Parallel-in-Time/PFASST.svg?branch=development
