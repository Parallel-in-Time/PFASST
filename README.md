PFASST                                                                                   {#mainpage}
======

The PFASST project is a C++ implementation of the *parallel full approximation scheme in space and
time* (PFASST) algorithm, which in turn is a time-parallel algorithm for solving ODEs and PDEs.  It
also contains basic implementations of the *spectral deferred correction* (SDC) and *multi-level
spectral deferred correction* (MLSDC) algorithms.


News
----

* August 29, 2014: PFASST v0.2.0 release. Please sse the [release notes](#releases) for more 
  information.

* July 25, 2014: PFASST v0.1.0 released.  Please see the [release notes](#releases) for more
  information.


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
* \subpage #page_troubleshooting


Releases
--------

* **v0.2.0** MPI-PFASST Release (2014/08/29)

  Implementation of PFASST with MPI.
  DOI: [10.5281/zenodo.11047][DOI_v010]
  See \subpage #page_changelog "the Changelog" for details.

* **v0.1.0** First Release (2014/07/25)

  Initial release with basic implementations of SDC and MLSDC.
  DOI: [10.5281/zenodo.11047][DOI_v010]
  See \subpage #page_changelog "the Changelog" for details.

[DOI_v010]: http://dx.doi.org/10.5281/zenodo.11047

Release tags will be signed by one of the following PGP keys:

    0x9CF9601F 2011-07-28 Torbjörn Klatt <t.klatt@fz-juelich.de>
    Fingerprint DB8D EA65 F6A7 3DE0 E7EA F607 6CE8 B4B1 9CF9 601F


Build status

| Branch      | Status                              |
|-------------|-------------------------------------|
| Master      | ![master-status][]                  |
| Development | ![dev-status][]                     |

For details see [Travis][pfasst-travis].

A test coverage report is generated and published by [Jenkins][coverage-report].


[documentation]:   https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST%20%28Docu%29/doxygen/
[pfasst-travis]:   https://travis-ci.org/Parallel-in-Time/PFASST
[master-status]:   https://travis-ci.org/Parallel-in-Time/PFASST.svg?branch=master
[dev-status]:      https://travis-ci.org/Parallel-in-Time/PFASST.svg?branch=development
[coverage-report]: https://pint.fz-juelich.de/ci/job/PFASST_GCC_4-8_Coverage/Test_Coverage
