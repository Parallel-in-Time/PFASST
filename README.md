PFASST                                                                                   {#mainpage}
======

The PFASST project is a C++ implementation of the *parallel full approximation scheme in space and
time* (PFASST) algorithm, which in turn is a time-parallel algorithm for solving ODEs and PDEs.  It
also contains basic implementations of the *spectral deferred correction* (SDC) and *multi-level
spectral deferred correction* (MLSDC) algorithms.


News
----

* May XX, 2015: PFASST v0.5.0 release. Please see the [release notes](#releases) for more
  information.

* April 17, 2015: PFASST v0.4.0 release. Please see the [release notes](#releases) for more
  information.

* December 12, 2014: PFASST v0.3.0 release. Please see the [release notes](#releases) for more
  information.

* August 29, 2014: PFASST v0.2.0 release. Please see the [release notes](#releases) for more
  information.

* July 25, 2014: PFASST v0.1.0 released. Please see the [release notes](#releases) for more
  information.


References
----------

See \subpage #citelist .


Documentation
-------------

Doxygen generated documentation can be found [on the PinT server][documentation].
Currently, it features the following content:

* \subpage #page_building_installing
* \ref Examples
* \subpage #page_contributing
* \subpage #page_style_guide
* \subpage #page_troubleshooting

[documentation]: https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST_LATEST_STABLE_DOCU/doxygen


Releases
--------

* **v0.5.0** The MPI Bugfix Release (2015/05/XX)

  A few important fixes to MPI behaviour.  
  DOI: [][DOI_v050].  
  See [the Changelog](#page_changelog) for details.

* **v0.4.0** The Big Improvement (2015/04/17)

  A bunch of improvements and usability enhancements.  
  DOI: [10.6084/m9.figshare.1381721][DOI_v040].  
  See [the Changelog](#page_changelog) for details.

* **v0.3.0** The Big Cleanup Release (2014/12/12)

  A lot of internal cleanup and usability enhancements.  
  DOI: [10.5281/zenodo.13221][DOI_v030].  
  See [the Changelog](#page_changelog) for details.

* **v0.2.0** MPI-PFASST Release (2014/08/29)

  Implementation of PFASST with MPI.  
  DOI: [10.5281/zenodo.11517][DOI_v020].  
  See [the Changelog](#page_changelog) for details.

* **v0.1.0** First Release (2014/07/25)

  Initial release with basic implementations of SDC and MLSDC.  
  DOI: [10.5281/zenodo.11047][DOI_v010].  
  See [the Changelog](#page_changelog) for details.

[DOI_v010]: http://dx.doi.org/10.5281/zenodo.11047
[DOI_v020]: http://dx.doi.org/10.5281/zenodo.11517
[DOI_v030]: http://dx.doi.org/10.5281/zenodo.13221
[DOI_v040]: http://dx.doi.org/10.6084/m9.figshare.1381721
[DOI_v050]: 

Release tags will be signed by one of the following PGP keys:

    0xAD9F8DC6 2014-12-12 Torbjörn Klatt <t.klatt@fz-juelich.de>
    Fingerprint 6277 E9D9 7AA2 DBBE 5DE7 7498 756B F4D7 AD9F 8DC6

    1024D/9950EF2E 2001-11-22 Matthew Emmett <matt@emmett.ca>
    Fingerprint B09C 1425 1C47 B58E B3AC 2E74 6F22 8460 9950 EF2E


Build status
------------

| Branch      | Status                 | Details                           |
|-------------|------------------------|-----------------------------------|
| Master      | ![master-status-img][] | [Jenkins Job][master-status-link] |
| Development | ![dev-status-img][]    | [Jenkins Job][dev-status-link]    |

Results of static code analysis such as test coverage or SLOCcount metrics can be found at
[this Jenkins job][coverage-job].
Further metrics can be found on [Open Hub aka Ohloh][openhub].

[master-status-img]:  https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST_LATEST_STABLE/badge/icon
[master-status-link]: https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST_LATEST_STABLE
[dev-status-img]:     https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST_DEVELOPMENT/badge/icon
[dev-status-link]:    https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST_DEVELOPMENT
[coverage-job]:       https://pint.fz-juelich.de/ci/view/PFASST/job/PFASST_LATEST_STABLE_ANALYSIS/
[openhub]:            https://www.openhub.net/p/PFASST
