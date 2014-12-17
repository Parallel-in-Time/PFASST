# Van-der-Pol Oscillator                                                  {#page_examples_vanderpol}

This directory contains an implementations of solver for the Van-der-Pol oscillator using the PFASST framework.

The SDC sweeper is defined in `vdp_sweeper.hpp`.
As this simple equation does not require any special handling, all the SDC magic is derived and 
taken from pfasst::encap::IMEXSweeper::sweep().
