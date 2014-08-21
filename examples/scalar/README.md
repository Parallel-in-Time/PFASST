# Scalar                                                                     {#page_examples_scalar}

This directory contains a simple implementations of scalar ODE solver using the PFASST framework.

The SDC sweeper is defined in `scalar_sweeper.hpp`.
As this simple equation does not require any special handling, all the SDC magic is derived and 
taken from pfasst::encap::IMEXSweeper::sweep().
