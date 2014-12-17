# Boris SDC                                                                   {#page_examples_boris}

This directory contains an implementations of a modification to the Boris method to solve second order ODEs with the 
Velocity-Verlet scheme using the PFASST framework.

The sweeper with it's _Boris magic_ is implemented in `boris_sweeper.hpp`.
The physical properties of a testbed example with a panning trap are defined in `physics.hpp` and `simple_physics.hpp`,
while the data structure for the particles are defined in `particle.hpp` and `particle_3d.hpp`.
