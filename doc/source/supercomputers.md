# Usage on Supercomputers                                                     {#page_supercomputers}

## JUQUEEN

There are a few steps one has to complete before using _PFASST++_ on _JUQUEEN_.

1. make sure the only modules loaded are `gcc/4.8.1` and `fftw3/3.3.3`

2. prepare a custom installation of _Boost_

   Because at time of writing only 1.47.0 is available compiled with _XL/C++_, which does not
   support the C++11 features we learned to love.

  1. download a recent version of _Boost_ (i.e. 1.57.0)

  2. extract the archive and `cd` into it

  3. run 

         CXX=`which mpig++` CC=`which mpigcc` ./bootstrap.sh \
           --with-libraries=program_options --prefix=<INSTALL_PLACE>

     where `INSTALL_PLACE` is something like `$HOME/progs/juqueen`

  4. create file called `user-config.jam` with the following content:

         using gcc : 4.8.1 : mpig++ ;

  5. now compile

         CXX=`which mpig++` CC=`which mpigcc` ./b2 --prefix=<INSTALL_PLACE> \
           --reconfigure link=static stage

     and install

         CXX=`which mpig++` CC=`which mpigcc` ./b2 --prefix=<INSTALL_PLACE> \
           link=static install

3. go to the sources of _PFASST++_, create a build folder and step into it

4. run _CMake_

       cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchain_juqueen.cmake \
         -DCMAKE_CXX_COMPILER=`which mpig++` -DCMAKE_C_COMPILER=`which mpigcc` \
         -DBOOST_ROOT=<INSTALL_PLACE> -Dpfasst_WITH_MPI=ON \
         -Dpfasst_BUILD_TESTS=OFF -Dpfasst_BUILD_SHARED_LIBS=OFF ..

5. run make

       make

6. write your _LLSubmit_ configuration file

   for example for the Advection-Diffusion example:

       #@job_name              = MPI_PFASST_test
       #@comment               = "a little test of PFASST with MPI"
       #@output                = mpi_pfasst_test_$(jobid)_$(stepid).out
       #@error                 = mpi_pfasst_test_$(jobid)_$(stepid).err
       #@environment           = COPY_ALL
       #@job_type              = bluegene
       #@notification          = never
       #@wall_clock_limit      = 00:10:00
       #@bg_size               = 1
       #@bg_connectivity       = TORUS
       #@queue
       runjob --np 32 --ranks-per-node 32 : \
         <PATH_TO_BUILD_DIR>/examples/advection_diffusion/mpi_pfasst \
         -q --tend 0.64 --dt 0.01 --num_iter 8

## Edison

There are a few steps one has to complete before using _PFASST++_ on
_Edison_.  These steps have been tested with the GNU programming
environment (that is, with the `PrgEnv-gnu` module).

1. Load the following modules:

       module load cmake python fftw eigen3 boost

2. Go to the sources of _PFASST++_, create a `build` folder and step into it:

       cd PFASST
       mkdir build
       cd build

3. Run _CMake_

       cmake -Dpfasst_BUILD_TESTS=OFF -Dpfasst_BUILD_SHARED_LIBS=OFF ..

4. Run _make_

       make -j 4

5. Run an example in an interactive job:

       qsub -I -q debug -l mppwidth=4
       cd PFASST/build/examples/advection_diffusion
       aprun -n 4 ./mpi_pfasst -c

