# Trouble Shooting                                                           {#page_troubleshooting}

## I'm getting SEGFAULTs, how do I track them down?

Please try running your program through GDB and/or valgrind.

## How do I run my MPI program through GDB?

To run your MPI program through GDB:

    mpirun -n 4 xterm -e gdb -x my-gdb-commands my-pfasst-program

where `my-gdb-commands` is a plain text file containing GDB commands
to run your program, and `my-pfasst-program` is your executable.  A
very simple set of GDB commands for the `my-gdb-commands` is:

    catch throw
    run

The result is: `mpirun` will launch 4 `xterm` windows, each of which
will immediately run your program through `gdb`, which subsequently
runs your program.  Hopefully GDB will catch the SEGFAULT.

## How do I run my MPI program through valgrind?

To run your MPI program through valgrind:

    mpirun -n 4 xterm -e valgrind my-pfasst-program







