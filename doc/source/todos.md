# ToDos

- create some sort of traits hierarchy to encapsulate type dependencies and reduce inter-dependencies
  of different components (especially for better mocking)

      template<typename precision>
      struct controller_traits
      {
        typedef precision                                   precision_type;
        typedef encap::VectorEncapsulation<precision_type>  encap_type;
        typedef Sweeper<precision_type, encap_type>         sweeper_type;
        typedef Communicator                                comm_type;
        typedef Status<precision_type>                      status_type;
        // others ...
      };

  use it in the form:

      template<class ControllerTraits>
      class Controller
      {
        // ...
        const shared_ptr<typename ControllerTraits::status_type> get_status() const;
        // ...
      };

  See encapsulation_traits for starting example.

- Sweeper
  - rethink place for pre-/post-predict/sweep methods
    are these better placed in Controller? Controller does the communication, not the sweeper

- Transfer
  - need to write it and put it into controller
  - think about: bind each transfer to two sweepers (maybe with backreferences)?

- Controller
  - see Sweeper::post/pre-predict/sweep
  - advance_time must call advance on all sweepers
    (might need to alter signature of Sweeper::advance)

- Communicator
  - API for creating derived datatypes once
    - therefore read http://accu.org/index.php/journals/2085
    - idea: callback functions in config object for pre-communicator/controller-setup (and
      potentially other places)
  - MPI-speciallization for non-blocking point-to-point and one-sided communication

- check test coverage

- write API documentation

- coding scheme
  - write up guide
    - non-public member variables with leading underscore
    - accessor methods for non-public variables
      - no-logic-setter methods just return reference
      - consistency-checking setter methods prepend with `set_`
      - read-only accessors of the form `const T get_variable() const`
  - cleanup quadrature stuff to comply with new scheme

- Logging
  - reformat logging output (especially DEBUG level)
  - add compile-time-flag to disable logging alltogether
    (and make sure logging calls are not compiled at all)
  - automatic log-formatting of all PFASST types (Encap, Sweeper, Controller, Communicator, Status)
    (especially when they are managed in shared_ptr)
