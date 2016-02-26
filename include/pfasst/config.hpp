/**
 * @file pfasst/config.hpp
 * @since v0.3.0
 */
#ifndef _PFASST__CONFIG_HPP_
#define _PFASST__CONFIG_HPP_

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
using namespace std;

#ifdef WITH_MPI
  #include <mpi.h>
#endif

namespace pfasst
{
  /**
   * @since v0.3.0
   */
  namespace config
  {

    template<typename T>
    T string_as(const std::string s)
    {
      // XXX: use boost lexical_cast instead?
      T v;
      std::istringstream stream(s);
      stream >> v;
      return v;
    }

    /**
     * Get MPI rank during initialization.
     *
     * When running without MPI (ie, without using mpirun/mpiexec), returns 0.  When running with
     * MPI, returns the MPI rank.
     *
     * If the user is running with MPI and MPI_Init hasn't been called yet, this will return 0.  I
     * hope this is rare.
     */
    int get_rank()
    {
#ifdef WITH_MPI
      int initialized = 0, rank = 0;
      MPI_Initialized(&initialized);
      if (initialized) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      }
      return rank;
#else
      return 0;
#endif
    }

    /**
     * Runtime config options provider.
     *
     * This singleton provides easy access to command line parameters at runtime.
     *
     * @note This uses the Singleton Pattern, and hence pfasst::options::get_instance() is
     *    thread-safe with C++11.
     * @since v0.3.0
     * @ingroup Internals
     */
    class options
    {
      public:
        /// line width of help and usage information
        static const size_t LINE_WIDTH = 100;

      private:
        bool initialized = false;

        typedef map<string, string> option_values_map;
        typedef map<string, int> option_names_map;

        option_values_map option_values;
        option_names_map option_names;

        //! @{
        options();
        options(const options&) = delete;
        void operator=(const options&) = delete;
        //! @}

      public:
        //! @{
        /**
         * Accessor to the singleton instance.
         *
         * @returns singleton config::options instance
         */
        static options& get_instance();
        //! @}

        option_values_map& get_option_values() { return this->option_values; }
        option_names_map& get_option_names() { return this->option_names; }

        //! @{
        /**
         * Adds a new option.
         *
         * @param[in] group Name of the parameter group.
         * @param[in] option Name of the command line parameter.
         * @param[in] help Help message to be displayed in help and usage information.
         * @param[in] boolean True if option is boolean.
         */
        static void add_option(const string& group, const string& option, const string& help, bool boolean=false) {
          auto& names = options::get_instance().get_option_names();
          names.emplace(string(option), boolean);
        }

        /**
         * Initialize program options.
         *
         * This initializes `boost::program_options` with all previously added options and groups.
         */
        void init();
    };

    /**
     * Get value of specific type @p T
     *
     * @tparam T type of the retreived value
     * @param[in] name Name of the (long) option as defined with @p option in options::add_option()
     * @returns Value of type @p T of desired option.
     * @throws boost::bad_any_cast if option is not given.
     *
     * @see
     *   [boost::any::bad_any_cast]
     *   (http://www.boost.org/doc/libs/1_57_0/doc/html/boost/bad_any_cast.html)
     */
    template<typename T>
    inline T get_value(const string& name)
    {
      auto& opts = options::get_instance().get_option_values();
      return string_as<T>(opts.at(name));
    }

    /**
     * Get value of specific type @p T with default value.
     *
     * @tparam T type of the retreived value
     *
     * @overload
     */
    template<typename T>
    inline T get_value(const string& name, const T& default_val)
    {
      auto& opts = options::get_instance().get_option_values();
      return opts.count(name) > 0 ? string_as<T>(opts.at(name)) : default_val;
    }

    /**
     * Compile basic help and usage information.
     *
     * Depending on @p if_no_params and presence of given command line parameters the help and
     * usage information is compiled.
     * In case @p if_no_params is `true` and there are no parameters given on the command line or
     * @p if_no_params is `false` no matter whether parameters are given, the help message is
     * generated.
     *
     * @param[in] if_no_params flag governing compilation of help and usage information
     * @returns string containing basic help and usage information; string may be empty
     */
    static string print_help(bool if_no_params = false)
    {
      // bool no_params_given = options::get_instance().get_variables_map().empty();

      //      if (!if_no_params || (if_no_params && no_params_given)) {
        // stringstream s;
        // // s << options::get_instance().get_all_options() << endl;
        // s << "Logging options:" << endl
        //   << "  -v [ --verbose ]       activates maximum verbosity" << endl
        //   << "  --v=arg                activates verbosity upto verbose level `arg`" << endl
        //   << "                         (valid range: 0-9)" << endl
        //   << "  -vmodule=arg           actives verbose logging for specific module" << endl
        //   << "                         (see [1] for details)" << endl << endl
        //   << "[1]: https://github.com/easylogging/easyloggingpp#vmodule" << endl;
        // return s.str();
      //      } else {
        return string();
        //      }
    }

    /**
     * Read and parse command line parameters.
     *
     * @param[in] argc Number of command line arguments as provided by `%main(int, char**)`.
     * @param[in] argv List of command line arguments as provided by `%main(int, char**)`.
     * @param[in] exit_on_help Whether to exit the program after displaying help and usage
     *   information.
     *
     * @note Will call `std::exit` in case the `help` parameter has been provided and
     *   @p exit_on_help is `true`.
     */
    static inline void read_commandline(int argc, char* argv[], bool exit_on_help = true)
    {
      int i = 1;

      auto& names  = options::get_instance().get_option_names();
      auto& values = options::get_instance().get_option_values();

      while (i < argc) {
        string name = string(argv[i]);
        if (names.count(name)) {
          auto boolean = names.at(name);
          if (boolean) {
            values.emplace(name, "true");
            i += 1;
          } else {
            string value = string(argv[i+1]);
            values.emplace(name, value);
            i += 2;
          }
        } else {
          i += 1;
        }
      }

      /* XXX: handle "help"... */

      /*
      if (options::get_instance().get_variables_map().count("help")) {
        if (get_rank() == 0) {
          cout << print_help() << endl;
        }

        if (exit_on_help) {
#ifdef WITH_MPI
          MPI_Finalize();
#endif
          exit(0);
        }
      }
       */
    }

    /**
     * Initialize options detection and parsing.
     */
    static inline void init()
    {
      options::add_option("Global", "help", "display this help message", true);
      options::add_option("Global", "quiet", "don't log to stdout", true);
      options::add_option("Global", "nocolor", "disable colorful logging", true);

      options::add_option("Duration", "dt", "time step size");
      options::add_option("Duration", "tend", "final time of simulation");
      options::add_option("Duration", "num_steps", "number time steps");
      options::add_option("Duration", "num_iters", "number of iterations");

      options::add_option("Quadrature", "num_nodes", "number of quadrature nodes");

      options::add_option("Tolerances", "abs_res_tol", "absolute residual tolerance");
      options::add_option("Tolerances", "rel_res_tol", "relative residual tolerance");

      options::get_instance().init();
    }
  }  // ::pfasst::config
}  // ::pfasst

#include "pfasst/config_impl.hpp"

#endif  // _PFASST__CONFIG_HPP_
