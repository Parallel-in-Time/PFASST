/**
 * @file pfasst/config.hpp
 * @since v0.3.0
 */
#ifndef _PFASST__CONFIG_HPP_
#define _PFASST__CONFIG_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <map>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#ifdef WITH_MPI
  #include <mpi.h>
#endif


namespace pfasst
{
  static constexpr const char* VERSION = "v0.5.0-rc1-48-g4da93b0-dirty";

  /**
   * @since v0.3.0
   */
  namespace config
  {
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
      // if we're not running under "mpirun/mpiexec", just assume rank 0.
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
        po::options_description all_options;
        map<string, po::options_description> option_groups;
        po::variables_map variables_map;
        vector<string> unrecognized_args;
        bool initialized = false;

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
        po::variables_map& get_variables_map();
        po::options_description& get_all_options();
        vector<string>& get_unrecognized_args();
        //! @}

        //! @{
        /**
         * Adds a new boolean flag.
         *
         * @param[in] group string identifying the parameter group
         * @param[in] option Name of the command line parameter.
         *   It is possible to specify a long and optional short option name by comma-separation.
         *   Short option names are identified by being only a single character.
         *   They are automatically parsed as '`-[SHORT]`' by `boost::program_options` in contrast
         *   to '`--[LONG]`'.
         * @param[in] help help text to be displayed in the help and usage information
         */
        static void add_option(const string& group, const string& option, const string& help);

        /**
         * Adds a new parameter with an expected value of type @p T.
         *
         * @tparam T type of the specified parameter
         * @param[in] group string identifying the parameter group
         * @param[in] option Name of the command line parameter.
         *   It is possible to specify a long and optional short option name by comma-separation.
         *   Short option names are identified by being only a single character.
         *   They are automatically parsed as '`-[SHORT]`' by `boost::program_options` in contrast
         *   to '`--[LONG]`'.
         * @param[in] help help text to be displayed in the help and usage information
         *
         * @overload
         */
        template<typename T>
        static void add_option(const string& group, const string& option, const string& help);
        //! @}

        /**
         * Initialize program options.
         *
         * This initializes `boost::program_options` with all previously added options and groups.
         */
        void init();
    };

    /**
     * get value of specific type @p T
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
      return options::get_instance().get_variables_map()[name].as<T>();
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
      return options::get_instance().get_variables_map().count(name)
              ? options::get_instance().get_variables_map()[name].as<T>() : default_val;
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
      bool no_params_given = options::get_instance().get_variables_map().empty();

      if (!if_no_params || (if_no_params && no_params_given)) {
        stringstream s;
        s << options::get_instance().get_all_options() << endl;
        s << "Logging options:" << endl
          << "  -v [ --verbose ]       activates maximum verbosity" << endl
          << "  --v=arg                activates verbosity upto verbose level `arg`" << endl
          << "                         (valid range: 0-9)" << endl
          << "  -vmodule=arg           actives verbose logging for specific module" << endl
          << "                         (see [1] for details)" << endl << endl
          << "[1]: https://github.com/easylogging/easyloggingpp#vmodule" << endl;
        return s.str();
      } else {
        return string();
      }
    }

    /**
     * Read config parameters from file.
     *
     * A potential INI file might look like the following:
     *
     * @code
     * # Global
     * quiet=yes
     *
     * # Duration
     * dt=0.015625
     * num_steps=2
     * num_iters=10
     *
     * # Tolerances
     * rel_res_tol=1e-10
     *
     * # PFASST
     * num_blocks=1
     * @endcode
     *
     * @param[in] file_name name of the INI-like file containing config parameters;
     *   path/name may be relative
     * @throws invalid_argument if the given file could not be opened
     *
     * @see
     *   [Boost Program Options Documentation on supported INI-like file format]
     *   (http://www.boost.org/doc/libs/1_57_0/doc/html/program_options/overview.html#idp343292240)
     */
    static inline void read_config_file(const string& file_name)
    {
      ifstream ifs(file_name.c_str(), ios_base::in);
      if (!ifs) {
        throw invalid_argument("Config file '" + file_name + "' not found.");
      } else {
        po::store(po::parse_config_file(ifs, options::get_instance().get_all_options()),
                  options::get_instance().get_variables_map());
        po::notify(options::get_instance().get_variables_map());
      }
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
      po::parsed_options parsed = po::command_line_parser(argc, argv)
                                    .options(options::get_instance().get_all_options())
                                    .allow_unregistered().run();
      options::get_instance().get_unrecognized_args() = po::collect_unrecognized(parsed.options,
                                                                                 po::exclude_positional);
      po::store(parsed, options::get_instance().get_variables_map());
      po::notify(options::get_instance().get_variables_map());

      if (options::get_instance().get_variables_map().count("input")) {
        string input_file = config::get_value<string>("input");
        read_config_file(input_file);
      }

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
    }

    /**
     * Initialize options detection and parsing.
     *
     * Prepopulates following groups and parameters:
     *
     * Group      | Parameter     | Type
     * -----------|---------------|---------
     * Global     | `h`, `help`   | `bool`
     * Global     | `q`, `quiet`  | `bool`
     * Global     | `input`       | `string`
     * Global     | `log_prefix`  | `string`
     * Global     | `c`,`nocolor` | `bool`
     * Duration   | `dt`          | `double`
     * Duration   | `tend`        | `double`
     * Duration   | `num_iters`   | `size_t`
     * Tolerances | `abs_res_tol` | `double`
     * Tolerances | `rel_res_tol` | `double`
     */
    static inline void init()
    {
      options::add_option        ("Global", "help,h", "display this help message");
      options::add_option        ("Global", "quiet,q", "don't log to stdout");
      options::add_option<string>("Global", "input", "INI file with configuration options");
      options::add_option<string>("Global", "log_prefix", "a prefix for the log files");
      options::add_option        ("Global", "nocolor,c", "disable colorful logging");

      options::add_option<double>("Duration", "dt", "time step size");
      options::add_option<double>("Duration", "tend", "final time of simulation");
      options::add_option<size_t>("Duration", "num_steps", "number time steps");
      options::add_option<size_t>("Duration", "num_iters", "number of iterations");

      options::add_option<size_t>("Quadrature", "num_nodes", "number of quadrature nodes");

      options::add_option<double>("Tolerances", "abs_res_tol", "absolute residual tolerance");
      options::add_option<double>("Tolerances", "rel_res_tol", "relative residual tolerance");

      options::get_instance().init();
    }
  }  // ::pfasst::config
}  // ::pfasst

#include "pfasst/config_impl.hpp"

#endif  // _PFASST__CONFIG_HPP_
