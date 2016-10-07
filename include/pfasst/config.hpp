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
#include <vector>
using namespace std;

#ifdef WITH_MPI
  #include <mpi.h>
#endif

#ifdef PFASST_USE_LOCAL_CONFIG
  #include "pfasst/site_config.hpp"
#else
  #include "pfasst/version.hpp"
#endif


namespace pfasst
{
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
      private:
        map<string, string> argmap;

        //! @{
        options() {};
        options(const options&) = delete;
        void operator=(const options&) = delete;
        //! @}

      public:
        map<string, string>& get_argument_map() { return argmap; }

        //! @{
        /**
         * Accessor to the singleton instance.
         *
         * @returns singleton config::options instance
         */
        static options& get_instance()
        {
          static options instance;
          return instance;
        }
        // static options& get_instance();
        //! @}
    };

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
      auto& argmap = options::get_instance().get_argument_map();
      if (argmap.count(name) > 0) {
        istringstream ss(argmap.at(name));
        T v;
        ss >> v;
        return v;
      }
      return default_val;
    }

    pair<string,string> split(const string &s) {
      auto i = s.find('=');
      if (i != string::npos) {
        return make_pair(s.substr(0, i), s.substr(i+1, s.size()-i));
      }
      return make_pair("", "");
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
    static inline void read_commandline(int argc, char* argv[])
    {
      auto& argmap = options::get_instance().get_argument_map();
      vector<string> arguments(argv, argv + argc);
      for (auto& arg : arguments) {
        auto kv = split(arg);
        if (get<0>(kv).size() > 0)
            argmap.emplace(kv);
      }
    }

  }  // ::pfasst::config
}  // ::pfasst




#endif  // _PFASST__CONFIG_HPP_
