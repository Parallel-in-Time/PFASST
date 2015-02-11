#ifndef _PFASST__CONFIG_HPP_
#define _PFASST__CONFIG_HPP_

#include <fstream>
#include <string>
#include <map>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;


namespace pfasst
{
  namespace config
  {
    /**
     * @note This uses the Singleton Pattern, and hence pfasst::options::get_instance() is
     *    thread-safe with C++11.
     */
    class options
    {
      public:
        static const size_t LINE_WIDTH = 100;

      private:
        po::options_description all_options;
        map<string, po::options_description> option_groups;
        po::variables_map variables_map;
        vector<string> unrecognized_args;
        bool initialized = false;

        options();
        options(const options&) = delete;
        void operator=(const options&) = delete;

      public:
        static options& get_instance();
        po::variables_map& get_variables_map();
        po::options_description& get_all_options();
        vector<string>& get_unrecognized_args();
        static void add_option(const string& group, const string& option, const string& help);

        template<typename T>
        static void add_option(const string& group, const string& option, const string& help);

        void init();
    };

    template<typename T>
    static T get_value(const string& name, const T& default_val);

    template<typename T>
    static T get_value(const string& name);

    /**
     * @returns empty string if params are set and `if_no_params` is `true`
     */
    static string print_help(bool if_no_params = false);

    static inline void read_commandline(int argc, char* argv[], bool exit_on_help = true)
    {
      po::parsed_options parsed = po::command_line_parser(argc, argv)
                                    .options(options::get_instance().get_all_options())
                                    .allow_unregistered().run();
      options::get_instance().get_unrecognized_args() = po::collect_unrecognized(parsed.options,
                                                                                 po::exclude_positional);
      po::store(parsed, options::get_instance().get_variables_map());
      po::notify(options::get_instance().get_variables_map());

      if (options::get_instance().get_variables_map().count("help")) {
        cout << print_help() << endl;
        if (exit_on_help) exit(0);
      }
    }

    /**
     * @throws invalid_argument if the given file could not be opened
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

    static inline void init()
    {
      options::add_option("Global", "help,h", "display this help message");

      options::add_option<double>("Duration", "dt", "time step size");
      options::add_option<double>("Duration", "tend", "final time of simulation");
      options::add_option<size_t>("Duration", "num_iters", "number of iterations");

      options::add_option<double>("Tolerances", "abs_res_tol", "absolute residual tolerance");
      options::add_option<double>("Tolerances", "rel_res_tol", "relative residual tolerance");

      options::get_instance().init();
    }
  }  // ::pfasst::config
}  // ::pfasst

#include "pfasst/config_impl.hpp"

#endif  // _PFASST__CONFIG_HPP_
