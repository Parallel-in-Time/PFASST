#ifndef _PFASST__CONFIG_HPP_
#define _PFASST__CONFIG_HPP_

#include <algorithm>
#include <cassert>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
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

        options() {}
        options(const options&) = delete;
        void operator=(const options&) = delete;

      public:
        static options& get_instance()
        {
          static options instance;
          return instance;
        }

        po::variables_map& get_variables_map()
        {
          return this->variables_map;
        }

        po::options_description& get_all_options()
        {
          return this->all_options;
        }

        vector<string>& get_unrecognized_args()
        {
          return this->unrecognized_args;
        }

        static void add_option(string group, string option, string help)
        {
          auto& opts = get_instance();
          opts.option_groups.emplace(make_pair<string,
                                     po::options_description>(string(group),
                                                              po::options_description(string(group), LINE_WIDTH)));
          opts.option_groups[group].add_options()
            (option.c_str(), help.c_str());
        }

        template<typename T>
        static void add_option(string group, string option, string help)
        {
          auto& opts = get_instance();

          opts.option_groups.emplace(make_pair<string,
                                     po::options_description>(string(group),
                                                              po::options_description(string(group), LINE_WIDTH)));
          opts.option_groups[group].add_options()
            (option.c_str(), po::value<T>(), help.c_str());
        }

        void init()
        {
          if (!this->initialized) {
            for (auto const & kv : this->option_groups) {
              this->all_options.add(kv.second);
            }
          }
          this->initialized = true;
        }

    };

    template<typename T>
    inline static T get_value(const string& name, const T& default_val)
    {
      return options::get_instance().get_variables_map().count(name)
              ? options::get_instance().get_variables_map()[name].as<T>() : default_val;
    }

    template<typename T>
    inline static T get_value(const string& name)
    {
      return options::get_instance().get_variables_map()[name].as<T>();
    }

    /**
     * @returns empty string if params are set and `if_no_params` is `true`
     */
    inline static string print_help(bool if_no_params = false)
    {
      bool no_params_given = options::get_instance().get_variables_map().empty();

      if (!if_no_params || (if_no_params && no_params_given)) {
        stringstream s;
        s << options::get_instance().get_all_options() << endl;
        s << "Logging options:" << endl
          << "  -v [ --verbose ]       activates maximum verbosity" << endl
          << "  --v=arg                activates verbosity upto verbose level 2" << endl
          << "                         (valid range: 0-9)" << endl
          << "  -vmodule=arg           actives verbose logging for specific module" << endl
          << "                         (see [1] for details)" << endl;
        s << "[1]: https://github.com/easylogging/easyloggingpp#vmodule" << endl;
        return s.str();
      } else {
        return string();
      }
    }

    inline static void read_commandline(int argc, char* argv[], bool exit_on_help = true)
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
    inline static void read_config_file(string file_name)
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

    inline static void init()
    {
      options::add_option("Global", "help,h", "display this help message");

      options::add_option<double>("Duration", "dt", "time step size");
      options::add_option<double>("Duration", "tend", "final time of simulation");
      options::add_option<size_t>("Duration", "num_iters", "number of iterations");

      options::get_instance().init();
    }

  }  // ::pfasst::config

}  // ::pfasst

#endif  // _PFASST__CONFIG_HPP_
