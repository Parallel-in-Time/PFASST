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
     * @note This is using the Singleton Pattern.
     *    pfasst::Options::get_instance() is thread-safe with C++11.
     */
    class Options
    {
      public:
        //! Maximum width of a line in the help output
        static const size_t LINE_WIDTH = 100;

      private:
        po::options_description all_options;
        vector<string> group_names;
        map<string, po::options_description> option_groups;
        map<string, function<void(po::options_description&)>> init_functions;
        po::variables_map vars_map;
        vector<string> unrecognized_args;
        bool initialized = false;

      private:
        //! ctor is private to make it a singleton
        Options()
        {}
        Options(const Options&) = delete;
        void operator=(const Options&) = delete;

      public:
        static Options& get_instance()
        {
          static Options instance;
          return instance;
        }

        po::variables_map& get_variables_map()
        { return this->vars_map; }

        po::options_description& get_all_options()
        { return this->all_options; }

        vector<string>& get_unrecognized_args()
        { return this->unrecognized_args; }

        /**
         * @throws invalid_argument if there is already an init function or options group with the 
         *     given name registered
         */
        void register_init_function(const string& name,
                                    const function<void(po::options_description&)>& fnc,
                                    int index = -1)
        {
          assert(fnc);
          if (!this->option_groups.count(name) && !this->init_functions.count(name)) {
            this->option_groups.emplace(make_pair<string,
                                                  po::options_description>(string(name),
                                                                           po::options_description(name,
                                                                                                   LINE_WIDTH)));
            this->init_functions[name] = fnc;
            if (index == -1) {
              this->group_names.push_back(name);
            } else {
              this->group_names.insert(this->group_names.begin() + index, name);
            }
          } else {
            throw invalid_argument("There is already an init function named '" + name + "' registered.");
          }
        }

        void init()
        {
          if (!this->initialized) {
            for(string name : this->group_names) {
              auto fnc = this->init_functions[name];
              // calling init function
              fnc(this->option_groups[name]);
              this->all_options.add(this->option_groups[name]);
            }
            this->initialized = true;
          }
        }
    };


    template<typename T>
    static const T get_value(const string& name, const T& default_val)
    {
      return Options::get_instance().get_variables_map().count(name) 
              ? Options::get_instance().get_variables_map()[name].as<T>() : default_val;
    }


    template<typename T>
    static const T get_value(const string& name)
    {
      return Options::get_instance().get_variables_map()[name].as<T>();
    }


    static bool no_params_given()
    {
      return Options::get_instance().get_variables_map().empty();
    }


    static string pretty_print()
    {
      stringstream s;
      s << Options::get_instance().get_all_options();
      return s.str();
    }


    /**
     * @returns empty string if params are set and `if_no_params` is `true`
     */
    static string print_help(bool if_no_params = false)
    {
      if (!if_no_params || (if_no_params && no_params_given())) {
        return pretty_print();
      } else {
        return string();
      }
    }


    static void init_global_options(po::options_description& opts)
    {
      opts.add_options()
        ("help,h", "display this help message");
    }


    static void init_config()
    {
      Options::get_instance()
        .register_init_function("Global Options",
                                function<void(po::options_description&)>(init_global_options),
                                0);

      Options::get_instance().init();
    }


    static void read_commandline(int argc, char* argv[], bool exit_on_help = true)
    {
      po::parsed_options parsed = po::command_line_parser(argc, argv)
                                    .options(Options::get_instance().get_all_options())
                                    .allow_unregistered().run();
      Options::get_instance().get_unrecognized_args() = po::collect_unrecognized(parsed.options,
                                                                                 po::exclude_positional);
      po::store(parsed, Options::get_instance().get_variables_map());
      po::notify(Options::get_instance().get_variables_map());

      if (Options::get_instance().get_variables_map().count("help")) {
        cout << print_help() << endl;
        if (exit_on_help) exit(0);
      }
    }


    /**
     * @throws invalid_argument if the given file could not be opened
     */
    static void read_config_file(string file_name)
    {
      ifstream ifs(file_name.c_str(), ios_base::in);
      if (!ifs) {
        throw invalid_argument("Config file '" + file_name + "' not found.");
      } else {
        po::store(po::parse_config_file(ifs, Options::get_instance().get_all_options()),
                  Options::get_instance().get_variables_map());
        po::notify(Options::get_instance().get_variables_map());
      }
    }

  }  // ::pfasst::config
}  // ::pfasst

#endif  // _PFASST__CONFIG_HPP_
