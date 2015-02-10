#include "pfasst/config.hpp"

#include <fstream>
#include <sstream>
#include <utility>
using namespace std;


namespace pfasst
{
  namespace config
  {
    options::options()
    {}

    options& options::get_instance()
    {
      static options instance;
      return instance;
    }

    po::variables_map& options::get_variables_map()
    {
      return this->variables_map;
    }

    po::options_description& options::get_all_options()
    {
      return this->all_options;
    }

    vector<string>& options::get_unrecognized_args()
    {
      return this->unrecognized_args;
    }

    void options::add_option(const string& group, const string& option, const string& help)
    {
      auto& opts = get_instance();
      opts.option_groups.emplace(make_pair<string,
                                 po::options_description>(string(group),
                                                          po::options_description(string(group),
                                                                                  LINE_WIDTH)));
      opts.option_groups[group].add_options()
                                  (option.c_str(), help.c_str());
    }

    template<typename T>
    void options::add_option(const string& group, const string& option, const string& help)
    {
      auto& opts = get_instance();

      opts.option_groups.emplace(make_pair<string,
                                 po::options_description>(string(group),
                                                          po::options_description(string(group),
                                                                                  LINE_WIDTH)));
      opts.option_groups[group].add_options()
                                  (option.c_str(), po::value<T>(), help.c_str());
    }

    void options::init()
    {
      if (!this->initialized) {
        for (auto const & kv : this->option_groups) {
          this->all_options.add(kv.second);
        }
      }
      this->initialized = true;
    }


    template<typename T>
    static T get_value(const string& name, const T& default_val)
    {
      return options::get_instance().get_variables_map().count(name)
              ? options::get_instance().get_variables_map()[name].as<T>() : default_val;
    }

    template<typename T>
    static T get_value(const string& name)
    {
      return options::get_instance().get_variables_map()[name].as<T>();
    }

    static string print_help(bool if_no_params)
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
  }  // ::pfasst::config
}  // ::pfasst
