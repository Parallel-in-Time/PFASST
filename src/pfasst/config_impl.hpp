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

    /**
     * @todo Make config::options::add_option() fail when called called after
     *   config::options::init().
     */
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

    /**
     * @todo Make config::options::add_option() fail when called called after
     *   config::options::init().
     */
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

    /**
     * @todo Make config::options::init() fail when called twice.
     */
    void options::init()
    {
      if (!this->initialized) {
        for (auto const& kv : this->option_groups) {
          this->all_options.add(kv.second);
        }
      }
      this->initialized = true;
    }
  }  // ::pfasst::config
}  // ::pfasst
