#ifndef _PFASST__LOGGING_HPP_
#define _PFASST__LOGGING_HPP_

#include "site_config.hpp"

#include <ctime>
#include <string>
using namespace std;

#include "config.hpp"

// enable easy logging of STL containers
#define _ELPP_STL_LOGGING
// disable creation of default log file
#define _ELPP_NO_DEFAULT_LOG_FILE
// enable passing `--logging-flags` via command line
#define _ELPP_LOGGING_FLAGS_FROM_ARG

#ifndef NDEBUG
  #define _ELPP_DEBUG_ASSERT_FAILURE
  #define _ELPP_STACKTRACE_ON_CRASH
#endif

#include <easylogging++.h>

#ifndef PFASST_LOGGER_INITIALIZED
  // initialize easyloggingpp
  // FIXME: this might already be called by code using PFASST++
  _INITIALIZE_EASYLOGGINGPP
  #define PFASST_LOGGER_INITIALIZED
#endif

namespace pfasst {
  namespace log {

    static void load_default_config()
    {
      el::Loggers::configureFromGlobal(PFASST_LOGGING_DEFAULT_CONF_FILE);
    }

    static void set_logging_flags()
    {
      el::Loggers::addFlag(el::LoggingFlag::NewLineForContainer);
      el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
      el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
      el::Loggers::addFlag(el::LoggingFlag::MultiLoggerSupport);
      el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);
      el::Loggers::addFlag(el::LoggingFlag::AutoSpacing);
    }

    static void set_logfile_from_options()
    {
      time_t now = time(nullptr);
      char time_str[31];
      strftime(time_str, sizeof(time_str), "%Y-%m-%dT%H-%M-%S%z", localtime(&now));
      el::Loggers::reconfigureLogger("data_values",
                                     el::ConfigurationType::Filename,
                                     pfasst::config::get_value<string>("data_values_out",
                                                                       string(time_str) + "_data_values.log"));
    }

    static void init_config_options(po::options_description& opts)
    {
      opts.add_options()
        ("data_values_out", po::value<string>(), "name of file to write 'data_values' to")
        ;
    }

    static void enable_config_options(size_t index = -1)
    {
      pfasst::config::Options::get_instance()
        .register_init_function("Logger",
                                function<void(po::options_description&)>(pfasst::log::init_config_options),
                                index);
    }
  }  // ::pfasst::log
}  // ::pfasst

#define PFASST_START_LOG(argc, argv)\
  _START_EASYLOGGINGPP(argc, argv);\
  pfasst::log::load_default_config();\
  pfasst::log::set_logging_flags();\
  pfasst::log::set_logfile_from_options();

#endif  // _PFASST__LOGGING_HPP_
