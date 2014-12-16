#ifndef _PFASST__LOGGING_HPP_
#define _PFASST__LOGGING_HPP_

#include "site_config.hpp"

#include <string>
using namespace std;

struct OUT
{
  public:
    static const string black;
    static const string red;
    static const string green;
    static const string yellow;
    static const string blue;
    static const string magenta;
    static const string cyan;
    static const string white;

    static const string bold;
    static const string underline;

    static const string reset;
};

const string OUT::black = "\033[30m";
const string OUT::red = "\033[31m";
const string OUT::green = "\033[32m";
const string OUT::yellow = "\033[33m";
const string OUT::blue = "\033[34m";
const string OUT::magenta = "\033[35m";
const string OUT::cyan = "\033[36m";
const string OUT::white = "\033[37m";
const string OUT::bold = "\033[1m";
const string OUT::underline = "\033[4m";
const string OUT::reset = "\033[0m";


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

#ifdef PFASST_NO_LOGGING
  #define _ELPP_DISABLE_LOGS
#endif

#include <pfasst/easylogging++.h>


#ifndef PFASST_LOGGER_INITIALIZED
  // initialize easyloggingpp
  // FIXME: this might already be called by code using PFASST++
  _INITIALIZE_EASYLOGGINGPP
  /**
   * guard symbol to ensure easylogging++ is only initialized once
   *
   * When this symbol is defined, it is expected that `_INITIALIZE_EASYLOGGINGPP` has been called once to initialize
   * easylogging++.
   *
   * \note In case the executable using PFASST++ is also using easylogging++ as its logging library and initializing it
   *     prior to including the PFASST++ headers, please define this symbol prior including any of the PFASST++ headers.
   */
  #define PFASST_LOGGER_INITIALIZED
#endif

#ifndef PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH
  //! precision of milliseconds to be printed
  #define PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH "4"
#endif

#ifndef VLOG_FUNC_START
  #define VLOG_FUNC_START(scope) \
    pfasst::log::stack_position++; \
    LOG(DEBUG) << std::string((pfasst::log::stack_position - 1) * 2, ' ') << "START:" << std::string(scope) + "::" + std::string(__func__) + "() "
#endif

#ifndef VLOG_FUNC_END
  #define VLOG_FUNC_END(scope) \
    pfasst::log::stack_position--; \
    LOG(DEBUG) << std::string(pfasst::log::stack_position * 2, ' ') << "DONE: " << std::string(scope) + "::" + std::string(__func__) + "()";
#endif

#ifndef LOG_PRECISION
  #define LOG_PRECISION 5
#endif

#define VLOG_INDENT(level) \
  VLOG(level) << string(pfasst::log::stack_position * 2, ' ')


namespace pfasst
{
  namespace log
  {
    static size_t stack_position;

    /**
     * sets default configuration for default loggers
     */
    inline static void load_default_config()
    {
      const string TIMESTAMP = OUT::white + "%datetime{%H:%m:%s,%g}" + OUT::reset + " ";
      const string LEVEL = "[%level]";
      const string VLEVEL = "[VERB%vlevel]";
      const string POSITION = "%fbase:%line";
      const string MESSAGE = "%msg";

      el::Configurations defaultConf;
      defaultConf.setToDefault();

      defaultConf.setGlobally(el::ConfigurationType::Format, "%msg");
      defaultConf.setGlobally(el::ConfigurationType::ToFile, "false");
      defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
      defaultConf.setGlobally(el::ConfigurationType::MillisecondsWidth, PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH);

      defaultConf.set(el::Level::Info, el::ConfigurationType::Format,
                      TIMESTAMP + OUT::blue + LEVEL + " " + MESSAGE + OUT::reset);

      defaultConf.set(el::Level::Debug, el::ConfigurationType::Format,
//                       TIMESTAMP + LEVEL + " " + POSITION + " " + MESSAGE + OUT::reset);
                      TIMESTAMP + OUT::white + LEVEL + " " + MESSAGE + OUT::reset);

      defaultConf.set(el::Level::Warning, el::ConfigurationType::Format,
                      TIMESTAMP + OUT::magenta + LEVEL + " " + MESSAGE + OUT::reset);

      defaultConf.set(el::Level::Error, el::ConfigurationType::Format,
                      TIMESTAMP + OUT::red + LEVEL + " " + MESSAGE + OUT::reset);

      defaultConf.set(el::Level::Fatal, el::ConfigurationType::Format,
                      TIMESTAMP + OUT::red + OUT::bold + LEVEL + " " + POSITION + " " + MESSAGE + OUT::reset);

      defaultConf.set(el::Level::Verbose, el::ConfigurationType::Format,
                      TIMESTAMP + VLEVEL + " " + MESSAGE + OUT::reset);

      el::Loggers::reconfigureAllLoggers(defaultConf);
    }

    /**
     * sets some default flags for easylogging++
     *
     * Current defaults are:
     *
     * - NewLineForContainer
     * - LogDetailedCrashReason
     * - DisableApplicationAbortOnFatalLog
     * - ColoredTerminalOutput
     * - MultiLoggerSupport
     * - AutoSpacing
     *
     * \see https://github.com/easylogging/easyloggingpp#logging-flags
     */
    inline static void set_logging_flags()
    {
//       el::Loggers::addFlag(el::LoggingFlag::NewLineForContainer);
      el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
      el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
      el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
      el::Loggers::addFlag(el::LoggingFlag::MultiLoggerSupport);
      el::Loggers::addFlag(el::LoggingFlag::CreateLoggerAutomatically);
      el::Loggers::addFlag(el::LoggingFlag::AutoSpacing);
    }

#ifdef NDEBUG
    inline static void test_logging_levels() {}
#else
    inline static void test_logging_levels()
    {
      cout << "### Example of different Logging Levels:" << endl;
      LOG(INFO) << "info";
      LOG(DEBUG) << "debug";
      LOG(WARNING) << "warning";
      LOG(ERROR) << "error";
      LOG(FATAL) << "fatal error";
      LOG(TRACE) << "trace";
      for (size_t level = 0; level <= 9; ++level) {
        VLOG(level) << "verbose level" << level;
      }
      cout << "### End Example Logging Levels" << endl << endl;
    }
#endif

    /**
     * starts easylogging++ with given arguments and loads configuration
     *
     * Usually, you want to pass the command line arguments given to `main()` in here and let easylogging++ figure
     * out what it needs.
     *
     * \param[in] argc number of command line arguments
     * \param[in] argv command line arguments
     */
    inline static void start_log(int argc, char** argv)
    {
      _START_EASYLOGGINGPP(argc, argv);
      set_logging_flags();
      load_default_config();
      pfasst::log::stack_position = 0;
    }
  }  // ::pfasst::log
}  // ::pfasst

#endif  // _PFASST__LOGGING_HPP_
