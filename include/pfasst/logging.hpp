#ifndef _PFASST__LOGGING_HPP_
#define _PFASST__LOGGING_HPP_

#include "pfasst/site_config.hpp"

#include <string>
using namespace std;

#ifdef WITH_MPI
  #include <mpi.h>
#endif

#include <boost/algorithm/string.hpp>

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
#define ELPP_STL_LOGGING
// disable creation of default log file
#define ELPP_NO_DEFAULT_LOG_FILE
// enable passing `--logging-flags` via command line
#define ELPP_LOGGING_FLAGS_FROM_ARG

#ifndef NDEBUG
  #define ELPP_DEBUG_ASSERT_FAILURE
  #define ELPP_STACKTRACE_ON_CRASH
#endif

#ifdef PFASST_NO_LOGGING
  #define ELPP_DISABLE_LOGS
#endif

#include <pfasst/easylogging++.h>


#ifndef PFASST_LOGGER_INITIALIZED
  // initialize easyloggingpp
  // FIXME: this might already be called by code using PFASST++
  INITIALIZE_EASYLOGGINGPP
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
    VLOG(9) << std::string((pfasst::log::stack_position - 1) * 2, ' ') << "START: " << std::string(scope) << "::" << std::string(__func__) << "()"
#endif

#ifndef VLOG_FUNC_END
  #define VLOG_FUNC_END(scope) \
    pfasst::log::stack_position--; \
    VLOG(9) << std::string(pfasst::log::stack_position * 2, ' ') << "DONE:  " << std::string(scope) << "::" << std::string(__func__) << "()";
#endif

#ifndef LOG_PRECISION
  #define LOG_PRECISION 5
#endif

#define LOG_INDENT string(pfasst::log::stack_position * 2, ' ')

//! length of logger ID to print
#define LOGGER_ID_LENGTH 8

namespace pfasst
{
  /**
   * Logging Facilities for PFASST++
   *
   * As PFASST++ is using easylogging++ as the logging backend, there are six distinced logging levels plus ten
   * additional verbose logging levels.
   *
   * To achieve consistent logging throughout PFASST++ and its examples, we agree on the following conventions for using
   * the different levels:
   *
   *   - `INFO` is used for general important messages to the user
   *   - `DEBUG` is ment for developping purposes and is only active when compiled without `-DNDEBUG`
   *   - `VLOG` - the verbose logging levels are used as follows:
   *     - 0 to 8
   *     - 9 for function enter and exit messages (\see VLOG_FUNC_START and \see VLOG_FUNC_END )
   */
  namespace log
  {
    static size_t stack_position;

    /**
     * \brief provides convenient way of adding additional named loggers
     * \details With this function one can easily create additional named loggers distinctable by the `id`.
     *   The first \ref LOGGER_ID_LENGTH characters of the ID (as uppercase) will be included in every line of the log.
     *   The ID is used in the actual logging calls.
     *
     *   \code{.cpp}
     *       add_custom_logger("MyCustomLogger")
     *       // somewhere else in the code
     *       CLOG(INFO, "MyCustomLogger") << "a logging message";
     *   \endcode
     *
     *   This results in the log line (for the default value of \ref LOGGER_ID_LENGTH):
     *
     *       <TIME> [MYCUSTOM, INFO ] a logging message
     * \note Please make sure to use `CLOG` (and `CVLOG` for verbose logging) to be able to specify a specific logger.
     *   Otherwise the default logger will be used.
     * \param[in] id The ID of the logger. This is used in logging calls.
     */
    inline static void add_custom_logger(const string& id)
    {
      const string TIMESTAMP = OUT::white + "%datetime{%H:%m:%s,%g}" + OUT::reset + " ";
      const string LEVEL = "%level]";
      const string VLEVEL = "VERB%vlevel]";
      const string POSITION = "%fbase:%line";
      const string MESSAGE = "%msg";

      const string INFO_COLOR = OUT::blue;
      const string DEBG_COLOR = "";
      const string WARN_COLOR = OUT::magenta;
      const string ERRO_COLOR = OUT::red;
      const string FATA_COLOR = OUT::red + OUT::bold;
      const string VERB_COLOR = OUT::white;

      const size_t id_length = id.size();
      string id2print = id.substr(0, LOGGER_ID_LENGTH);
      boost::to_upper(id2print);
      if (id_length < LOGGER_ID_LENGTH) {
        id2print.append(LOGGER_ID_LENGTH - id_length, ' ');
      }

      el::Logger* logger = el::Loggers::getLogger(id);
      el::Configurations* conf = logger->configurations();
      conf->setGlobally(el::ConfigurationType::MillisecondsWidth,
                        PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH);
#ifdef WITH_MPI
      int initialized = 0;
      MPI_Initialized(&initialized);
      assert((bool)initialized);
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      conf->setGlobally(el::ConfigurationType::ToFile, "true");
      conf->setGlobally(el::ConfigurationType::ToStandardOutput, "false");
      conf->setGlobally(el::ConfigurationType::Filename,
                        string("mpi_run_") + to_string(rank) + string(".log"));
#else
      conf->setGlobally(el::ConfigurationType::ToFile, "false");
      conf->setGlobally(el::ConfigurationType::ToStandardOutput, "true");
#endif
      conf->set(el::Level::Info, el::ConfigurationType::Format,
                TIMESTAMP + INFO_COLOR + "[" + id2print + ", " + LEVEL  + " " + MESSAGE + OUT::reset);
      conf->set(el::Level::Debug, el::ConfigurationType::Format,
                TIMESTAMP + DEBG_COLOR + "[" + id2print + ", " + LEVEL  + " " + POSITION + " " + MESSAGE + OUT::reset);
      conf->set(el::Level::Warning, el::ConfigurationType::Format,
                TIMESTAMP + WARN_COLOR + "[" + id2print + ", " + LEVEL  + " " + MESSAGE + OUT::reset);
      conf->set(el::Level::Error, el::ConfigurationType::Format,
                TIMESTAMP + ERRO_COLOR + "[" + id2print + ", " + LEVEL  + " " + MESSAGE + OUT::reset);
      conf->set(el::Level::Fatal, el::ConfigurationType::Format,
                TIMESTAMP + FATA_COLOR + "[" + id2print + ", " + LEVEL  + " " + POSITION + " " + MESSAGE + OUT::reset);
      conf->set(el::Level::Verbose, el::ConfigurationType::Format,
                TIMESTAMP + VERB_COLOR + "[" + id2print + ", " + VLEVEL + " " + MESSAGE + OUT::reset);
      el::Loggers::reconfigureLogger(logger, *conf);
    }

    /**
     * sets default configuration for default loggers
     */
    inline static void load_default_config()
    {
      el::Configurations defaultConf;
      defaultConf.setToDefault();

#ifdef WITH_MPI
      int initialized = 0;
      MPI_Initialized(&initialized);
      assert((bool)initialized);
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      defaultConf.setGlobally(el::ConfigurationType::ToFile, "true");
      defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "false");
      defaultConf.setGlobally(el::ConfigurationType::Filename,
                              string("mpi_run_") + to_string(rank) + string(".log"));
#else
      defaultConf.setGlobally(el::ConfigurationType::ToFile, "false");
      defaultConf.setGlobally(el::ConfigurationType::ToStandardOutput, "true");
#endif
      defaultConf.setGlobally(el::ConfigurationType::MillisecondsWidth, PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH);
      el::Loggers::reconfigureAllLoggers(defaultConf);

      add_custom_logger("default");
      add_custom_logger("Controller");
      add_custom_logger("Sweeper");
      add_custom_logger("Encap");
      add_custom_logger("Quadrature");
      add_custom_logger("User");
    }

    /**
     * sets some default flags for easylogging++
     *
     * Current defaults are:
     *
     * - LogDetailedCrashReason
     * - DisableApplicationAbortOnFatalLog
     * - ColoredTerminalOutput
     * - MultiLoggerSupport
     *
     * \see https://github.com/easylogging/easyloggingpp#logging-flags
     */
    inline static void set_logging_flags()
    {
      el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
      el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
      el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
      el::Loggers::addFlag(el::LoggingFlag::MultiLoggerSupport);
      el::Loggers::addFlag(el::LoggingFlag::CreateLoggerAutomatically);
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
        VLOG(level) << "verbosity level " << level;
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
      START_EASYLOGGINGPP(argc, argv);
      set_logging_flags();
      load_default_config();
      pfasst::log::stack_position = 0;
    }
  }  // ::pfasst::log
}  // ::pfasst

#endif  // _PFASST__LOGGING_HPP_
