/**
 * @file pfasst/logging.hpp
 * @since v0.3.0
 * 
 * @section logging_macros Logging Macros
 *
 * The most important and most frequently used logging macros provided by easylogging++ are wrapped
 * in the `ML_*` macros to avoid the evaluation of string conversion operations when logging is
 * disabled.
 *
 * Without these wrappers, using the easylogging++ logging macros will lead to evaluation of string
 * operations even with logging disabled at compile time.
 *
 * Instead of using
 *
 *     LOG(INFO) << "This is a logging line with string conversion " << my_variable;
 *
 * one should now use
 *
 *     ML_LOG(INFO, "This is a logging line with string conversion " << my_variable);
 */
#ifndef _PFASST__LOGGING_HPP_
#define _PFASST__LOGGING_HPP_

#include <iomanip>
#ifndef NDEBUG
  #include <iostream>  // used by pfasst::log::test_logging_levels()
#endif
#include <string>
using namespace std;

#ifdef WITH_MPI
  #include <sstream>
  #include <mpi.h>
#endif

#include <boost/algorithm/string.hpp>

#include "pfasst/config.hpp"


/**
 * Some general variables for coloured terminal output.
 *
 * @since v0.3.0
 */
struct OUT
{
  public:
    //! @{
    static const string black;
    static const string red;
    static const string green;
    static const string yellow;
    static const string blue;
    static const string magenta;
    static const string cyan;
    static const string white;
    //! @}

    //! @{
    static const string bold;
    static const string underline;
    //! @}

    //! resets both, colour and style
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
  #define ML_NOLOG
#endif

#include <pfasst/easylogging++.h>

//! @{
#ifndef ML_NOLOG
  /**
   * same as `LOG(level, x)` from easylogging++
   * 
   * @see @ref logging_macros
   */
  #define ML_LOG(level, x) LOG(level) << x

  /**
   * same as `CLOG(level, logger, x)` from easylogging++
   *
   * @see @ref logging_macros
   */
  #define ML_CLOG(level, logger_id, x) CLOG(level, logger_id) << x

  /**
   * same as `CLOG_IF(condition, level, logger, x)` from easylogging++
   *
   * @see @ref logging_macros
   */
  #define ML_CLOG_IF(condition, level, logger_id, x) CLOG_IF(condition, level, logger_id) <<  x

  /**
   * same as `CVLOG(verbosity, logger, x)` from easylogging++
   *
   * @see @ref logging_macros
   */
  #define ML_CVLOG(verbose_level, logger_id, x) CVLOG(verbose_level, logger_id) << x

  /**
   * same as `CVLOG_IF(condition, verbosity, logger, x)` from easylogging++
   *
   * @see @ref logging_macros
   */
  #define ML_CVLOG_IF(condition, verbose_level, logger_id, x) CVLOG_IF(condition, verbose_level, logger_id) << x
#else
  #define ML_LOG(level, x)
  #define ML_CLOG(level, logger_id, x)
  #define ML_CLOG_IF(condition, level, logger_id, x)
  #define ML_CVLOG(verbose_level, logger_id, x)
  #define ML_CVLOG_IF(condition, verbose_level, logger_id, x)
#endif
//! @}


#ifndef PFASST_LOGGER_INITIALIZED
  // initialize easyloggingpp
  // FIXME: this might already be called by code using PFASST++
  INITIALIZE_EASYLOGGINGPP
  /**
   * Guard symbol to ensure easylogging++ is only initialized once.
   *
   * When this symbol is defined, it is expected that `_INITIALIZE_EASYLOGGINGPP` has been called
   * once to initialize easylogging++.
   *
   * @note In case the executable using PFASST++ is also using easylogging++ as its logging library
   *   and initializing it prior to including the PFASST++ headers, please define this symbol prior
   *   including any of the PFASST++ headers.
   */
  #define PFASST_LOGGER_INITIALIZED
#endif

#ifndef PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH
  //! precision of milliseconds to be printed
  #define PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH "2"
#endif

#ifndef VLOG_FUNC_START
  /**
   * @deprecated The procedure to use verbose logging this way should be removed in favour of using
   *   `CVLOG(level, logger)` directly.
   */
  #define VLOG_FUNC_START(scope) \
    pfasst::log::stack_position++; \
    VLOG(9) << std::string((pfasst::log::stack_position - 1) * 2, ' ') << "START: " << std::string(scope) << "::" << std::string(__func__) << "()"
#endif

#ifndef VLOG_FUNC_END
  /**
   * @deprecated The procedure to use verbose logging this way should be removed in favour of using
   *   `CVLOG(level, logger)` directly.
   */
  #define VLOG_FUNC_END(scope) \
    pfasst::log::stack_position--; \
    VLOG(9) << std::string(pfasst::log::stack_position * 2, ' ') << "DONE:  " << std::string(scope) << "::" << std::string(__func__) << "()";
#endif

#ifndef LOG_PRECISION
  #define LOG_PRECISION 5
#endif

/**
 * Utility macro for creating identation depending on current stack position.
 *
 * This uses pfasst::log::stack_position to create a spacing string twice as long.
 */
#define LOG_INDENT \
  string(pfasst::log::stack_position * 2, ' ')

/**
 * Length of logger ID to print.
 *
 * Longer logger IDs will usually get cut off.
 *
 * @see pfasst::log::add_custom_logger()
 */
#define LOGGER_ID_LENGTH 10


namespace pfasst
{
  /**
   * Logging facilities for PFASST++.
   *
   * As PFASST++ is using easylogging++ as the logging backend, there are six distinced logging
   * levels plus ten additional verbose logging levels.
   *
   * To achieve consistent logging throughout PFASST++ and its examples, we agree on the following
   * conventions for using the different levels:
   *
   *   - `INFO` is used for general important messages to the user
   *   - `DEBUG` is ment for developping purposes and is only active when compiled without `-DNDEBUG`
   *   - `VLOG` - the verbose logging levels are used as follows:
   *     - 0 to 8
   *     - 9 for function enter and exit messages (cfg. @ref VLOG_FUNC_START and @ref VLOG_FUNC_END)
   *
   * @see [easylogging++](https://github.com/easylogging/easyloggingpp)
   */
  namespace log
  {
    /**
     * @details Global static helper variable to achieve nested indentation of logging output across
     *   different scopes.
     *
     * @since v0.3.0
     */
    static size_t stack_position;

    /**
     * internal flag identifying whether the default logger has been initialized.
     *
     * @since v1.0.0
     */
    static bool initialized = false;

    /**
     * Formats the local MPI world rank as a string.
     *
     * The local rank number as given by pfasst::config::get_rank() is filled from the left with the
     * @p fill character.
     *
     * @param[in] width width of the formatted rank
     * @param[in] fill  character used for left-filling
     * @returns formatted rank
     *
     * @since v1.0.0
     *
     * @ingroup Internals
     */
    inline string format_mpi_rank(const char fill = ' ')
    {
      ostringstream frmter;
      frmter << std::setw(4) << std::setfill(fill) << pfasst::config::get_rank();
      return frmter.str();
    }

    /**
     * Composes the log file name for the local rank.
     *
     * The rank-local log file name is composed of a potentially given _log_prefix_ as given from
     * the command line.
     * In case MPI is enabled this prefix gets extended by `_mpi-rank-<RANK>` where `<RANK>` is
     * the formatted MPI rank as determined by format_mpi_rank().
     * The file extension is always `.log`.
     *
     * @returns rank-local log file name
     *
     * @since v1.0.0
     *
     * @ingroup Internals
     */
    inline string get_log_file_name()
    {
      string log_name = config::get_value<string>("log_prefix", "");
#ifdef WITH_MPI
      if (log_name.size() > 0) {
        log_name += "_";
      }
      log_name += "mpi-rank-" + format_mpi_rank('0');
#endif
      log_name += ".log";
      return log_name;
    }

    /**
     * Sets global logging options.
     *
     * @param[in,out] conf         configuration to extend and set
     * @param[in]     default_conf optional, default configuration set to be used to retreive
     *   default options
     *
     * @since v1.0.0
     *
     * @ingroup Internals
     */
    inline void set_global_logging_options(el::Configurations* conf,
                                           const el::Configurations* default_conf = nullptr)
    {
      string milliseconds_width, to_stdout;
      if (default_conf) {
        el::Configurations* default_conf_nc = const_cast<el::Configurations*>(default_conf);
        milliseconds_width = default_conf_nc->get(el::Level::Info,
                                                  el::ConfigurationType::MillisecondsWidth)->value();
        to_stdout = default_conf_nc->get(el::Level::Info,
                                         el::ConfigurationType::ToStandardOutput)->value();
      } else {
        milliseconds_width = PFASST_LOGGER_DEFAULT_GLOBAL_MILLISECOND_WIDTH;
        to_stdout = (pfasst::config::options::get_instance().get_variables_map()
                                                            .count("quiet")) ? "false" : "true";
      }

      conf->setGlobally(el::ConfigurationType::ToStandardOutput, to_stdout);
      conf->setGlobally(el::ConfigurationType::Filename, get_log_file_name());
    }

    /**
     * Provides convenient way of adding additional named loggers.
     *
     * With this function one can easily create additional named loggers distinctable by the `id`.
     * The first @ref LOGGER_ID_LENGTH characters of the ID (as uppercase) will be included in
     * every line of the log.
     * The ID is used in the actual logging calls.
     *
     * @code{.cpp}
     * add_custom_logger("MyCustomLogger")
     * // somewhere else in the code
     * CLOG(INFO, "MyCustomLogger") << "a logging message";
     * @endcode
     *
     * This results in the log line (for the default value of @ref LOGGER_ID_LENGTH):
     *
     *     <TIME> [MYCUSTOM, INFO ] a logging message
     *
     * In case MPI is enabled, the line will also contain the MPI rank as given by
     * pfasst::config::get_rank():
     *
     *     <TIME> [MYCUSTOM, INFO , MPI  <rank> ] a logging message
     *
     * @param[in] id the ID of the logger; this is used in logging calls
     *
     * @note Please make sure to use `CLOG` (and `CVLOG` for verbose logging) to be able to specify
     *   a specific logger.
     *   Otherwise the default logger will be used.
     *
     * @since v0.3.0
     */
    inline static void add_custom_logger(const string& id)
    {
      bool colorize = pfasst::config::options::get_instance().get_variables_map()
                                                             .count("nocolor") ? false : true;

      const string INFO_COLOR = (colorize) ? OUT::blue : "";
      const string DEBG_COLOR = (colorize) ? "" : "";
      const string WARN_COLOR = (colorize) ? OUT::magenta : "";
      const string ERRO_COLOR = (colorize) ? OUT::red : "";
      const string FATA_COLOR = (colorize) ? OUT::red + OUT::bold : "";
      const string VERB_COLOR = (colorize) ? OUT::white : "";
      const string TIMESTAMP_COLOR = (colorize) ? OUT::white : "";
      const string RESET = (colorize) ? OUT::reset : "";

      const string TIMESTAMP = TIMESTAMP_COLOR + "%datetime{%H:%m:%s,%g}" + RESET + " ";
      const string LEVEL = "%level";
      const string VLEVEL = "VERB%vlevel";
      const string POSITION = "%fbase:%line";
      const string MESSAGE = "%msg";
#ifdef WITH_MPI
      const string MPI_RANK = ", MPI " + format_mpi_rank();
#else
      const string MPI_RANK = "";
#endif

      const size_t id_length = id.size();
      string id2print = id.substr(0, LOGGER_ID_LENGTH);
      boost::to_upper(id2print);
      if (initialized) {
        ML_CLOG(DEBUG, "default", "initializing custom logger '" << id << "' as '" << id2print << "'");
      }
      if (id_length < LOGGER_ID_LENGTH) {
        id2print.append(LOGGER_ID_LENGTH - id_length, ' ');
      }

      el::Logger* logger = el::Loggers::getLogger(id);
      el::Configurations* conf = logger->configurations();
      const el::Configurations* default_conf = el::Loggers::defaultConfigurations();
      set_global_logging_options(conf, default_conf);

      conf->set(el::Level::Info, el::ConfigurationType::Format,
                TIMESTAMP + INFO_COLOR + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + MESSAGE + RESET);
      conf->set(el::Level::Debug, el::ConfigurationType::Format,
                TIMESTAMP + DEBG_COLOR + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + POSITION + " " + MESSAGE + RESET);
      conf->set(el::Level::Warning, el::ConfigurationType::Format,
                TIMESTAMP + WARN_COLOR + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + MESSAGE + RESET);
      conf->set(el::Level::Error, el::ConfigurationType::Format,
                TIMESTAMP + ERRO_COLOR + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + MESSAGE + RESET);
      conf->set(el::Level::Fatal, el::ConfigurationType::Format,
                TIMESTAMP + FATA_COLOR + "[" + id2print + ", " + LEVEL  + MPI_RANK + "] " + POSITION + " " + MESSAGE + RESET);
      conf->set(el::Level::Verbose, el::ConfigurationType::Format,
                TIMESTAMP + VERB_COLOR + "[" + id2print + ", " + VLEVEL + MPI_RANK + "] " + MESSAGE + RESET);
      el::Loggers::reconfigureLogger(logger, *conf);
    }

    /**
     * Sets default configuration for default loggers.
     *
     * @since v0.3.0
     *
     * @ingroup Internals
     */
    inline static void load_default_config()
    {
      el::Configurations defaultConf;
      defaultConf.setToDefault();

      set_global_logging_options(&defaultConf);

      el::Loggers::setDefaultConfigurations(defaultConf, true);

      add_custom_logger("default");
      pfasst::log::initialized = true;
      add_custom_logger("Controller");
      add_custom_logger("Communicator");
      add_custom_logger("Sweeper");
      add_custom_logger("Encap");
      add_custom_logger("Quadrature");
      add_custom_logger("User");
    }

    /**
     * Sets some default flags for easylogging++.
     *
     * Current defaults are:
     *
     * - LogDetailedCrashReason
     * - DisableApplicationAbortOnFatalLog
     * - ColoredTerminalOutput
     * - MultiLoggerSupport
     * - CreateLoggerAutomatically
     *
     * @see
     *   [easylogging++ documentation on logging flags]
     *   (https://github.com/easylogging/easyloggingpp#logging-flags)
     *
     * @since v0.3.0
     *
     * @ingroup Internals
     */
    inline static void set_logging_flags()
    {
      el::Loggers::addFlag(el::LoggingFlag::LogDetailedCrashReason);
      el::Loggers::addFlag(el::LoggingFlag::DisableApplicationAbortOnFatalLog);
#ifdef NO_COLOR
      el::Loggers::removeFlag(el::LoggingFlag::ColoredTerminalOutput);
#else
      el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
#endif
      el::Loggers::addFlag(el::LoggingFlag::MultiLoggerSupport);
      el::Loggers::addFlag(el::LoggingFlag::CreateLoggerAutomatically);
    }

    /**
     * @fn inline static void test_logging_levels()
     *
     * Prints example logging for all logging levels for the default logger.
     *
     * @note When compiling with `-DNDEBUG` this function is a no-op.
     *
     * @since v0.3.0
     *
     * @ingroup Internals
     */
#ifdef NDEBUG
    inline static void test_logging_levels() {}
#else
    inline static void test_logging_levels()
    {
      cout << "### Example of different Logging Levels:" << endl;
      ML_LOG(INFO, "info");
      ML_LOG(DEBUG, "debug");
      ML_LOG(WARNING, "warning");
      ML_LOG(ERROR, "error");
      ML_LOG(FATAL, "fatal error");
      ML_LOG(TRACE, "trace");
      for (size_t level = 0; level <= 9; ++level) {
        ML_CVLOG(level, "default", "verbosity level " << level);
      }
      cout << "### End Example Logging Levels" << endl << endl;
    }
#endif

    /**
     * Starts easylogging++ with given arguments and loads configuration.
     *
     * Usually, you want to pass the command line arguments given to `%main(int, char**)` in here
     * and let easylogging++ figure out what it needs.
     *
     * @param[in] argc number of command line arguments as provided by `%main(int, char**)`
     * @param[in] argv command line arguments as provided by `%main(int, char**)`
     *
     * @note This function should not be called directly from user code but through pfasst::init.
     *
     * @since v0.3.0
     *
     * @ingroup Internals
     */
    inline static void start_log(int argc, char** argv)
    {
      START_EASYLOGGINGPP(argc, argv);
      set_logging_flags();
      load_default_config();
      pfasst::log::stack_position = 0;
      ML_CLOG(INFO, "default", "PFASST++ version " << pfasst::VERSION);
    }
  }  // ::pfasst::log
}  // ::pfasst

#endif  // _PFASST__LOGGING_HPP_
