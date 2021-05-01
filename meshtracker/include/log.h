#ifndef LOG_H
#define LOG_H

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>

namespace logging = boost::log;

static void init_logging
(
  const int         & verbosity, 
  const std::string &     fname
)
{
   logging::register_simple_formatter_factory<logging::trivial::severity_level, char>("Severity");


    logging::add_file_log(
      boost::log::keywords::file_name  = fname,
      boost::log::keywords::format     = "[%TimeStamp%] %Message%",
      boost::log::keywords::auto_flush = true // potential performance increase when false
    );

    // Severity Levels:
    // (0) trace,
    // (1) debug,
    // (2) info,
    // (3) warning,
    // (4) error,
    // (5) fatal

    logging::core::get()->set_filter
    (
        logging::trivial::severity >= logging::trivial::severity_level(verbosity)
    );

    logging::add_common_attributes();
}

#endif