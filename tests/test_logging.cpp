#include "fixtures/test_helpers.hpp"

#include <pfasst/logging.hpp>
using namespace pfasst::log;


TEST(Formatting, format_mpi_rank) {
  int rank = 0;
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0) {
    EXPECT_THAT(format_mpi_rank(), StrEq("   0"));
  } else if (rank == 1) {
    EXPECT_THAT(format_mpi_rank(), StrEq("   1"));
  } else if (rank == 10) {
    EXPECT_THAT(format_mpi_rank(), StrEq("  10"));
  }
}

TEST(Formatting, log_file_name) {
  int rank = 0;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == 0) {
#ifdef WITH_MPI
    EXPECT_THAT(get_log_file_name(), StrEq("mpi-rank-0000.log"));
#else
    EXPECT_THAT(get_log_file_name(), StrEq(".log"));
#endif

#ifdef WITH_MPI
    EXPECT_THAT(get_log_file_name("prefix"), StrEq("prefix_mpi-rank-0000.log"));
#else
    EXPECT_THAT(get_log_file_name("prefix"), StrEq("prefix.log"));
#endif
  } else if (rank == 1) {
    EXPECT_THAT(get_log_file_name(), StrEq("mpi-rank-0001.log"));
    EXPECT_THAT(get_log_file_name("prefix"), StrEq("prefix_mpi-rank-0001.log"));
  } else if (rank == 10) {
    EXPECT_THAT(get_log_file_name(), StrEq("mpi-rank-0010.log"));
    EXPECT_THAT(get_log_file_name("prefix"), StrEq("prefix_mpi-rank-0010.log"));
  }
}


TEST(Colourizing, non_bold_colours)
{
  LOG(INFO) << OUT::reset << OUT::black << "black" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::red << "red" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::green << "green" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::yellow << "yellow" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::blue << "blue" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::magenta << "magenta" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::cyan << "cyan" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::white << "white" << OUT::reset;
}

TEST(Colourizing, bold_colours)
{
  LOG(INFO) << OUT::reset << OUT::bold << OUT::black << "bold black" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::bold << OUT::red << "bold red" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::bold << OUT::green << "bold green" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::bold << OUT::yellow << "bold yellow" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::bold << OUT::blue << "bold blue" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::bold << OUT::magenta << "bold magenta" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::bold << OUT::cyan << "bold cyan" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::bold << OUT::white << "bold white" << OUT::reset;
}

TEST(Colourizing, underline_colours)
{
  LOG(INFO) << OUT::reset << OUT::underline << OUT::black << "underlined black" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::underline << OUT::red << "underlined red" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::underline << OUT::green << "underlined green" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::underline << OUT::yellow << "underlined yellow" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::underline << OUT::blue << "underlined blue" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::underline << OUT::magenta << "underlined magenta" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::underline << OUT::cyan << "underlined cyan" << OUT::reset;
  LOG(INFO) << OUT::reset << OUT::underline << OUT::white << "underlined white" << OUT::reset;
}

TEST(Colourizing, reset_formatting)
{
  LOG(INFO) << OUT::reset << OUT::red << "red" << OUT::reset << " reset " << OUT::blue << "blue" << OUT::reset;
}


TEST_MAIN()
