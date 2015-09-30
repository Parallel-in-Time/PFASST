PFASST   := $(abspath ../..)
STACK    ?= $(PFASST)/stack
CXX      := $(STACK)/bin/mpic++
CXXFLAGS ?= -std=c++11 -I$(PFASST)/include -I$(PFASST)/src -I$(STACK)/include -I$(STACK)/include/eigen3

LDFLAGS  ?= -L$(STACK)/lib -Wl,-rpath=$(STACK)/lib -lboost_program_options -lfftw3

all: vanilla_sdc mpi_pfasst serial_mlsdc

vanilla_sdc: vanilla_sdc.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

serial_mlsdc: serial_mlsdc.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

mpi_pfasst: mpi_pfasst.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)
