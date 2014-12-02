PFASST   := $(abspath ../..)
CXX      := mpic++
STACK    ?= $(PFASST)/stack
CXXFLAGS ?= -std=c++11 -I$(PFASST)/include -I$(STACK)/include -I$(STACK)/include/eigen3

LDFLAGS  ?= -L$(STACK)/lib -Wl,-rpath=$(STACK)/lib -lboost_program_options -lfftw3

all: vanilla_sdc mpi_pfasst serial_mlsdc

vanilla_sdc: vanilla_sdc.cpp
	g++ -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

serial_mlsdc: serial_mlsdc.cpp
	g++ -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

mpi_pfasst: mpi_pfasst.cpp
	mpic++ -o $@ $(CXXFLAGS) $^ $(LDFLAGS)
