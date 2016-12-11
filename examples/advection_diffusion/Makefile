PFASST   := $(abspath ../..)
CXX      := mpic++
CXXFLAGS ?= -std=c++11 -g -I$(PFASST)/include -I/usr/include/eigen3
LDFLAGS  ?= -lfftw3

all: vanilla_sdc mpi_pfasst serial_mlsdc

vanilla_sdc: vanilla_sdc.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

serial_mlsdc: serial_mlsdc.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)

mpi_pfasst: mpi_pfasst.cpp
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)
