PFASST   := $(abspath ../..)
CXX      := mpic++
CXXFLAGS ?= -std=c++11 -g -I$(PFASST)/include -I/usr/include/eigen3

all: heat

heat: heat.cpp ex-02-lib.o
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS)
