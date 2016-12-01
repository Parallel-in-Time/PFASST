PFASST   := $(abspath ../..)
EIGEN3   ?= $(PFASST)/include/eigen-eigen-b9cd8366d4e8
CXX      ?= g++
CXXFLAGS += -std=c++11 -g -I$(PFASST)/include -I$(EIGEN3)

quad: multiprecision.cpp
	$(CXX) -o $@ $(CXXFLAGS) $<
