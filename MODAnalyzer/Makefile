CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11 -g -Wno-deprecated-declarations

FASTJETINC = $(shell fastjet-config --cxxflags)
PYTHIA8INC = $(shell pythia8-config --cxxflags)

FASTJETLIB = $(shell fastjet-config --libs)
PYTHIA8LIB = $(shell pythia8-config --libs)
LIBRARIES += $(PYTHIA8LIB) $(FASTJETLIB) -lgmp -lCGAL
INCLUDE   += $(PYTHIA8INC) $(FASTJETINC) -Iinclude

$(shell mkdir -p tmp bin)

.PHONY: clean

all:
	make dijets_analyzer

dijets_analyzer: tmp/Main.o tmp/CmdLine.o tmp/uMODAnalyzer.o tmp/DijetsAnalyzer.o tmp/JetMatcher.o tmp/Utils.o
	$(CXX) -o bin/$@ $^ $(CXXFLAGS) $(LIBRARIES)

tmp/Main.o: src/Main.cc
	$(CXX) $(INCLUDE) -o $@ $< -c $(CXXFLAGS)

tmp/CmdLine.o: src/CmdLine.cc
	$(CXX) $(INCLUDE) -o $@ $< -c $(CXXFLAGS)

tmp/uMODAnalyzer.o: src/uMODAnalyzer.cc
	$(CXX) $(INCLUDE) -o $@ $< -c $(CXXFLAGS)

tmp/DijetsAnalyzer.o: src/DijetsAnalyzer.cc
	$(CXX) $(INCLUDE) -o $@ $< -c $(CXXFLAGS)

tmp/JetMatcher.o: src/JetMatcher.cc
	$(CXX) $(INCLUDE) -o $@ $< -c $(CXXFLAGS)

tmp/Utils.o: src/Utils.cc
	$(CXX) $(INCLUDE) -o $@ $< -c $(CXXFLAGS)

clean:
	rm -rfv bin tmp
