CXXFLAGS= -std=c++0x

all: rpa dens_test

rpa: *.cc *.h 
	g++ rpa_lib.cc -DRPA_MAIN -o rpa

densities.o: *.h ../*.h *.cc

dens_test: ../generatormt.o ../nucleus_data.o ../elements.o

