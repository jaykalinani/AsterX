CXX=g++
CXXFLAGS=
LIBS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5
OBJ=Main.o RIT_EOS_standalone.o
DEPS=RIT_EOS_standalone.hh

Main: $(OBJ)
	$(CXX) $(OBJ) $(LIBS) -o Main
