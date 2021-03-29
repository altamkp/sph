CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
HDRS = sph.h
OBJS = cwMain.o sph.o sphInit.o sphPara.o sphDisp.o
LIBS = -lboost_program_options -lmpi_cxx -lmpi -llapack -lblas

%.o : %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< -I/usr/lib/x86_64-linux-gnu/openmpi/include

sph: $(OBJS)
	$(CXX) -o $@ $^ $(LIBS)

all: sph

.PHONY: clean
	target

clean:
	-rm -f *.o sph
