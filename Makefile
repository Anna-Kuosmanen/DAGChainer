CPPFLAGS = -Wall -std=c++11
CPP= g++
INCLUDES= -I /cs/work/home/aekuosma/lemon/include
LIBRARY= -L /cs/work/home/aekuosma/lemon/lib

OBJS= MC-MPC/decomposer/MPC.o MC-MPC/decomposer/decomposition.o MC-MPC/util/utils.o SequenceGraph.o Vertex.o ColinearSolver.o RMaxQTree.o BruteForceSolver.o

all: pipeline

pipeline: $(OBJS) pipeline.o
	$(CPP) $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $(LIBRARY) -o pipeline pipeline.o $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .cpp
.cpp.o: ; $(CPP) $(CPPFLAGS) $(INCLUDES) $(LIBRARY) -MMD -c $*.cpp -o $@

clean:
	rm -f *.o *.d ../pipeline

