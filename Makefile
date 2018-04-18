## Change these dirs to correspond to the locations of the libraries
LEMON_DIR=/cs/work/home/aekuosma/lemon/
SDSL_DIR=/cs/work/home/aekuosma/software/sdsl-lite
DIVSUFSORT_DIR=/cs/work/home/aekuosma/software
GCSA_DIR=/cs/work/home/aekuosma/software/gcsa2

CC= g++
CPPFLAGS = -Wall -std=c++11 -g -O3 -DNDEBUG -fopenmp -pthread
INCLUDES= -I$(LEMON_DIR)/include -I/$(GCSA_DIR)/include/ -I$(SDSL_DIR)/include/ -I$(DIVSUFSORT_DIR)/include
LIBRARY= -L$(LEMON_DIR)/lib -L$(GCSA_DIR)/lib -L$(SDSL_DIR)/lib -L$(DIVSUFSORT_DIR)/lib
LIBS= -lsdsl -ldivsufsort -ldivsufsort64
OBJS= SequenceGraph.o Vertex.o ColinearSolver.o RMaxQTree.o SamAlignment.o SamReader.o GenomeReader.o FastaReader.o MaximalExactMatch.o
MPCOBJS= MC-MPC/decomposer/MPC.o MC-MPC/decomposer/decomposition.o MC-MPC/util/utils.o
GCSAOBJS= $(GCSA_DIR)/gcsa.cpp $(GCSA_DIR)/lcp.cpp $(GCSA_DIR)/algorithms.o $(GCSA_DIR)/dbg.o $(GCSA_DIR)/files.o $(GCSA_DIR)/internal.o $(GCSA_DIR)/path_graph.o $(GCSA_DIR)/support.o $(GCSA_DIR)/utils.o

all: transcriptPipeline

transcriptPipeline: $(OBJS) $(MPCOBJS) $(GCSAOBJS) transcriptPipeline.o
	$(CC) $(CXXFLAGS) $(CPPFLAGS) -o transcriptPipeline transcriptPipeline.o $(OBJS) $(MPCOBJS) $(GCSAOBJS) $(INCLUDES) $(LIBRARY) $(LIBS)


.SUFFIXES:
.SUFFIXES: .o .cpp
.cpp.o: ; $(CC) $(CPPFLAGS) -MMD -c $*.cpp -o $@ $(INCLUDES) $(LIBRARY) $(LIBS)

clean:
	rm -f *.o *.d transcriptPipeline

