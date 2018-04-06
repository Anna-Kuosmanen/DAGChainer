# Uncomment the next two lines to use parallel processing (OpenMP)
#PARALLEL_FLAGS = -DPARALLEL_SUPPORT -fopenmp
#PARALLEL_LIB = -lgomp

CC=/usr/bin/g++
LEMONPATH=lemon/
LIBLEMON=$(LEMONPATH)lib/libemon.a
#BAMTOOLSPATH=bamtools/
#LIBBAMTOOLS=$(BAMTOOLSPATH)lib/libbamtools.a
LIBFLAGS= #-lbamtools 
CPPFLAGS = -Wall -g -O3 -DNDEBUG -std=c++0x $(PARALLEL_FLAGS) -I $(LEMONPATH)include/ -L $(LEMONPATH)lib/ #-I $(BAMTOOLSPATH)include/ -L $(BAMTOOLSPATH)lib/
OBJS= ../../SamAlignment.cpp ../../SamReader.cpp

all: traphlor runTraphlor

traphlor: $(OBJS) traphlor.o
	$(CC) $(CPPFLAGS) -o traphlor traphlor.cpp $(LIBFLAGS) 

runTraphlor: $(OBJS) runTraphlor.o
	$(CC) $(CPPFLAGS) -o ../runTraphlor runTraphlor.cpp $(OBJS) $(LIBFLAGS) $(PARALLEL_LIB)

.SUFFIXES:
.SUFFIXES: .o .cpp
.cpp.o: ; $(CC) $(CPPFLAGS) -MMD -c $*.cpp -o $@ $(INCLUDES) $(LIBRARY)



clean:
	rm -f traphlor ../runTraphlor







