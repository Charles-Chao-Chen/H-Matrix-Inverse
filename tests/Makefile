CC = g++

EIGEN   = /usr/include/eigen3/
INCLUDE = -I ../include -I $(EIGEN)

CFLAGS  = -g -Wall -O3 $(INCLUDE) -DDEBUG
HEADER  = ../include/hmap.hpp

EXE = solve

OBJS =  solve.o \
	../src/hmat.o ../src/node.o ../src/vec2.o

%.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(EXE) : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

clean :
	rm $(EXE) $(OBJS) *~ ../src/*~ ../include/*~

